# -*- coding: utf-8 -*-
"""
Title: 
    ventilation_recon_3d.py

Agenda: 
    A .py script that uses compressed sensing methods to reconstruct 
non-cartesian MRI data.

Author: 
    Joseph Plummer - joseph.plummer@cchmc.org 

Creation date: 
    2023-10-30
    
Modification date: 
    --
    
"""

# %% Import packages
from tqdm import tqdm
import nibabel as nib
import cv2
from sklearn.cluster import KMeans
from functions import tv, optpoly
from functions import convexalg, util
import sigpy.mri as mr
import sigpy.plot as pl
import sigpy as sp
import numpy as np
import matplotlib.pyplot as plt
import random
import os
plt.style.use('dark_background')

# %% Reconstruction settings

# Device settings
devnum = 0  # -1 for CPU, 0 for GPU
device = sp.Device(devnum)
xp = device.xp
def mvd(x): return sp.to_device(x, device)
def mvc(x): return sp.to_device(x, sp.cpu_device)


# Density compensation
use_dcf = False

# Subsample
subsample_pc = 100  # percentage

# Regularization
lamda_l2 = 1e-5
lamda_1 = 1e-8
lamda_2 = 1e-8
num_iters = 50

# %% Import MRI k-space data/coordinates

# Filenames
data_file = "data/ventilation-3d/ksp.npy"
coords_file = "data/ventilation-3d/coord.npy"

# Load data (channels, projections, samples)
data = np.array(np.load(data_file))

# Visualize
fig, (ax1, ax2) = plt.subplots(2, sharex=False)
# First few points of kspace_data2 have zero signal due to LNA issue
ax1.plot(abs(data[0, :, 0]), color='r')
ax1.set_ylabel('$k_0$ intensity')
ax1.set_title("magnetization decay")
ax1.set_xlabel("projection number")
ax2.plot(abs(data[0, 0, :]), color='b')
ax2.set_ylabel('k-space intensity')
ax2.set_xlabel('sample number')
ax2.set_title("fid")
fig.tight_layout(pad=1.0, w_pad=0.5, h_pad=1.0)

# Load coords (projections, samples, dimensions)
coords = np.array(np.load(coords_file))

# Verify inputs
N_dimensions = np.shape(coords)[2]  # Number of dimensions
N_samp = np.shape(data)[2]  # Samples per projection
N_proj = np.shape(data)[1]  # Number of projections
N_channels = np.shape(data)[0]  # Number of channels

# Plot trajectories
plt.figure()
ax = plt.axes(projection='3d')
N_visual = N_proj//20  # Number of projections you want to show, inefficient for large n
color = iter(plt.cm.viridis(np.linspace(0, 1, N_visual)))
for i in np.linspace(0, N_proj-1, N_visual):
    i = int(i)
    c = next(color)
    ax.scatter(coords[i, :, 0],
               coords[i, :, 1],
               coords[i, :, 2],
               color=c, s=0.5, marker='.')
# ax.axis('off')
color_tuple = (0.0, 0.0, 0.0, 0.0)
ax.xaxis.set_pane_color(color_tuple)
ax.yaxis.set_pane_color(color_tuple)
ax.zaxis.set_pane_color(color_tuple)
ax.xaxis.line.set_color(color_tuple)
ax.yaxis.line.set_color(color_tuple)
ax.zaxis.line.set_color(color_tuple)
ax.set_zlabel('$k_z$', fontsize=18)
ax.set_ylabel('$k_y$', fontsize=18)
ax.set_xlabel('$k_x$', fontsize=18)
ax.set_yticklabels([])
ax.set_xticklabels([])
ax.set_zticklabels([])
plt.show()

# %% Reconstruction settings

# Image reconstruction sizes
recon_resolution = 120
scan_resolution = 120

# Make manipulations to reconstructed volume
alpha = 1  # Overgridding factor (alpha = 0.5 means reconstructing on 2x grid)
resize_factor = recon_resolution / scan_resolution
image_size = int(resize_factor * scan_resolution / alpha)
image_shape = (int(image_size), int(image_size), int(image_size))
coords_resize = image_size * coords / resize_factor

# %% Subset the data to experiment with undersampling

# Select only a subset of trajectories and data
N_proj_subset = int(N_proj * subsample_pc / 100)
subset = range(0, N_proj_subset)
data_subset = data[:, subset, :]
coords_subset = coords_resize[subset, ...]

# %% Calculate density compensation

# Pipe-Menon DCF settings
beta = 8  # Default = 8
width = 3  # Default = 4, choose ~3 for high SNR, choose 5~8 for sharpness

try:
    dcf = mr.pipe_menon_dcf(coords_subset, img_shape=np.dot(2, image_shape),
                            beta=beta, width=width, device=device)
    dcf /= np.max(dcf)
except:
    print("Error: DCF was not calculated. ")
    raise

# Visualize
fig = plt.figure()
plt.plot(sp.to_device(dcf[0, ...], -1), 'b')
plt.xlabel('Sample number')
plt.ylabel('DCF')
plt.suptitle('Sample Density Compensation Function')
plt.title(r'$ \beta $ = %.2g, width = %.2g' % (beta, width))

# Move to CPU
dcf = mvc(dcf)

# %% Generate a hyperpolarized weighting for the inverse NUFFT

# Simulate T2* decay
t2_star_xe = 18  # ms
readout = 10  # ms
dwell_time = readout/N_samp
relaxation = np.zeros((N_samp,))
for i in range(N_samp):
    relaxation[i] = np.exp(-i*dwell_time/t2_star_xe)

# Estimate magnetization decay due to RF excitations
k0 = abs(data_subset[0, :, 0])
k = np.zeros((N_proj_subset, N_samp))
for i in range(N_proj_subset):
    k[i, :] = k0[i] * relaxation

# Visualize
fig, (ax1, ax2) = plt.subplots(2, sharex=False)
# First few points of kspace_data2 have zero signal due to LNA issue
ax1.plot(abs(k[:, 0]), color='g')
ax1.set_ylabel('$signal intensity')
ax1.set_title("magnetization decay from excitation")
ax1.set_xlabel("projection number")
ax2.plot(abs(k[0, :]), color='y')
ax2.set_ylabel('signal intensity')
ax2.set_xlabel('sample number')
ax2.set_title("magnetization decay from T$_2^*$")
fig.tight_layout(pad=1.0, w_pad=0.5, h_pad=1.0)

# %% Estimate sensitivity map

# Use J-SENSE methods to calculate a sensitivity map
try:
    if N_channels != 1:
        # Works for multi-channel data
        sens_map = mr.app.JsenseRecon(data, coord=coords_resize, weights=dcf,
                                      mps_ker_width=16, ksp_calib_width=24, lamda=0, img_shape=image_shape, device=device).run()
        pl.ImagePlot(sens_map, z=0, mode='m',
                     title='Estimated sensitivity maps using J-Sense',
                     colormap="gray")
        sens_map = np.array(sens_map, np.complex64)
        print("Sensitivity map estimated using J-SENSE.")
    else:
        raise Exception(
            "N_channels = 1, single channel data used. J-SENSE not possible.")
except:
    # For single channel data, the sensitivity map can be assumed to be all ones
    sens_map = np.ones((1, image_size, image_size, image_size), dtype=int)
    print("Sensitivity map of all ones was assumed.")


# %% Inverse NUFFT reconstruction

# Move to GPU
data = mvd(data_subset)
coords = mvd(coords_subset)
dcf = mvd(dcf)
sens_map = mvd(sens_map)

# Compute linear operators
S = sp.linop.Multiply(image_shape, sens_map)
F = sp.linop.NUFFT(sens_map.shape,
                   coord=coords,
                   oversamp=1.25,  # Recommended 1-1.5
                   width=4,  # Recommended 2-4
                   toeplitz=True)
D = sp.linop.Multiply(F.oshape, dcf**0.5)

# Compute a single x = A.H b operation (i.e. inverse NUFFT)
K = sp.linop.Multiply(F.oshape, k**0)
A_dcf = K * D * F * S
b_dcf = data * dcf**0.5 / xp.linalg.norm(data * dcf**0.5)
img_nufft = mvc(abs(A_dcf.H * b_dcf))
pl.ImagePlot(np.rot90(img_nufft[:, :, int(recon_resolution/3):int(recon_resolution/1.5):2], 2),
             x=0, y=1, z=2,
             title="$(DF)^H(y)$")

# Compute a single x = A.H b operation (i.e. inverse NUFFT)
data_normalized = data
K = sp.linop.Multiply(F.oshape, k**0)
for i in range(N_proj_subset):
    data_normalized[:, i, :] = data[:, i, :] * np.mean(k0) / k0[i]
A_dcf = K * D * F * S
b_dcf = data_normalized * dcf**0.5 / xp.linalg.norm(data_normalized * dcf**0.5)
img_nufft_norm = mvc(abs(A_dcf.H * b_dcf))
pl.ImagePlot(np.rot90(img_nufft_norm[:, :, int(recon_resolution/3):int(recon_resolution/1.5):2], 2),
             x=0, y=1, z=2,
             title="$(DF)^H(y_n)$")

# Compute a single x = A.H b operation (i.e. inverse NUFFT)
K = sp.linop.Multiply(F.oshape, k**1)
A_dcf = K * D * F * S
b_dcf = data * dcf**0.5 / xp.linalg.norm(data * dcf**0.5)
img_nufft_model = mvc(abs(A_dcf.H * b_dcf))
pl.ImagePlot(np.rot90(img_nufft_model[:, :, int(recon_resolution/3):int(recon_resolution/1.5):2], 2),
             x=0, y=1, z=2,
             title="$(RDF)^H(y)$")


# %% Iterative NUFFT reconstruction

data = mvd(data)

# Declare if density compensated
K = sp.linop.Multiply(F.oshape, k**0)
if use_dcf:
    A = K * D * F * S
    data_norm = xp.linalg.norm(data * dcf**0.5)
    b = data * dcf**0.5 / data_norm
else:
    A = K * F * S
    data_norm = xp.linalg.norm(data)
    b = data / data_norm

# Calculate maximum eigenvalue
LL = sp.app.MaxEig(A.N, dtype=xp.complex64,
                   device=device).run() * 1.01

# Divide by svd eigenvalue to ensure maximum eigenvalue falls between [0,1]
A *= np.sqrt(1/LL)

# Solve ||Ax - b||^2 + lamda||x||^2 using conjugate gradient descent
img_cg = abs(mvc(convexalg.cg(num_iters=50,
                              ptol=1e-5,
                              A=A, b=b,
                              lamda=lamda_l2,  # l2-norm constraint for additional regularization
                              verbose=True,
                              draw_output=True,
                              save=None)))

pl.ImagePlot(np.rot90(img_cg[:, :, int(recon_resolution/3):int(recon_resolution/1.5):2], 2),
             x=0, y=1, z=2,
             title="$||Fx-y||^2_2$")

# Declare if density compensated
K = sp.linop.Multiply(F.oshape, k**0)
data_normalized = np.zeros_like(data_normalized)
for i in range(N_proj_subset):
    data_normalized[:, i, :] = data[:, i, :] * np.mean(k0) / k0[i]
if use_dcf:
    A = K * D * F * S
    data_norm = xp.linalg.norm(data_normalized * dcf**0.5)
    b = data_normalized * dcf**0.5 / data_norm
else:
    A = K * F * S
    data_norm = xp.linalg.norm(data_normalized)
    b = data_normalized / data_norm

# Calculate maximum eigenvalue
LL = sp.app.MaxEig(A.N, dtype=xp.complex64,
                   device=device).run() * 1.01

# Divide by svd eigenvalue to ensure maximum eigenvalue falls between [0,1]
A *= np.sqrt(1/LL)

# Solve ||Ax - b||^2 + lamda||x||^2 using conjugate gradient descent
img_cg_normalized = abs(mvc(convexalg.cg(num_iters=50,
                                         ptol=1e-5,
                                         A=A, b=b,
                                         lamda=lamda_l2,  # l2-norm constraint for additional regularization
                                         verbose=True,
                                         draw_output=True,
                                         save=None)))

pl.ImagePlot(np.rot90(img_cg_normalized[:, :, int(recon_resolution/3):int(recon_resolution/1.5):2], 2),
             x=0, y=1, z=2,
             title="$||Fx-y_n||^2_2$")

# Declare if density compensated
K = sp.linop.Multiply(F.oshape, k**1)
if use_dcf:
    A = K * D * F * S
    data_norm = xp.linalg.norm(data * dcf**0.5)
    b = data * dcf**0.5 / data_norm
else:
    A = K * F * S
    data_norm = xp.linalg.norm(data)
    b = data / data_norm

# Calculate maximum eigenvalue
LL = sp.app.MaxEig(A.N, dtype=xp.complex64,
                   device=device).run() * 1.01

# Divide by svd eigenvalue to ensure maximum eigenvalue falls between [0,1]
A *= np.sqrt(1/LL)

# Solve ||Ax - b||^2 + lamda||x||^2 using conjugate gradient descent
img_cg_model = abs(mvc(convexalg.cg(num_iters=50,
                                    ptol=1e-5,
                                    A=A, b=b,
                                    lamda=lamda_l2,  # l2-norm constraint for additional regularization
                                    verbose=True,
                                    draw_output=True,
                                    save=None)))

pl.ImagePlot(np.rot90(img_cg_model[:, :, int(recon_resolution/3):int(recon_resolution/1.5):2], 2),
             x=0, y=1, z=2,
             title=r"$||RFx-y||^2_2$")


# Define regularizing linear operators
W1 = sp.linop.Wavelet(S.ishape, wave_name="db4")
def g1(x): return lamda_1 * xp.linalg.norm(W1(x).ravel(), ord=1)


prox_g1 = sp.prox.UnitaryTransform(
    sp.prox.L1Reg(W1.oshape, lamda_1), W1)

prox_g2 = tv.ProxTV(S.ishape, lamda_2)
def g2(x): return lamda_2 * xp.linalg.norm(prox_g2.G(x))


lst_g = [g1, g2]
lst_proxg = [prox_g1, prox_g2]

# Solve ||Ax - b||^2 + lamda_1 |Wx| + lamda_2 |Gx| using ADMM
img_cs_cg = mvc(abs(convexalg.admm(num_iters=num_iters,
                                   ptol=1e-5,
                                   A=A, b=b,
                                   num_normal=15,
                                   lst_proxg=lst_proxg,
                                   rho=5e2,
                                   lst_g=lst_g,
                                   method="cg",
                                   verbose=True,
                                   draw_output=True)))
pl.ImagePlot(np.rot90(img_cs_cg[:, :, int(recon_resolution/3):int(recon_resolution/1.5):2], 2),
             x=0, y=1, z=2,
             title=r"$||RFx-y||^2_2 + g(x)$")


# %% Show the plots

plot_slice = image_size//2

plt.figure()
plt.imshow(np.rot90(img_nufft[1:image_size-1,
           1:image_size-1, plot_slice], 3), cmap="gray")
plt.colorbar()
plt.title('Reconstruction using inverse NUFFT \n (no norm)')
plt.axis('off')

plt.figure()
plt.imshow(np.rot90(
    img_cg_normalized[1:image_size-1, 1:image_size-1, plot_slice], 3), cmap="gray")
plt.colorbar()
plt.title('Reconstruction using iterative \n (norm decay)')
plt.axis('off')

plt.figure()
plt.imshow(np.rot90(img_cg_model[1:image_size-1,
           1:image_size-1, plot_slice], 3), cmap="gray")
plt.colorbar()
plt.title('Reconstruction using iterative NUFFT \n (model decay)')
plt.axis('off')

plt.figure()
plt.imshow(
    np.rot90(img_cs_cg[1:image_size-1, 1:image_size-1, plot_slice], 3), cmap="gray")
plt.colorbar()
plt.title('Reconstruction using GM iterative NUFFT \n (model decay)')
plt.axis('off')


# %% Try a signal mask using K means

# Choose an image
img = cv2.normalize(img_nufft[..., plot_slice],
                    None, 0, 1.0, cv2.NORM_MINMAX, dtype=cv2.CV_32F)

# Perform clusting
N_clusters = 2
kmeans = KMeans(n_clusters=N_clusters, random_state=0).fit(np.reshape(
    # Reshape to (voxelcount,1) shape array
    img, (np.shape(np.ravel(img))[0], 1)))
clustered = kmeans.cluster_centers_[kmeans.labels_]

# Binarize
if N_clusters == 2:
    clustered[clustered == np.percentile(clustered, 75)] = 1
    clustered[clustered == np.percentile(clustered, 25)] = 0

# Reshape clusters into images
img_cluster = np.zeros(shape=(N_clusters, np.shape(
    img)[0], np.shape(img)[1]))
labels = kmeans.labels_
for n in range(N_clusters):
    cluster_temp = []
    for i in range(len(labels)):
        if (labels[i]) == n:
            cluster_temp.append(float(clustered[i]))
        else:
            cluster_temp.append(1)
    img_cluster[n, ...] = np.array(
        cluster_temp).reshape(img.shape)

# Calculate norm of each cluster*img_full to see which one is larger (and which one is signal versus noise)
norm_cluster0 = np.linalg.norm(img_cluster[0, ...]*img)
norm_cluster1 = np.linalg.norm(img_cluster[1, ...]*img)

if norm_cluster0 >= norm_cluster1:
    img_signal = img_cluster[0, ...]
    img_noise = img_cluster[1, ...]
else:
    img_signal = img_cluster[1, ...]
    img_noise = img_cluster[0, ...]

# Binarize noise mask
# img_noise[img_noise != 1] = 0
img_noise = np.uint8([img_signal != 1])[0, ...]

# Erode the noise image
kernel = np.ones((7, 7), np.uint8)
img_noise = cv2.erode(img_noise, kernel, iterations=2)

# Dilate/erode the signal image (optional)
kernel = np.ones((3, 3), np.uint8)
img_signal = cv2.erode(img_signal, kernel, iterations=1)
img_signal = cv2.dilate(img_signal, kernel, iterations=1)

plt.figure()
plt.imshow(img_signal, cmap="gray")
plt.title('Signal mask \n using k-means')
plt.axis('off')

plt.figure()
plt.imshow(img_noise, cmap="gray")
plt.title('Noise mask \n using eroded k-means mask')
plt.axis('off')


# %% Visualize an example

# Choose image
img = img_cg_model[..., plot_slice]


def gradient_img(img):
    """gradient_img

    Args:
        img (Array): NxN image array
    """
    img = cv2.normalize(img, None, 0, 1.0, cv2.NORM_MINMAX, dtype=cv2.CV_32F)
    # set the kernel size, depending on whether we are using the Sobel
    # operator of the Scharr operator, then compute the gradients along
    # the x and y axis, respectively
    ksize = -1
    gX = cv2.Sobel(img, ddepth=cv2.CV_32F, dx=1, dy=0, ksize=ksize)
    gY = cv2.Sobel(img, ddepth=cv2.CV_32F, dx=0, dy=1, ksize=ksize)
    # the gradient magnitude images are now of the floating point data
    # type, so we need to take care to convert them back a to unsigned
    # 8-bit integer representation so other OpenCV functions can operate
    # on them and visualize them
    gX = cv2.convertScaleAbs(gX)
    gY = cv2.convertScaleAbs(gY)
    # combine the gradient representations into a single image
    img_grad = cv2.addWeighted(gX, 0.5, gY, 0.5, 0)

    print("The L2 norm of the Sobel gradient image is: " +
          str(np.linalg.norm(img_grad)))

    return img_grad


# Generate images
img_slice = cv2.normalize(img, None, 0, 1.0, cv2.NORM_MINMAX, dtype=cv2.CV_32F)
img_snr_slice = util.calculate_snr(img, img[img_noise == 1])

# Use total variation operator
G = sp.linop.FiniteDifference(img.shape)
img_grad_slice = G.N * img
img_grad_slice = cv2.normalize(
    img_grad_slice, None, 0, 1.0, cv2.NORM_MINMAX, dtype=cv2.CV_32F)
print("The L2 norm of G.N * img is: " + str(np.linalg.norm(img_grad_slice)))

# Use homemade function for visualization
img_grad_slice = gradient_img(img)

# Plot
plt.figure(figsize=(9, 3), dpi=100)
fig, ((ax1, ax2, ax3)) = plt.subplots(1, 3, sharey=True, gridspec_kw={'wspace': 0, 'hspace': 0},
                                      squeeze=True)
img_plot1 = ax1.imshow(np.rot90(img_slice, 3), cmap="gray", vmin=0, vmax=1)
img_plot2 = ax2.imshow(np.rot90(img_snr_slice, 3), cmap="jet", vmin=0, vmax=30)
img_plot3 = ax3.imshow(np.rot90(img_grad_slice, 3),
                       cmap="gnuplot", vmin=0, vmax=1.5)
ax1.axis('off')
ax2.axis('off')
ax3.axis('off')
plt.show()


# %% Generate a summary figure

# SNR limits
vmax_snr = np.percentile([util.calculate_snr(img_nufft, img_nufft[img_noise[..., plot_slice] == 1]),
                          util.calculate_snr(
                              img_cg, img_cg[img_noise[..., plot_slice] == 1]),
                          util.calculate_snr(img_cg_model, img_cg_model[img_noise[..., plot_slice] == 1])], 99.9)

# Set up subplot
plt.figure()
fig, axes = plt.subplots(nrows=6, ncols=3, figsize=(9, 18), dpi=100)

img = img_nufft[..., plot_slice]
ax = axes[0, 0]
im = ax.imshow(np.rot90(img, 3), cmap="gray")
ax.set_title(r'$(FD)^Hy$')
fig.colorbar(im, ax=ax)
ax.axis('off')
ax = axes[0, 1]
snr = ax.imshow(np.rot90(util.calculate_snr(
    img, img[img_noise == 1]), 3), cmap="jet", vmax=vmax_snr)
ax.set_title('snr, mean = %0.2f' % np.mean(
    util.calculate_snr(img, img[img_noise == 1])[img_noise != 1]))
fig.colorbar(snr, ax=ax)
ax.axis('off')
ax = axes[0, 2]
img_grad_slice = G.N * img
img_grad_slice = cv2.normalize(
    img_grad_slice, None, 0, 1.0, cv2.NORM_MINMAX, dtype=cv2.CV_32F)
grad = ax.imshow(np.rot90(gradient_img(img), 3),
                 cmap="gnuplot", vmin=0, vmax=1.5)
ax.set_title('grad, l2-norm = %0.2f' % np.linalg.norm(img_grad_slice))
fig.colorbar(grad, ax=ax)
ax.axis('off')

img = img_nufft_norm[..., plot_slice]
ax = axes[1, 0]
im = ax.imshow(np.rot90(img, 3), cmap="gray")
ax.set_title(r'$(FD)^Hy_n$')
fig.colorbar(im, ax=ax)
ax.axis('off')
ax = axes[1, 1]
snr = ax.imshow(np.rot90(util.calculate_snr(
    img, img[img_noise == 1]), 3), cmap="jet", vmax=vmax_snr)
ax.set_title('snr, mean = %0.2f' % np.mean(
    util.calculate_snr(img, img[img_noise == 1])[img_noise != 1]))
fig.colorbar(snr, ax=ax)
ax.axis('off')
ax = axes[1, 2]
img_grad_slice = G.N * img
img_grad_slice = cv2.normalize(
    img_grad_slice, None, 0, 1.0, cv2.NORM_MINMAX, dtype=cv2.CV_32F)
grad = ax.imshow(np.rot90(gradient_img(img), 3),
                 cmap="gnuplot", vmin=0, vmax=1.5)
ax.set_title('grad, l2-norm = %0.2f' % np.linalg.norm(img_grad_slice))
fig.colorbar(grad, ax=ax)
ax.axis('off')

img = img_cg[..., plot_slice]
ax = axes[2, 0]
im = ax.imshow(np.rot90(img, 3), cmap="gray")
ax.set_title(r'$||Fx-y||^2_2$')
fig.colorbar(im, ax=ax)
ax.axis('off')
ax = axes[2, 1]
snr = ax.imshow(np.rot90(util.calculate_snr(
    img, img[img_noise == 1]), 3), cmap="jet", vmax=vmax_snr)
ax.set_title('snr, mean = %0.2f' % np.mean(
    util.calculate_snr(img, img[img_noise == 1])[img_noise != 1]))
fig.colorbar(snr, ax=ax)
ax.axis('off')
ax = axes[2, 2]
img_grad_slice = G.N * img
img_grad_slice = cv2.normalize(
    img_grad_slice, None, 0, 1.0, cv2.NORM_MINMAX, dtype=cv2.CV_32F)
grad = ax.imshow(np.rot90(gradient_img(img), 3),
                 cmap="gnuplot", vmin=0, vmax=1.5)
ax.set_title('grad, l2-norm = %0.2f' % np.linalg.norm(img_grad_slice))
fig.colorbar(grad, ax=ax)
ax.axis('off')

img = img_cg_normalized[..., plot_slice]
ax = axes[3, 0]
im = ax.imshow(np.rot90(img, 3), cmap="gray")
ax.set_title(r'$||Fx-y_n||^2_2$')
fig.colorbar(im, ax=ax)
ax.axis('off')
ax = axes[3, 1]
snr = ax.imshow(np.rot90(util.calculate_snr(
    img, img[img_noise == 1]), 3), cmap="jet", vmax=vmax_snr)
ax.set_title('snr, mean = %0.2f' % np.mean(
    util.calculate_snr(img, img[img_noise == 1])[img_noise != 1]))
fig.colorbar(snr, ax=ax)
ax.axis('off')
ax = axes[3, 2]
img_grad_slice = G.N * img
img_grad_slice = cv2.normalize(
    img_grad_slice, None, 0, 1.0, cv2.NORM_MINMAX, dtype=cv2.CV_32F)
grad = ax.imshow(np.rot90(gradient_img(img), 3),
                 cmap="gnuplot", vmin=0, vmax=1.5)
ax.set_title('grad, l2-norm = %0.2f' % np.linalg.norm(img_grad_slice))
fig.colorbar(grad, ax=ax)
ax.axis('off')

img = img_cg_model[..., plot_slice]
ax = axes[4, 0]
im = ax.imshow(np.rot90(img, 3), cmap="gray")
ax.set_title(r'$||RFx-y||^2_2$')
fig.colorbar(im, ax=ax)
ax.axis('off')
ax = axes[4, 1]
snr = ax.imshow(np.rot90(util.calculate_snr(
    img, img[img_noise == 1]), 3), cmap="jet", vmax=vmax_snr)
ax.set_title('snr, mean = %0.2f' % np.mean(
    util.calculate_snr(img, img[img_noise == 1])[img_noise != 1]))
fig.colorbar(snr, ax=ax)
ax.axis('off')
ax = axes[4, 2]
img_grad_slice = G.N * img
img_grad_slice = cv2.normalize(
    img_grad_slice, None, 0, 1.0, cv2.NORM_MINMAX, dtype=cv2.CV_32F)
grad = ax.imshow(np.rot90(gradient_img(img), 3),
                 cmap="gnuplot", vmin=0, vmax=1.5)
ax.set_title('grad, l2-norm = %0.2f' % np.linalg.norm(img_grad_slice))
fig.colorbar(grad, ax=ax)
ax.axis('off')


img = img_cs_cg[..., plot_slice]
ax = axes[5, 0]
im = ax.imshow(np.rot90(img, 3), cmap="gray")
ax.set_title(r'$||RFx-y||^2_2 + g(x)$')
fig.colorbar(im, ax=ax)
ax.axis('off')
ax = axes[5, 1]
snr = ax.imshow(np.rot90(util.calculate_snr(
    img, img[img_noise == 1]), 3), cmap="jet", vmax=vmax_snr)
ax.set_title('snr, mean = %0.2f' % np.mean(
    util.calculate_snr(img, img[img_noise == 1])[img_noise != 1]))
fig.colorbar(snr, ax=ax)
ax.axis('off')
ax = axes[5, 2]
img_grad_slice = G.N * img
img_grad_slice = cv2.normalize(
    img_grad_slice, None, 0, 1.0, cv2.NORM_MINMAX, dtype=cv2.CV_32F)
grad = ax.imshow(np.rot90(gradient_img(img), 3),
                 cmap="gnuplot", vmin=0, vmax=1.5)
ax.set_title('grad, l2-norm = %0.2f' % np.linalg.norm(img_grad_slice))
fig.colorbar(grad, ax=ax)
ax.axis('off')
