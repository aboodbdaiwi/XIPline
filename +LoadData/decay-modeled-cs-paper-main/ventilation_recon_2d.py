# -*- coding: utf-8 -*-
"""
Title: 
    ventilation_recon_2d.py

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
import cv2
from sklearn.cluster import KMeans
from functions import tv, optpoly
from functions import convexalg, util
import sigpy.mri as mr
import sigpy.plot as pl
import sigpy as sp
import numpy as np
import matplotlib.pyplot as plt
from skimage.transform import resize
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

# Proximal gradient descent settings
if use_dcf:
    num_iters = 50
    lamda_1 = 8e-5
    lamda_2 = 1e-2
    lamda_l2 = 1e-5
    rho = 1
    ptol = 1e-1
    num_normal = 4
else:
    # More iterations required without left-preconditioning (i.e. DCF)
    num_iters = 200
    lamda_1 = 8e-7  # Default 5e-7
    lamda_2 = 1e-5  # Default 1e-7
    lamda_l2 = 1e-5
    rho = 5e1
    ptol = 1e-2
    num_normal = 4

# %% Import MRI k-space data/coordinates

# Filenames
data_file = "data/ventilation-2d/ksp.npy"
coords_file = "data/ventilation-2d/coord.npy"

# Load data (channels, projections, samples)
data = np.array(np.load(data_file))

# Select a slice
imslice = 8
data = data[imslice, :, :]

# Add noise to data
# noise_mean, noise_sd = 0, np.max(abs(data))*0.005
# noise = util.generate_complex_normal(noise_mean, noise_sd, np.shape(data))
# data += noise

# Visualize
fig, (ax1, ax2) = plt.subplots(2, sharex=False)
# First few points of kspace_data2 have zero signal due to LNA issue
ax1.plot(abs(data[:, 0]), color='r')
ax1.set_ylabel('$k_0$ intensity')
ax1.set_title("magnetization decay")
ax1.set_xlabel("projection number")
ax2.plot(abs(data[0, :]), color='b')
ax2.set_ylabel('k-space intensity')
ax2.set_xlabel('sample number')
ax2.set_title("fid")
fig.tight_layout(pad=1.0, w_pad=0.5, h_pad=1.0)

# Load coords (projections, samples, dimensions)
coords = np.array(np.load(coords_file))

# Verify inputs
N_dimensions = np.shape(coords)[2]  # Number of dimensions
N_samp = np.shape(data)[1]  # Samples per projection
N_proj = np.shape(data)[0]  # Number of projections
N_channels = 1  # Number of coil channels

# Plot trajectories:
plt.figure(figsize=(5, 5), dpi=100)
N_visual = N_proj
color = iter(plt.cm.viridis(np.linspace(0, 1, N_visual)))
for i in np.linspace(0, N_visual-1, N_visual):
    i = int(i)
    c = next(color)
    plt.scatter(coords[i, :, 0],
                coords[i, :, 1], color=c, s=1, marker='.')
plt.ylabel('$k_y$', fontsize=18)
plt.xlabel('$k_x$', fontsize=18)
plt.xticks([])
plt.yticks([])
plt.axis('square')
plt.show()

# %% Reconstruction settings

# Image reconstruction sizes
scan_resolution = 120
recon_resolution = int(scan_resolution * 3)

# Make manipulations to reconstructed volume
alpha = 1  # Overgridding factor (alpha = 0.5 means reconstructing on 2x grid)
resize_factor = recon_resolution / scan_resolution
image_size = int(resize_factor * scan_resolution / alpha)
image_shape = (int(image_size), int(image_size))
coords_resize = image_size * coords / resize_factor

# %% Subset the data to experiment with undersampling

# Select only a subset of trajectories and data
N_proj_subset = int(N_proj * subsample_pc / 100)
subset = range(0, N_proj_subset)
data_subset = data[subset, :]
coords_subset = coords_resize[subset, ...]

# %% Calculate density compensation

# Pipe-Menon DCF settings
beta = 8  # Default = 8
width = 4  # Default = 4, choose ~3 for high SNR, choose 5~8 for sharpness

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


# %% Generate a hyperpolarized weighting for the inverse NUFFT

# Simulate T2* decay
t2_star_xe = 18  # ms
readout = 10  # ms
dwell_time = readout/N_samp
relaxation = np.zeros((N_samp,))
for i in range(N_samp):
    relaxation[i] = np.exp(-i*dwell_time/t2_star_xe)

# Estimate magnetization decay due to RF excitations
k0 = abs(data_subset[:, 0])
k = np.zeros((N_proj_subset, N_samp), dtype=complex)
for i in range(N_proj_subset):
    k[i, :] = k0[i] * relaxation

# Apply an optional frequency offset
phase_k0 = np.angle(np.mean(data[:, 0]))
phase_k1 = np.angle(np.mean(data[:, 1]))
phase_acc = phase_k1 - phase_k0
print("The phase accumulation that built up between the first two points: " + str(phase_acc))
freq_offset = phase_acc/(2*np.pi*dwell_time*1e-3)
print("The estimated frequency offset from the phase accumulation is: " + str(freq_offset))

freq_offset = -15  # -25 for ksp12 # -10 for ksp1 # -40 for ksp5 # +30 # for ksp9 # Hz
for i in range(N_samp):
    scale_factor = np.exp(-2*np.pi*1j*i*freq_offset*dwell_time*1e-3)
    k[:, i] *= scale_factor

# Visualize
fig, (ax1, ax2) = plt.subplots(2, sharex=False)
# First few points of kspace_data2 have zero signal due to LNA issue
ax1.plot(abs(k[:, 0]), color='g')
ax1.set_ylabel('signal intensity')
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
        sens_map = mr.app.JsenseRecon(data_subset, coord=coords_resize, weights=dcf,
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
    sens_map = np.ones((image_size, image_size), dtype=int)
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
A_dcf = K * D * F
b_dcf = data * dcf**0.5 / xp.linalg.norm(data * dcf**0.5)
img_nufft = mvc(abs(A_dcf.H * b_dcf))

# Visualize reconstructed phantom:
plt.figure()
plt.imshow(img_nufft, cmap="gray")
plt.colorbar()
plt.title('Reconstruction using inverse NUFFT \n (no norm)')
plt.axis('off')

# Compute a single x = A.H b operation (i.e. inverse NUFFT)
data_normalized = np.zeros_like(data)
K = sp.linop.Multiply(F.oshape, k**0)
for i in range(N_proj_subset):
    data_normalized[i, :] = data[i, :] / k0[i]
A_dcf = K * D * F
b_dcf = data_normalized * dcf**0.5 / xp.linalg.norm(data_normalized * dcf**0.5)
img_nufft_norm = mvc(abs(A_dcf.H * b_dcf))

# Visualize reconstructed phantom:
plt.figure()
plt.imshow(img_nufft_norm, cmap="gray")
plt.colorbar()
plt.title('Reconstruction using inverse NUFFT \n (norm data)')
plt.axis('off')

# Compute a single x = A.H b operation (i.e. inverse NUFFT)
K = sp.linop.Multiply(F.oshape, k**1)
A_dcf = K * D * F
b_dcf = data * dcf**0.5 / xp.linalg.norm(data * dcf**0.5)
img_nufft_model = mvc(abs(A_dcf.H * b_dcf))

# Visualize reconstructed phantom:
plt.figure()
plt.imshow(img_nufft_model, cmap="gray")
plt.colorbar()
plt.title('Reconstruction using inverse NUFFT \n (model decay)')
plt.axis('off')


# %% Iterative NUFFT reconstruction

data = mvd(data)

# Declare if density compensated
K = sp.linop.Multiply(F.oshape, k**0)
if use_dcf:
    A = K * D * F
    data_norm = xp.linalg.norm(data * dcf**0.5)
    b = data * dcf**0.5 / data_norm
else:
    A = K * F
    data_norm = xp.linalg.norm(data)
    b = data / data_norm

# Calculate maximum eigenvalue
LL = sp.app.MaxEig(A.N, dtype=xp.complex64,
                   device=device).run() * 1.01

# Divide by svd eigenvalue to ensure maximum eigenvalue falls between [0,1]
A *= np.sqrt(1/LL)

# Solve ||Ax - b||^2 + lamda||x||^2 using conjugate gradient descent
img_cg = abs(mvc(convexalg.cg(num_iters=20,
                              ptol=0,
                              A=A, b=b,
                              lamda=lamda_l2,  # l2-norm constraint for additional regularization
                              verbose=True,
                              draw_output=True,
                              save=None)))

plt.figure()
plt.imshow(img_cg, cmap="gray")
plt.colorbar()
plt.title('Reconstruction using iterative NUFFT \n (no norm)')
plt.axis('off')

# Declare if density compensated
K = sp.linop.Multiply(F.oshape, k**0)
data_normalized = np.zeros_like(data)
for i in range(N_proj_subset):
    data_normalized[i, :] = data[i, :] / k0[i]

if use_dcf:
    A = K * D * F
    data_norm = xp.linalg.norm(data_normalized * dcf**0.5)
    b = data_normalized * dcf**0.5 / data_norm
else:
    A = K * F
    data_norm = xp.linalg.norm(data_normalized)
    b = data_normalized / data_norm

# Calculate maximum eigenvalue
LL = sp.app.MaxEig(A.N, dtype=xp.complex64,
                   device=device).run() * 1.01

# Divide by svd eigenvalue to ensure maximum eigenvalue falls between [0,1]
A *= np.sqrt(1/LL)

# Solve ||Ax - b||^2 + lamda||x||^2 using conjugate gradient descent
img_cg_normalized = abs(mvc(convexalg.cg(num_iters=20,
                                         ptol=0,
                                         A=A, b=b,
                                         lamda=lamda_l2,  # l2-norm constraint for additional regularization
                                         verbose=True,
                                         draw_output=True,
                                         save=None)))

plt.figure()
plt.imshow(img_cg_normalized, cmap="gray")
plt.colorbar()
plt.title('Reconstruction using iterative NUFFT \n (norm data)')
plt.axis('off')

# Declare if density compensated
K = sp.linop.Multiply(F.oshape, k**1)
if use_dcf:
    A = K * D * F
    data_norm = xp.linalg.norm(data * dcf**0.5)
    b = data * dcf**0.5 / data_norm
else:
    A = K * F
    data_norm = xp.linalg.norm(data)
    b = data / data_norm

# Calculate maximum eigenvalue
LL = sp.app.MaxEig(A.N, dtype=xp.complex64,
                   device=device).run() * 1.01

# Divide by svd eigenvalue to ensure maximum eigenvalue falls between [0,1]
A *= np.sqrt(1/LL)

# Solve ||Ax - b||^2 + lamda||x||^2 using conjugate gradient descent
img_cg_model = abs(mvc(convexalg.cg(num_iters=20,
                                    ptol=0,
                                    A=A, b=b,
                                    lamda=lamda_l2,  # l2-norm constraint for additional regularization
                                    verbose=True,
                                    draw_output=True,
                                    save=None)))

plt.figure()
plt.imshow(img_cg_model, cmap="gray")
plt.colorbar()
plt.title('Reconstruction using iterative NUFFT \n (model norm)')
plt.axis('off')


# Define regularizing linear operators and their proximal operators
W = sp.linop.Wavelet(image_shape, wave_name="db4")
def g1(x): return lamda_1 * xp.linalg.norm(W(x).ravel(), ord=1)


prox_g1 = sp.prox.UnitaryTransform(
    sp.prox.L1Reg(W.oshape, lamda_1), W)

prox_g2 = tv.ProxTV(A.ishape, lamda_2)
def g2(x): return lamda_2 * xp.linalg.norm(prox_g2.G(x))


# Make list of objectives and proximal operators
lst_g = [g1, g2]
lst_proxg = [prox_g1, prox_g2]

# Solve ||Ax - b||^2 + lamda_1 |Wx| + lamda_2 |Gx| using ADMM and CG
img_cs_cg = mvc(abs(convexalg.admm(num_iters=num_iters,
                                   ptol=ptol,
                                   num_normal=num_normal,
                                   A=A, b=b,
                                   lst_proxg=lst_proxg,
                                   rho=rho,
                                   lst_g=lst_g,
                                   method="cg",
                                   verbose=True,
                                   draw_output=True)))

plt.figure()
plt.imshow(img_cs_cg, cmap="gray")
plt.colorbar()
plt.title('Reconstruction using CS: ADMM w/ CG')
plt.axis('off')

# Solve ||Ax - b||^2 + lamda_1 |Wx| + lamda_2 |Gx| using ADMM and PI
img_cs_pi = mvc(abs(convexalg.admm(num_iters=num_iters,
                                   ptol=ptol,
                                   num_normal=num_normal,
                                   A=A, b=b,
                                   lst_proxg=lst_proxg,
                                   rho=rho,
                                   lst_g=lst_g,
                                   method="pi",
                                   verbose=True,
                                   draw_output=True)))
plt.figure()
plt.imshow(img_cs_pi, cmap="gray")
plt.colorbar()
plt.title('Reconstruction using CS: ADMM w/ PI')
plt.axis('off')


# %% Try a signal mask using K means

# Choose an image
img = cv2.normalize(img_cg_model, None, 0, 1.0,
                    cv2.NORM_MINMAX, dtype=cv2.CV_32F)

# Perform clusting
N_clusters = 3
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
img_noise = np.uint8([img_signal != 1])[0, ...]
img_signal = np.uint8([img_signal])[0, ...]  # This line rounds to 0s and 1s

# Erode the noise image
kernel = np.ones((7, 7), np.uint8)
img_noise = cv2.erode(img_noise, kernel, iterations=5)

# Dilate the signal image (optional)
kernel = np.ones((3, 3), np.uint8)
img_signal = cv2.dilate(img_signal, kernel, iterations=1)

plt.figure()
plt.imshow(img_signal, cmap="gray")
plt.title('Signal mask \n using dilated K-means mask')
plt.axis('off')

plt.figure()
plt.imshow(img_noise, cmap="gray")
plt.title('Noise mask \n using eroded K-means mask')
plt.axis('off')

# %% Visualize an example


def gradient_img(img):
    """gradient_img

    Args:
        img (Array): NxN image array
    """
    img = cv2.normalize(img, None, 0, 3.0, cv2.NORM_MINMAX, dtype=cv2.CV_32F)
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
          str(np.linalg.norm(img_grad[img_signal == 1])))
    print("The mean of the Sobel gradient image is: " +
          str(np.mean(img_grad[img_signal == 1])))

    return img_grad


# Choose image
img = img_cg_model

# Generate images
img_slice = cv2.normalize(img, None, 0, 1.0, cv2.NORM_MINMAX, dtype=cv2.CV_32F)
img_snr_slice = util.calculate_snr(img, np.ravel(img[img_noise == 1]))

# Use total variation operator
G = sp.linop.FiniteDifference(S.oshape)
img_grad_slice = G.N * img
img_grad_slice = cv2.normalize(
    img_grad_slice, None, 0, 1.0, cv2.NORM_MINMAX, dtype=cv2.CV_32F)
print("The L2 norm of G.N * img is: " +
      str(np.linalg.norm(img_grad_slice[img_signal == 1])))

# Use homemade function for visualization
img_grad_slice = gradient_img(img)

# Flip image for correct dimension
img_slice = np.flipud(img_slice)
img_snr_slice = np.flipud(img_snr_slice)
img_grad_slice = np.flipud(img_grad_slice)

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

# Plot images individually
plt.figure(figsize=(2, 2), dpi=300)
plt.imshow(np.rot90(img_slice, 3), cmap="gray", vmin=0, vmax=1)
plt.axis('off')
plt.show()

x_min = int(image_size*0.2)
x_max = int(image_size*0.5)
y_min = int(image_size*0.15)
y_max = int(image_size*0.45)
plt.figure(figsize=(2, 2), dpi=300)
plt.imshow(np.rot90(img_grad_slice, 3)[x_min:x_max, y_min:y_max],
           cmap="gray", vmin=0, vmax=10)
# The vmin/vmax limits depend on the limits of the cv.normalize function inside gradient_img()
plt.axis('off')
plt.show()

plt.figure(figsize=(2, 2), dpi=300)
plt.imshow(np.rot90(img_slice, 3)[x_min:x_max, y_min:y_max],
           cmap="gray", vmin=0, vmax=1)
# The vmin/vmax limits depend on the limits of the cv.normalize function inside gradient_img()
plt.axis('off')
plt.show()

plt.figure(figsize=(2, 2), dpi=300)
plt.imshow(np.rot90(img_snr_slice, 3), cmap="jet", vmin=0, vmax=47)
plt.axis('off')
plt.show()
plt.close()

plt.figure(figsize=(2, 2), dpi=300)
im = plt.imshow(np.rot90(np.log(img_slice), 3), cmap="jet", vmin=-4, vmax=0)
plt.axis('off')
# plt.show()
# plt.close()

cbar = fig.colorbar(im, location='bottom')

cbar.outline.set_edgecolor('black')
cbar.set_ticks([-4, -2, 0])
cbar.set_ticklabels(["-4.0", "-2.0", "0.0"])
cbar.ax.tick_params(labelsize=16)
# plt.figure(figsize=(10, 10), dpi=300)
# im = plt.imshow(img_weights_slice,
#                 cmap="gnuplot", vmin=0, vmax=1)
# plt.axis('off')


# %% Generate a summary figure

# SNR limits
vmax_snr = np.percentile([util.calculate_snr(img_nufft, img_nufft[img_noise == 1]),
                          util.calculate_snr(img_cg, img_cg[img_noise == 1]),
                          util.calculate_snr(img_cg_model, img_cg_model[img_noise == 1])], 99.9)

# Set up subplot
plt.figure()
fig, axes = plt.subplots(nrows=6, ncols=3, figsize=(9, 18), dpi=100)

img = img_nufft
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

img = img_nufft_norm
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

img = img_cg
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

img = img_cg_normalized
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

img = img_cg_model
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


img = img_cs_cg
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

# %%
