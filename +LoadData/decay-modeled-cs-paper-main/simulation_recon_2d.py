# -*- coding: utf-8 -*-
"""
Title: 
    simulation_recon_2d.py

Agenda: 
    A .py script that uses compressed sensing methods to reconstruct 
non-cartesian MRI data.

Author: 
    Joseph Plummer - joseph.plummer@cchmc.org 

Creation date: 
    2023-04-28
    
Modification date: 
    2024-
    
"""

# %% Import packages
import cv2
from sklearn.cluster import KMeans
from functions import tv, optpoly, convexalg, phantom, util
import sigpy.mri as mr
import sigpy.plot as pl
import sigpy as sp
import numpy as np
import matplotlib.pyplot as plt
import random
from skimage.metrics import structural_similarity as ssim
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

# Image size
image_size = 128

# Choose spiral or radial sampling
sampling_type = str('spiral')
# sampling_type = str('radial')

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
    lamda_1 = 9e-7
    lamda_2 = 5e-6
    lamda_l2 = 1e-5
    rho = 5e1
    ptol = 1e-2
    num_normal = 4

# %% Set up a simulated phantom

# Phantom
x = phantom.lungphant_2d(image_size)

# Arbitrarily define mask region
img_mask = np.zeros_like(x)
img_mask[x > np.percentile(x, 10)] = 1

# Simulate defects: randomly generate 10 circles with random radii and positions
np.random.seed(123)
for i in range(10):
    # Generate a random radius between 5 and 20 pixels
    r = np.random.randint(3, 10)

    # Generate a random position between (r, r) and (image_size-r, image_size-r)
    rx = np.random.randint(0, image_size//2-r)
    ry = np.random.randint(image_size//2, image_size-r)

    # Create a meshgrid of x and y values
    xx, yy = np.mgrid[:image_size, :image_size]

    # Create a boolean mask of pixels inside the circle
    circle = (xx - rx) ** 2 + (yy - ry) ** 2 <= r ** 2

    # Set the pixels inside the circle to a random value between 0.5 and 1.0
    x[circle] = np.random.uniform(0.0, 0.5)
    x *= img_mask

plt.imshow(x, cmap="gray", vmin=0, vmax=1)
plt.title('Simulated phantom')
plt.colorbar()
plt.axis('off')

# %% Set up coordinates

if sampling_type == 'spiral':
    # Set up spiral sampling
    N_proj = int(image_size//3)  # Number of projections
    traj = mr.spiral(fov=0.3,  # field of view in meters.
                     N=image_size,  # effective matrix shape.
                     # undersampling factor in freq encoding direction.
                     f_sampling=0.39,
                     R=0.33,  # undersampling factor.
                     ninterleaves=N_proj,  # number of spiral interleaves.
                     alpha=3,  # variable density factor.
                     gm=40E-3,  # maximum gradient apltitude (T/m).
                     sm=180,  # maximum slew rate (T/m/s).
                     gamma=73.997E6)  # gyromagnetic ratio in rad/T/s.
    traj = np.reshape(traj, (N_proj, int(len(traj)/N_proj), 2))
    subset = int(N_proj * subsample_pc / 100)
    random.seed(123)
    subset_range = random.sample(
        range(0, N_proj), subset)  # Randomize subsampling
    traj = traj[subset_range, ...]
    N_proj = np.shape(traj)[0]  # Update to adjust with subsampling

elif sampling_type == 'radial':
    N_proj = int(np.pi * image_size)  # Number of projections
    # Sample slightly more times to prevent ring at edge of FOV
    N_samp = int(0.75 * image_size)
    N_dimensions = 2
    traj = mr.radial(coord_shape=(N_proj, N_samp, N_dimensions),  # e.g. (512, 96, 2) = 512 proj, 96 samples, 2 dimensions
                     img_shape=np.shape(x),
                     golden=True)
    subset = int(N_proj * subsample_pc / 100)
    traj = traj[:subset, ...]
    N_proj = np.shape(traj)[0]  # Update to adjust with subsampling

# Plot trajectories:
plt.figure(figsize=(5, 5), dpi=100)
N_visual = np.shape(traj)[0]
color = iter(plt.cm.viridis(np.linspace(0, 1, N_visual)))
for i in np.linspace(0, N_visual-1, N_visual):
    i = int(i)
    c = next(color)
    plt.scatter(traj[i, :, 0],
                traj[i, :, 1], color=c, s=1, marker='.')
plt.ylabel('$k_y$', fontsize=18)
plt.xlabel('$k_x$', fontsize=18)
plt.xticks([])
plt.yticks([])
plt.axis('square')
plt.show()

# Verify inputs
N_dimensions = np.shape(traj)[2]  # Number of dimensions
N_samp = np.shape(traj)[1]  # Samples per projection
N_proj = np.shape(traj)[0]  # Number of projections

# %% Generate k-space data using NUFFT

# Make NUFFT linear operator
F = sp.linop.NUFFT(np.shape(x),
                   traj,
                   oversamp=1.25,
                   width=4,
                   toeplitz=True)

# Simulate T2* decay
t2_star_xe = 18  # ms
readout = 22  # ms
dwell_time = readout/N_samp
relaxation = np.zeros((N_samp,))
for i in range(N_samp):
    relaxation[i] = np.exp(-i*dwell_time/t2_star_xe)

# Simulate magnetization decay
flip_angle = np.arctan(np.sqrt(2/(N_proj)))*1.8
data = np.zeros((N_proj, N_samp), dtype=complex)
for i in range(N_proj):
    img_temp = x * np.cos(flip_angle)**i
    data_temp = F * img_temp
    data[i, :] = data_temp[i, :]

# Multiply by T2* relaxation vector
data *= relaxation

# Visualize
fig, (ax1, ax2, ax3, ax4) = plt.subplots(
    4, sharex=False, figsize=(6, 6), dpi=100)
ax1.plot(abs(data[:, 0]), color='r')
ax1.set_ylabel(r'$M_{z}$')
ax1.set_title("magnetization decay")
ax1.set_xlabel("projection number")
ax2.plot(abs(data[0, :]), color='b')
ax2.set_ylabel(r'$M_{xy}$')
ax2.set_xlabel('sample number')
ax2.set_title("fid")

# Manuscript
fig, (ax1) = plt.subplots(
    1, sharex=False, figsize=(6, 3), dpi=100)
ax1.plot(abs(util.generate_complex_normal(
    0, np.max(abs(data))*0.01, np.shape(data[0, :]))), color='y')
ax1.set_ylabel('intensity')
ax1.set_title("noise")
ax1.set_xlabel("sample number")

fig, (ax1) = plt.subplots(
    1, sharex=False, figsize=(6, 3), dpi=100)
ax1.plot(abs(data[:, 0]), color='r', linewidth=4)
ax1.set_ylabel('$k_0$ intensity')
ax1.set_title("magnetization")
ax1.set_xlabel("projection number")

# Normalize data to max k0
mean_k0 = np.mean(abs(data[:, 0]))

# Add noise to data
noise_mean, noise_sd = 0, np.max(abs(data))*0.001
noise = util.generate_complex_normal(noise_mean, noise_sd, np.shape(data))

data += noise
ax3.plot(abs(data[:, 0]), color='r')
ax3.set_ylabel('k0 intensity')
ax3.set_title("normalized magnetization decay with noise added")
ax3.set_xlabel("projection number")
ax4.plot(abs(data[0, :]), color='b')
ax4.set_ylabel('k-space intensity')
ax4.set_xlabel('sample number')
ax4.set_title("fid with noise added")
fig.tight_layout(pad=1.0, w_pad=0.5, h_pad=1.0)

# %% Calculate density compensation

# Pipe-Menon DCF settings
beta = 8  # Default = 8
width = 3  # Default = 4, choose ~3 for high SNR, choose 5~8 for sharpness

try:
    dcf = mr.pipe_menon_dcf(traj, img_shape=np.dot(2, np.shape(x)),
                            beta=beta, width=width, device=device)
    D = sp.linop.Multiply(F.oshape, dcf**0.5)
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
k = np.zeros((N_proj, N_samp))
for i in range(N_proj):
    k[i, :] = np.cos(flip_angle)**i * relaxation

# Visualize
fig, (ax1, ax2) = plt.subplots(2, sharex=False)
# First few points of kspace_data2 have zero signal due to LNA issue
ax1.plot(abs(k[:, 0]), color='g')
ax1.set_ylabel(r'$M_{z}$')
ax1.set_title("magnetization decay from excitation and $T_1$")
ax1.set_xlabel("projection number")
ax2.plot(abs(k[0, :]), color='y')
ax2.set_ylabel(r'$M_{xy}$')
ax2.set_xlabel('sample number')
ax2.set_title("magnetization decay from $T_2^*$")
fig.tight_layout(pad=1.0, w_pad=0.5, h_pad=1.0)

# %% Inverse NUFFT reconstruction

data = mvd(data)

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
for i in range(N_proj):
    data_normalized[i, :] = data[i, :] * mean_k0 / data[i, 0]
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
img_cg = abs(mvc(convexalg.cg(num_iters=100,
                              ptol=0,
                              A=A, b=b,
                              lamda=lamda_l2,  # l2-norm constraint for additional regularization
                              verbose=True,
                              draw_output=True,
                              save=None)))

plt.figure()
plt.imshow(img_cg, cmap="gray")
plt.colorbar()
plt.title('Reconstruction using iterative NUFFT')
plt.axis('off')

# Declare if density compensated
K = sp.linop.Multiply(F.oshape, k**0)
data_normalized = np.zeros_like(data_normalized)
for i in range(N_proj):
    data_normalized[i, :] = data[i, :] * mean_k0 / data[i, 0]
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
img_cg_normalized = abs(mvc(convexalg.cg(num_iters=100,
                                         ptol=0,
                                         A=A, b=b,
                                         lamda=lamda_l2,  # l2-norm constraint for additional regularization
                                         verbose=True,
                                         draw_output=True,
                                         save=None)))

plt.figure()
plt.imshow(img_cg_normalized, cmap="gray")
plt.colorbar()
plt.title('Reconstruction using iterative NUFFT')
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
img_cg_model = abs(mvc(convexalg.cg(num_iters=100,
                                    ptol=0,
                                    A=A, b=b,
                                    lamda=lamda_l2,  # l2-norm constraint for additional regularization
                                    verbose=True,
                                    draw_output=True,
                                    save=None)))

plt.figure()
plt.imshow(img_cg_model, cmap="gray")
plt.colorbar()
plt.title('Reconstruction using iterative NUFFT')
plt.axis('off')

# Define regularizing linear operators and their proximal operators
W = sp.linop.Wavelet(x.shape, wave_name="db4")
def g1(x): return lamda_1 * xp.linalg.norm(W(x).ravel(), ord=1)


prox_g1 = sp.prox.UnitaryTransform(
    sp.prox.L1Reg(W.oshape, lamda_1), W)

prox_g2 = tv.ProxTV(A.ishape, lamda_2)
def g2(x): return lamda_2 * xp.linalg.norm(prox_g2.G(x))


# Make list of objectives and proximal operators
lst_g = [g1, g2]
lst_proxg = [prox_g1, prox_g2]
# lst_g = [g2]
# lst_proxg = [prox_g2]

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
# TODO: the above method uses some advanced preconditioning to solve ||Ax-b||^2 faster. 
# It makes negligible difference to this script, but its worth borrowing the code for high res/high channel-number 3D applications. 
# If interested in optimizing, you need to adjust rho, num_iters, num_normal.
# Read more about this here: https://doi.org/10.1137/22M1530355

# %% Normalize image


def norm_img(image):
    """
    Normalize an image between 0 and 1.

    Parameters:
        image (numpy.ndarray): The input image as a NumPy array.

    Returns:
        numpy.ndarray: The normalized image.
    """
    # Find the minimum and maximum pixel values in the image
    min_value = np.min(image)
    max_value = np.max(image)

    # Normalize the image to the range [0, 1]
    normalized_image = (image - min_value) / (max_value - min_value)

    return normalized_image

# %% Summary plots for as shown in paper


plt.style.use('default')

# Choose a line profile index
line_idx = 40
def lineprof(x): return x[line_idx, :]


img_tmp = img_cs_cg
# img_tmp = img_cs_pi
# img_tmp = img_cg
# img_tmp = img_cg_normalized
# img_tmp = img_cg_model
# img_tmp = img_nufft
# img_tmp = img_nufft_norm
img_tmp = cv2.normalize(img_tmp, None, 0, 1.0,
                        cv2.NORM_MINMAX, dtype=cv2.CV_32F)
x = cv2.normalize(x, None, 0, 1.0, cv2.NORM_MINMAX, dtype=cv2.CV_32F)
color = "green"

# Mean SNR
snr_mean = np.mean(util.calculate_snr(
    img_tmp, img_tmp[img_mask == 0])[img_mask != 0])
# text_label = str(r'$|\Delta|=$' + str(np.round(np.linalg.norm((x-img_tmp), 1), 1)) +
#                  r', $\overline{SNR}=$ ' + str(np.round(snr_mean, 1)))
text_label = str(
    r'$RMSE=$' + str(np.round(util.calculate_rmse(x, img_tmp), 3)))


# Set up subplot
plt.figure(dpi=300)
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 6),
                         gridspec_kw={'wspace': 0.0, 'hspace': 0},
                         squeeze=True)

ax = axes[0]
im = ax.imshow(img_tmp, cmap="gray")
ax.plot((0, image_size), (line_idx, line_idx),
        color, linestyle="solid", linewidth=3)
# ax.set_title('object')
# fig.colorbar(im, ax=ax)
ax.axis('off')
ax = axes[1]
ax.plot(lineprof(norm_img(x)), 'k', linestyle="dashed",
        linewidth=4, label=r'$\hat{x}$')
ax.plot(lineprof(norm_img(img_tmp)), color, linewidth=2)
# ax.legend(loc='upper right', fontsize=10)
# ax.set_title('line profile')
ax.set_xlim([0, image_size])
ax.set_ylim([0, 0.75])
ax.tick_params(axis='both', which='major', labelsize=16)
ax.set_xticks([])
ax.set_yticks([])
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.axis('off')
ax.set_xlabel(
    r'$|diff|=$'+str(np.round(np.linalg.norm(x-img_tmp, 1), 1)), fontsize=16)
ax.text(0.5, 0.9, text_label,
        transform=ax.transAxes,
        horizontalalignment='center',
        verticalalignment='center',
        size=24)
