# -*- coding: utf-8 -*-
"""
Title: 
    gasex_recon_3d.py

Agenda: 
    A .py script that uses compressed sensing methods to reconstruct 
non-cartesian MRI data.

Author: 
    Joseph Plummer - joseph.plummer@cchmc.org 

Creation date: 
    2023-07-24
    
Modification date: 
    --
    
"""

# %% Import packages
from scipy.signal import decimate
import nibabel as nib
import os
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


# Folder
folder = "data/gasex-3d/"

# Density compensation
use_dcf = False
analytical = False  # A useful preconditioner for 3D radials at full resolution

# Normalize data
normalize_data = False

# Multiply F by signal weighting
gamma = 1  # 1 = full HP weighting on F, 0 = no weighting on F

# Subsample
subsample_pc_gas = 100  # percentage
subsample_pc_diss = 100

# Select type of under-sampling
subset_type = "first n"  # Useful for only using first high signal projections
subset_type = "every other n"  # Useful for preserving spread of magnetization decay

# DCF settings
beta = 8  # Default = 8
width = 4  # Default = 4, choose ~3 for high SNR, choose 5~8 for sharpness
crd_osamp = 4  # Rec: 1-4, set to 4 for low res recon
img_osamp = 1  # Rec: 1.25, Make higher than crd_osamp for sharpness/lower for signal

# Proximal gradient descent settings
if use_dcf:
    if analytical:
        num_iters = 5
        lamda_1 = 1e-6
        lamda_2 = 1e-6
        lamda_l2 = 1e-3
        lamda_3 = 1e-4
        zeta = 1
        lamda_4 = 5e-2
        eta = 1e-5
        rho = 1
        ptol = 1e-3
        num_normal = 5
    else:
        num_iters = 10
        lamda_1 = 5e-4
        lamda_2 = 1e-2
        lamda_l2 = 1e-2
        lamda_3 = 1e-4
        zeta = 1
        lamda_4 = 5e-2
        eta = 1e-5
        rho = 1
        ptol = 1e-3
        num_normal = 5
else:
    # More iterations required without left-preconditioning (i.e. DCF)
    # 64 matrix
    num_iters = 100  # Set to very high if ill-conditioned (i.e. under-sampled 3D radial)
    lamda_1 = 5e-7
    lamda_2 = 3e-5
    lamda_l2 = 1e-5
    lamda_3 = 5e-5  # 5e-5 for the radial CPIR scans
    zeta = 1
    lamda_4 = 1e-4  # TODO: Find optimal
    eta = 5e-3  # 1e-2 default
    rho = 5e2  # 5e2 for CG # TODO: Experiment with PI recon vs CG
    ptol = 1e-5
    num_normal = 20  # Set to larger value if ill-conditioned (i.e. under-sampled 3D radial)


# %% Import MRI k-space data/coordinates

data_file = "ksp_g.npy"
coords_file = "coord_g.npy"

# Load data (channels, projections, samples)
data = np.array(np.load(folder + data_file))

# Visualize
fig, (ax1, ax2) = plt.subplots(2, sharex=False)
ax1.plot(abs(data[0, :, 0]), color='r')
ax1.set_ylabel('$k_0$ intensity')
ax1.set_title("magnetization decay")
ax1.set_xlabel("projection number")
ax2.plot(abs(data[0, 0, :]), color='b')
ax2.plot(abs(data[0, 50, :]), color='b')
ax2.plot(abs(data[0, 100, :]), color='b')
ax2.set_ylabel('k-space intensity')
ax2.set_xlabel('sample number')
ax2.set_title("fid (three arbitrary projections)")
fig.tight_layout(pad=1.0, w_pad=0.5, h_pad=1.0)

# Normalize data to max k0
k0 = abs(data[0, :, 0])
data_normalized = np.zeros_like(data)
for i in range(np.shape(data)[1]):
    data_normalized[0, i, :] = data[0, i, :] / k0[i]

if normalize_data:
    data = data_normalized

# Load coords (projections, samples, dimensions)
coords = np.array(np.load(folder + coords_file))
print(np.max(coords))
print(abs(data[-1, 0, -1]))

# Flip x-y-z dimensions if necessary (use for gas recon)
coords[:, :, [2, 1]] = coords[:, :, [1, 2]]

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

# %% Downsample k-space data if oversampled along readout (Philips R59 does this)


def moving_average(data, window_size=4):
    return np.convolve(data, np.ones(window_size)/window_size, mode='same')


if N_samp == 232:

    for ii in range(N_proj):
        data[0, ii, :] = moving_average(
            data[0, ii, :], window_size=int(232/58))
        coords[ii, :, 0] = moving_average(
            coords[ii, :, 0], window_size=int(232/58))
        coords[ii, :, 1] = moving_average(
            coords[ii, :, 1], window_size=int(232/58))
        coords[ii, :, 2] = moving_average(
            coords[ii, :, 2], window_size=int(232/58))

    # Downsampling
    data = decimate(data, q=int(232/58))
    N_samp = data.shape[-1]
    coordsx = np.zeros((N_proj, N_samp))
    coordsy = np.zeros((N_proj, N_samp))
    coordsz = np.zeros((N_proj, N_samp))
    for ii in range(N_proj):
        coordsx[ii, :] = decimate(coords[ii, :, 0], q=int(232/58))
        coordsy[ii, :] = decimate(coords[ii, :, 1], q=int(232/58))
        coordsz[ii, :] = decimate(coords[ii, :, 2], q=int(232/58))

    coords = np.zeros((N_proj, N_samp, 3))
    coords[..., 0] = coordsx
    coords[..., 1] = coordsy
    coords[..., 2] = coordsz


# %% Reconstruction settings

# Image reconstruction sizes
pad_width = 3
scan_resolution = 64 + 2*pad_width 
recon_resolution = int(scan_resolution * 1)

# Make manipulations to reconstructed volume
# Overgridding factor (alpha = 0.5 means reconstructing on 2x grid)
alpha = 1
resize_factor = recon_resolution / scan_resolution
image_size = int(resize_factor * scan_resolution / alpha)
image_shape = (int(image_size), int(image_size), int(image_size))
coords_resize = image_size * coords / resize_factor

# %% Subset the data to experiment with undersampling

# Steady state
steady_state = 0

# Select only a subset of trajectories and data
if subset_type == "first n":
    subset = int(N_proj * subsample_pc_gas / 100)
    data_subset = data[:, steady_state:subset, :]
    coords_subset = coords_resize[steady_state:subset, ...]
else:
    subset = int(100 / subsample_pc_gas)
    data_subset = data[:, steady_state::subset, :]
    coords_subset = coords_resize[steady_state::subset, ...]

# Update number of projections
N_proj_subset = data_subset.shape[1]

# Update k0
k0_subset = abs(data_subset[0, :, 0])

# %% Deal with NaN values
data_subset = np.reshape(data_subset, [1, 1, N_samp*N_proj_subset])
coords_subset = np.reshape(coords_subset, [1, N_samp*N_proj_subset, 3])
# dcf = np.reshape(dcf,[1, N_samp*subset,1])

coords_subset = coords_subset[0, ~np.isnan(data_subset[0, 0, :]), :]
# dcf = dcf[0,~np.isnan(data_subset[0,:,0]),0]
data_subset = data_subset[0, 0, ~np.isnan(data_subset[0, 0, :])]

data_subset = np.expand_dims(data_subset, 0)
data_subset = np.expand_dims(data_subset, 1)
coords_subset = np.expand_dims(coords_subset, 0)

# %% Calculate density compensation

try:
    if analytical:
        dcf = ((coords_subset*crd_osamp)[..., 0]**2 +
               (coords_subset*crd_osamp)[..., 1] ** 2 +
               (coords_subset*crd_osamp)[..., 2]**2)**0.5
        dcf = dcf**2  # Square again so in same units as Pipe-menon DCF
    else:
        dcf = mr.pipe_menon_dcf(coords_subset*crd_osamp, img_shape=(np.dot(img_osamp, image_shape).astype(int)),
                                beta=beta, width=width, device=device,
                                max_iter=30)
        dcf /= np.max(dcf)

except:
    print("Error: DCF was not calculated. ")
    raise

# Visualize
fig = plt.figure()
plt.plot(sp.to_device(dcf[0, :N_samp], -1), 'b')
plt.xlabel('Sample number')
plt.ylabel('DCF')
plt.suptitle('Sample Density Compensation Function')
plt.title(r'$ \beta $ = %.2g, width = %.2g' % (beta, width))

# Move to CPU
dcf = mvc(dcf)


# %% Generate a hyperpolarized weighting for the inverse NUFFT

# Simulate T2* decay
t2_star_xe = 18  # ms
echo_time = 0.45  # ms
if "data/floret-kumc/" in folder:
    readout = 10  # ms
else:
    readout = 5 + 0.8  # ms
dwell_time = readout/N_samp
relaxation = np.zeros((N_samp,))
for i in range(N_samp):
    relaxation[i] = np.exp(-(i*dwell_time + echo_time)/t2_star_xe)

# Estimate magnetization decay due to RF excitations
k = np.zeros((N_proj_subset, N_samp))
for i in range(N_proj_subset):
    k[i, :] = k0_subset[i] * relaxation

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

# Reshape
k = np.reshape(k, [1, 1, N_samp*N_proj_subset])

# %% Estimate sensitivity map

sens_map = np.ones((1, image_size, image_size, image_size), dtype=int)
print("Sensitivity map of all ones was assumed.")

# %% Initialize linear algebra operators for the compressed sense formulation

# Set device
with device:

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
    K = sp.linop.Multiply(F.oshape, k**gamma)

    if use_dcf:
        A = K * D * F * S
    else:
        A = K * F * S

    # Calculate maximum eigenvalue
    LL = sp.app.MaxEig(A.N, dtype=xp.complex64,
                       device=device).run() * 1.01

    # Divide by svd eigenvalue to ensure maximum eigenvalue falls between [0,1]
    A = np.sqrt(1/LL) * A

    # Define regularizing linear operators and their proximal operators
    W = sp.linop.Wavelet(S.ishape, wave_name="db4")
    def g1(x): return lamda_1 * xp.linalg.norm(W(x).ravel(), ord=1)
    prox_g1 = sp.prox.UnitaryTransform(
        sp.prox.L1Reg(W.oshape, lamda_1), W)

    prox_g2 = tv.ProxTV(A.ishape, lamda_2)
    def g2(x): return lamda_2 * xp.linalg.norm(prox_g2.G(x))

    # Make list of objectives and proximal operators
    lst_g = [g1, g2]
    lst_proxg = [prox_g1, prox_g2]

# %% Perform image reconstructions

# Set device
with device:

    # Compute norm
    if use_dcf:
        data_norm = xp.linalg.norm(data * dcf**0.5)
    else:
        data_norm = xp.linalg.norm(data)

    # Correct for sample density and normalize
    if use_dcf:
        b = data * dcf**0.5
        b = b / data_norm
    else:
        b = data
        b = b / data_norm

    # Compute a single x = A.H b operation (i.e. inverse NUFFT)
    A_dcf = D * F * S
    b_dcf = data * dcf**0.5 / xp.linalg.norm(data * dcf**0.5)
    img_nufft_gas = mvc(A_dcf.H * b_dcf)
    pl.ImagePlot(np.rot90(img_nufft_gas[:, :, int(recon_resolution/3):int(recon_resolution/1.5):2], 2),
                 x=0, y=1, z=2,
                 title="Reconstruction using inverse NUFFT")

    # Make a point spread function to understand aliasing
    img_nufft_psf = mvc(A.H * xp.ones(np.shape(b)))

    # Solve ||Ax - b||^2 + lamda||x||^2 using conjugate gradient descent
    img_cg_gas = mvc(convexalg.cg(num_iters=num_iters,  # num_iters,
                                  ptol=ptol,
                                  A=A, b=b,
                                  lamda=lamda_l2,  # l2-norm constraint for additional regularization
                                  verbose=True,
                                  draw_output=True,
                                  save=None))
    pl.ImagePlot(np.rot90(img_cg_gas[:, :, int(recon_resolution/3):int(recon_resolution/1.5):2], 2),
                 x=0, y=1, z=2,
                 title="Reconstruction using CG")

    # Solve ||Ax - b||^2 + lamda_1 |Wx| + lamda_2 |Gx| using ADMM
    img_cs_cg_gas = mvc(convexalg.admm(num_iters=num_iters,
                                       ptol=ptol,
                                       A=A, b=b,
                                       num_normal=num_normal,
                                       lst_proxg=lst_proxg,
                                       rho=rho,
                                       lst_g=lst_g,
                                       method="cg",
                                       verbose=True,
                                       draw_output=True))
    pl.ImagePlot(np.rot90(img_cs_cg_gas[:, :, int(recon_resolution/3):int(recon_resolution/1.5):2], 2),
                 x=0, y=1, z=2,
                 title='Reconstruction using CS: ADMM w/ CG')

    LL_gas = LL
    norm_dcf_gas = xp.linalg.norm(data * dcf**0.5)
    norm_gas = xp.linalg.norm(data)

# %% Now set up the dissolved phase recon

data_file = "ksp_d.npy"
coords_file = "coord_d.npy"


# Load data (channels, projections, samples)
data = np.array(np.load(folder + data_file))

# Visualize
fig, (ax1, ax2) = plt.subplots(2, sharex=False)
ax1.plot(abs(data[0, :, 0]), color='r')
ax1.set_ylabel('$k_0$ intensity')
ax1.set_title("magnetization decay")
ax1.set_xlabel("projection number")
ax2.plot(abs(data[0, 0, :]), color='b')
ax2.plot(abs(data[0, 50, :]), color='b')
ax2.plot(abs(data[0, 100, :]), color='b')
ax2.set_ylabel('k-space intensity')
ax2.set_xlabel('sample number')
ax2.set_title("fid (three arbitrary projections)")
fig.tight_layout(pad=1.0, w_pad=0.5, h_pad=1.0)

# Normalize data to max k0
k0 = abs(data[0, :, 0])
data_normalized = np.zeros_like(data)
for i in range(np.shape(data)[1]):
    data_normalized[0, i, :] = data[0, i, :] / k0[i]

if normalize_data:
    data = data_normalized

# Load coords (projections, samples, dimensions)
coords = np.array(np.load(folder + coords_file))

# Edit dimensions of array for dissolved phase (radial only)
if "data/radial-" or "paper" in folder:
    coords[:, :, [2, 1]] = coords[:, :, [1, 2]]

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
coords_resize = image_size * coords / resize_factor

# %% Subset the data to experiment with undersampling

# Steady state
steady_state = 0

# Select only a subset of trajectories and data
if subset_type == "first n":
    subset = int(N_proj * subsample_pc_diss / 100)
    data_subset = data[:, steady_state:subset, :]
    coords_subset = coords_resize[steady_state:subset, ...]
else:
    subset = int(100 / subsample_pc_diss)
    data_subset = data[:, steady_state::subset, :]
    coords_subset = coords_resize[steady_state::subset, ...]

# Update number of projections
N_proj_subset = data_subset.shape[1]

# Update k0
k0_subset = abs(data_subset[0, :, 0])

# %% Deal with NaN values
data_subset = np.reshape(data_subset, [1, 1, N_samp*N_proj_subset])
coords_subset = np.reshape(coords_subset, [1, N_samp*N_proj_subset, 3])
# dcf = np.reshape(dcf,[1, N_samp*subset,1])

nan_idx = ~np.isnan(data_subset[0, 0, :])
coords_subset = coords_subset[0, nan_idx, :]
# dcf = dcf[0,~np.isnan(data_subset[0,:,0]),0]
data_subset = data_subset[0, 0, nan_idx]

data_subset = np.expand_dims(data_subset, 0)
data_subset = np.expand_dims(data_subset, 1)
coords_subset = np.expand_dims(coords_subset, 0)

# %% Calculate density compensation

try:
    if analytical:
        dcf = ((coords_subset*crd_osamp)[..., 0]**2 +
               (coords_subset*crd_osamp)[..., 1] ** 2 +
               (coords_subset*crd_osamp)[..., 2]**2)**0.5
        dcf = dcf**2  # Square again so in same units as Pipe-menon DCF
    else:
        dcf = mr.pipe_menon_dcf(coords_subset*crd_osamp, img_shape=(np.dot(img_osamp, image_shape).astype(int)),
                                beta=beta, width=width, device=device,
                                max_iter=30)
        dcf /= np.max(dcf)
except:
    print("Error: DCF was not calculated. ")
    raise

# Visualize
fig = plt.figure()
plt.plot(sp.to_device(dcf[0, :N_samp], -1), 'b')
plt.xlabel('Sample number')
plt.ylabel('DCF')
plt.suptitle('Sample Density Compensation Function for Dissolved Phase')
plt.title(r'$ \beta $ = %.2g, width = %.2g' % (beta, width))

# Move to CPU
dcf = mvc(dcf)

# %% Generate a hyperpolarized weighting for the inverse NUFFT

# Simulate T2* decay
t2_star_xe = 1  # ms
readout = 0.5*0.8  # ms
dwell_time = readout/N_samp
echo_time = 0.45 + (steady_state*dwell_time)  # ms
relaxation = np.zeros((N_samp,))
for i in range(N_samp):
    relaxation[i] = np.exp(-(i*dwell_time + echo_time)/t2_star_xe)

# Estimate magnetization decay due to RF excitations
k0_subset = k0_subset[~np.isnan(k0_subset)]
k = np.zeros((N_proj_subset, N_samp))
for i in range(N_proj_subset):
    k[i, :] = k0_subset[i] * relaxation

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

k = np.reshape(k, [1, 1, N_samp*N_proj_subset])
k = k[0, 0, nan_idx]


# %% Initialize linear algebra operators for the compressed sense formulation

# Set device
with device:

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
                       width=2,  # Recommended 2-4
                       toeplitz=True)
    D = sp.linop.Multiply(F.oshape, dcf**0.5)
    K = sp.linop.Multiply(F.oshape, k**gamma)

    if use_dcf:
        A = K * D * F * S
    else:
        A = K * F * S

    # Calculate maximum eigenvalue
    LL = sp.app.MaxEig(A.N, dtype=xp.complex64,
                       device=device).run() * 1.01

    # Divide by svd eigenvalue to ensure maximum eigenvalue falls between [0,1]
    A = np.sqrt(1/LL) * A

    # Make list of objectives and proximal operators
    # lst_g = [g3]
    # lst_proxg = [prox_g3]
    # lst_g = [g4]
    # lst_proxg = [prox_g4]

    # Compute norm
    if use_dcf:
        data_norm = xp.linalg.norm(data * dcf**0.5)
    else:
        data_norm = xp.linalg.norm(data)

    # Correct for sample density and normalize
    if use_dcf:
        b = data * dcf**0.5
        b = b / data_norm
    else:
        b = data
        b = b / data_norm

    # Compute a single x = A.H b operation (i.e. inverse NUFFT)
    A_dcf = D * F * S
    b_dcf = data * dcf**0.5 / xp.linalg.norm(data * dcf**0.5)
    img_nufft_dissolved = mvc(A_dcf.H * b_dcf)
    pl.ImagePlot(np.rot90(img_nufft_dissolved[:, :, int(recon_resolution/3):int(recon_resolution/1.5):2], 2),
                 x=0, y=1, z=2,
                 title="Reconstruction using inverse NUFFT")

    # Solve ||Ax - b||^2 + lamda||x||^2 using conjugate gradient descent
    img_cg_dissolved = mvc(convexalg.cg(num_iters=num_iters,
                                        ptol=ptol,
                                        A=A, b=b,
                                        lamda=lamda_l2,  # l2-norm constraint for additional regularization
                                        verbose=True,
                                        draw_output=True,
                                        save=None))
    pl.ImagePlot(np.rot90(img_cg_dissolved[:, :, int(recon_resolution/3):int(recon_resolution/1.5):2], 2),
                 x=0, y=1, z=2,
                 title="Reconstruction using CG")
    pl.ImagePlot(np.rot90(np.transpose(img_cg_dissolved[:, :, :]), axes=(0, 2), k=-1),
                 x=0, y=2, z=1,
                 title='Reconstruction using CG')

    # Solve ||Ax - b||^2 + lamda_1 |Wx| + lamda_2 |Gx| using ADMM
    img_cs_cg_dissolved = mvc(convexalg.admm(num_iters=num_iters,
                                             ptol=ptol,
                                             A=A, b=b,
                                             num_normal=num_normal,
                                             lst_proxg=lst_proxg,
                                             rho=rho,
                                             lst_g=lst_g,
                                             method="cg",
                                             verbose=True,
                                             draw_output=True))
    pl.ImagePlot(np.rot90(img_cs_cg_dissolved[:, :, int(recon_resolution/3):int(recon_resolution/1.5):2], 2),
                 x=0, y=1, z=2,
                 title='Reconstruction using CS: ADMM w/ CG')
    pl.ImagePlot(np.rot90(np.transpose(img_cs_cg_dissolved[:, :, :]), axes=(0, 2), k=-1),
                 x=0, y=2, z=1,
                 title='Reconstruction using CS: ADMM w/ CG')

    LL_dissolved = LL
    norm_dcf_dissolved = xp.linalg.norm(data * dcf**0.5)
    norm_dissolved = xp.linalg.norm(data)


# %% Crop edges of image

crop = True

if crop:
    # Pad width (some radial scans have noise collected at edges due to undersampling)
    pad = pad_width  # pixels

    img_nufft_gas = img_nufft_gas[pad:-pad, pad:-pad, pad:-pad]
    img_cg_gas = img_cg_gas[pad:-pad, pad:-pad, pad:-pad]
    img_cs_cg_gas = img_cs_cg_gas[pad:-pad, pad:-pad, pad:-pad]
    img_nufft_dissolved = img_nufft_dissolved[pad:-pad, pad:-pad, pad:-pad]
    img_cg_dissolved = img_cg_dissolved[pad:-pad, pad:-pad, pad:-pad]
    img_cs_cg_dissolved = img_cs_cg_dissolved[pad:-pad, pad:-pad, pad:-pad]

    # Update recon_resolution
    recon_resolution = img_nufft_gas.shape[0]
    G = sp.linop.FiniteDifference(img_nufft_gas.shape)

    print("Images were cropped by " + str(pad) + " pixels from the edges.")


# %% Image segmentation and noise estimation (automatic)

img = img_nufft_gas

# Perform clustering
N_clusters = 2
kmeans = KMeans(n_clusters=N_clusters, random_state=0).fit(np.reshape(
    # Reshape to (voxelcount,1) shape array
    abs(img), (np.shape(np.ravel(img))[0], 1)))
clustered = kmeans.cluster_centers_[kmeans.labels_]

# Binarize
if N_clusters == 2:
    clustered[clustered == np.percentile(clustered, 75)] = 1
    clustered[clustered == np.percentile(clustered, 25)] = 0

# Reshape clusters into images
img_cluster = np.zeros(shape=(N_clusters, np.shape(
    img)[0], np.shape(img)[1], np.shape(img)[2]))
labels = kmeans.labels_
for n in range(N_clusters):
    cluster_temp = []
    for i in range(len(labels)):
        if (labels[i]) == n:
            cluster_temp.append(float(clustered[i]))
        else:
            cluster_temp.append(1)
    img_cluster[n, :, :, :] = np.array(cluster_temp).reshape(img.shape)

# Calculate norm of each cluster*img_cs to see which one is larger (and which one is signal versus noise)
norm_cluster0 = np.linalg.norm(img_cluster[0, :, :, :]*abs(img))
norm_cluster1 = np.linalg.norm(img_cluster[1, :, :, :]*abs(img))

if norm_cluster0 >= norm_cluster1:
    img_signal = img_cluster[0, :, :, :]
    img_noise = img_cluster[1, :, :, :]
else:
    img_signal = img_cluster[1, :, :, :]
    img_noise = img_cluster[0, :, :, :]

# Binarize noise mask
img_noise[img_noise != 1] = 0

# Erode the noise image
kernel = np.ones((7, 7), np.uint8)
img_noise = cv2.erode(img_noise, kernel, iterations=3)

# Dilate the signal image (optional)
# kernel = np.ones((3, 3), np.uint8)
kernel = np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]], np.uint8)
img_signal = cv2.dilate(img_signal, kernel, iterations=1)

# Visualize
sp.plot.ImagePlot(img_signal + 1 - img_noise, x=0, y=1, z=2,
                  title="Signal and noise mask, estimated using k-means",
                  vmin=0, vmax=2, colormap="YlGnBu_r")

# %% Make a beautiful plot showing magnitude, SNR, and phase

# Select slice
plot_slice = int(recon_resolution*0.53)

# Translation and rotation
dx = 0
dy = 0
rot = 1


img = np.flipud(util.ps(util.ir(img_cs_cg_gas, rot), dx, dy))
img = cv2.normalize(abs(img), None, 0, 1.0, cv2.NORM_MINMAX, dtype=cv2.CV_32F)

# Generate images
img_slice = np.fliplr(np.rot90(abs(img[plot_slice, :, :]), 3))
img_slice /= np.percentile(np.ravel(img_slice), 99.9)
img_snr_slice = np.fliplr(np.rot90(util.calculate_snr(
    abs(img[plot_slice, :, :]), abs(img[img_noise == 1])), 3))
img_grad_slice = np.fliplr(
    np.rot90(abs(G.N*img)[plot_slice, :, :], 3))
# img_grad_slice = np.fliplr(
#     np.rot90((eta/(np.sqrt((eta**2 + abs(G.N*img)**2))))[plot_slice, :, :], 3))
# img_weights_slice = np.fliplr(
#     np.rot90(abs(np.flipud(util.ps(util.ir(mvc(weights), rot), dx, dy)))[plot_slice, :, :], 3))

# Plot
plt.figure(figsize=(10, 30), dpi=300)
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, sharey=True, gridspec_kw={'wspace': 0, 'hspace': 0},
                                    squeeze=True)
img_plot1 = ax1.imshow(img_slice,
                       cmap="gray", vmin=0, vmax=1)
# img_plot2 = ax2.imshow(img_snr_slice,
#                        cmap="jet", vmin=0, vmax=20)
img_plot2 = ax2.imshow(np.log(img_slice),
                       cmap="gnuplot", vmin=-3, vmax=0)
img_plot3 = ax3.imshow(img_grad_slice,
                       cmap="gnuplot", vmin=0, vmax=1)
ax1.axis('off')
ax2.axis('off')
ax3.axis('off')

plt.figure(figsize=(10, 20), dpi=300)
fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True, gridspec_kw={'wspace': 0, 'hspace': 0},
                               squeeze=True)
img_plot1 = ax1.imshow(img_slice,
                       cmap="gray", vmin=0, vmax=1)
img_plot2 = ax2.imshow(img_grad_slice,
                       cmap="gnuplot", vmin=0, vmax=1)
ax1.axis('off')
ax2.axis('off')

# cbar = fig.colorbar(img_plot3, location='bottom')
# cbar.outline.set_edgecolor('black')
# cbar.set_ticks([0, 0.5, 1])
# cbar.set_ticklabels(["0.0", "0.5", "1.0"])
# cbar.ax.tick_params(labelsize=16)
# plt.figure(figsize=(10, 10), dpi=300)
# im = plt.imshow(img_weights_slice,
#                 cmap="gnuplot", vmin=0, vmax=1)
# plt.axis('off')

img = np.flipud(util.ps(util.ir(img_cs_cg_gas, rot), dx, dy))
# img = np.flipud(util.ps(util.ir(img_nufft_gas, rot), dx, dy))

# Generate images
img_slice = np.fliplr(np.rot90(abs(img[plot_slice, :]), 3))
img_slice /= np.percentile(np.ravel(img_slice), 99.9)
img_snr_slice = np.fliplr(np.rot90(util.calculate_snr(
    abs(img[plot_slice, :, :]), abs(img[img_noise == 1])), 3))
img_phase_slice = np.fliplr(
    np.rot90(np.angle(img[plot_slice, :, :]), 3))

# Plot
plt.figure(figsize=(10, 30), dpi=300)
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, sharey=True, gridspec_kw={'wspace': 0, 'hspace': 0},
                                    squeeze=True)
img_plot1 = ax1.imshow(img_slice,
                       cmap="gray", vmin=0, vmax=1)
img_plot2 = ax2.imshow(img_snr_slice,
                       cmap="jet", vmin=0, vmax=20)
img_plot3 = ax3.imshow(img_phase_slice,
                       cmap="twilight", vmin=-3.14, vmax=3.14)
ax1.axis('off')
ax2.axis('off')
ax3.axis('off')

plt.figure(figsize=(10, 30), dpi=300)
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, sharey=True, gridspec_kw={'wspace': 0, 'hspace': 0},
                                    squeeze=True)
img_plot1 = ax1.imshow(img_slice,
                       cmap="gray", vmin=0, vmax=1)
img_plot2 = ax2.imshow(np.log(img_slice),
                       cmap="gnuplot", vmin=-3, vmax=0)
img_plot3 = ax3.imshow(img_phase_slice,
                       cmap="twilight", vmin=-3.14, vmax=3.14)
ax1.axis('off')
ax2.axis('off')
ax3.axis('off')


plt.figure(figsize=(10, 10), dpi=300)
fig, (ax1) = plt.subplots(1, 1, sharey=True, gridspec_kw={'wspace': 0, 'hspace': 0},
                          squeeze=True)
img_plot1 = ax1.imshow(img_slice,
                       cmap="gray", vmin=0, vmax=1)
# img_plot1 = ax1.imshow(img_grad_slice,
#                        cmap="gnuplot", vmin=0, vmax=1)
# img_plot1 = ax1.imshow(img_phase_slice,
#                        cmap="twilight", vmin=-3.14, vmax=3.14)
ax1.axis('off')


# %% Generate a ratio image

plot_slice = plot_slice
img = np.flipud(util.ps(util.ir(img_nufft_dissolved/img_nufft_gas)))
corr_factor = mvc((norm_dcf_dissolved)/(norm_dcf_gas))
img_slice = np.fliplr(np.rot90(abs(img[plot_slice, :, :]), 3))

plt.figure(figsize=(10, 10), dpi=300)
fig, (ax1) = plt.subplots(1, 1, sharey=True, gridspec_kw={'wspace': 0, 'hspace': 0},
                          squeeze=True)
img_plot1 = ax1.imshow(img_slice * corr_factor,
                       cmap="twilight_shifted", vmin=0, vmax=1)
# fig.colorbar(img_plot1)
ax1.axis('off')

img = np.flipud(util.ps(util.ir(img_cs_cg_dissolved/img_cs_cg_gas)))
corr_factor = mvc((norm_dissolved)/(norm_gas))
img_slice = np.fliplr(np.rot90(abs(img[plot_slice, :, :]), 3))

plt.figure(figsize=(10, 10), dpi=300)
fig, (ax1) = plt.subplots(1, 1, sharey=True, gridspec_kw={'wspace': 0, 'hspace': 0},
                          squeeze=True)
img_plot1 = ax1.imshow(img_slice * corr_factor,
                       cmap="twilight_shifted", vmin=0, vmax=1)
# fig.colorbar(img_plot1)
ax1.axis('off')


# %% Optional: Apply an artificial filter to the image to remove undersampling artifacts at edges

fade = False

if fade:
    def create_3d_fading_circle(size, radius):
        # Create a grid of coordinates
        x, y, z = np.meshgrid(np.arange(size[1]), np.arange(
            size[0]), np.arange(size[2]))

        # Calculate the distance from the center
        distance = np.sqrt((x - size[1] / 2)**2 +
                           (y - size[0] / 2)**2 + (z - size[2] / 2)**2)

        # Create a fading effect using the distance
        fading_effect = 1 - \
            np.clip((distance - radius) /
                    (size[1] / 2 - radius), 1e-9, 1 - (1e-9))

        # Create the 3D fading circle
        fading_circle_3d = np.ones(size) * fading_effect

        return 1 - fading_circle_3d

    # Set the size of the 3D array and the radius of the circle
    array_size_3d = image_shape
    circle_radius_3d = recon_resolution//2 + scan_resolution//3

    # Create the fading circle
    fading_circle = create_3d_fading_circle(array_size_3d, circle_radius_3d)

    sp.plot.ImagePlot(fading_circle, x=0, y=1, z=2,
                      title="Artificial filter", colormap='gray')

    def filter(x): return fading_circle*x

    img_nufft_gas = filter(abs(img_nufft_gas))
    img_cg_gas = filter(abs(img_cg_gas))
    img_cs_cg_gas = filter(abs(img_cs_cg_gas))
    img_nufft_dissolved = filter(abs(img_nufft_dissolved))
    img_cg_dissolved = filter(abs(img_cg_dissolved))
    img_cs_cg_dissolved = filter(abs(img_cs_cg_dissolved))
else:
    print("Edges of images were not faded automatically.")


# %% Make an image summary


# Select slice
plot_slice = int(recon_resolution*0.55)

# Set up subplot
plt.figure()
fig, axes = plt.subplots(nrows=3, ncols=2, figsize=(6, 9))
if use_dcf:
    fig.suptitle(folder + '\n max_iter=%i, \u03C1=%.1e, $p_{tol}$=%.1e, $A^TA$=%i,  \n \u03BB$_1$=%.1e, \u03BB$_2$=%.1e, \u03BB$_3$=%.1e, \zeta=%.1f, \n \u03BB$_4$=%.1e, $\eta$=%.1e' % (
        num_iters, rho, ptol, num_normal, lamda_1, lamda_2, lamda_3, zeta, lamda_4, eta), y=1.03)
else:
    fig.suptitle(folder + '\nmax_iter=%i, \u03C1=%.1e, $p_{tol}$=%.1e, $A^TA$=%i, \n \u03BB$_1$=%.1e, \u03BB$_2$=%.1e, \u03BB$_3$=%.1e, \zeta=%.1f, \n \u03BB$_4$=%.1e, $\eta$=%.1e' % (
        num_iters, rho, ptol, num_normal, lamda_1, lamda_2, lamda_3, zeta, lamda_4, eta), y=1.03)

ax = axes[0, 0]
im = ax.imshow(
    abs(np.flipud(np.rot90(img_nufft_gas[:, plot_slice, :]))), cmap="gray")
ax.set_title(r'$F^HD(b)$')
fig.colorbar(im, ax=ax)
ax.axis('off')
ax = axes[0, 1]
im = ax.imshow(
    abs(np.flipud(np.rot90(img_nufft_dissolved[:, plot_slice, :]))), cmap="gray")
ax.set_title(r'$F^HD(b)$')
fig.colorbar(im, ax=ax)
ax.axis('off')

ax = axes[1, 0]
im = ax.imshow(
    abs(np.flipud(np.rot90(img_cg_gas[:, plot_slice, :]))), cmap="gray")
ax.set_title(r'$||Ax-b||^2_2$')
fig.colorbar(im, ax=ax)
ax.axis('off')
ax = axes[1, 1]
im = ax.imshow(
    abs(np.flipud(np.rot90(img_cg_dissolved[:, plot_slice, :]))), cmap="gray")
ax.set_title(r'$||Ax-b||^2_2$')
fig.colorbar(im, ax=ax)
ax.axis('off')

ax = axes[2, 0]
im = ax.imshow(
    abs(np.flipud(np.rot90(img_cs_cg_gas[:, plot_slice, :]))), cmap="gray")
ax.set_title(r'$||Ax-b||^2_2 + g(x)$ (cg)')
fig.colorbar(im, ax=ax)
ax.axis('off')
ax = axes[2, 1]
im = ax.imshow(
    abs(np.flipud(np.rot90(img_cs_cg_dissolved[:, plot_slice, :]))), cmap="gray")
ax.set_title(r'$||Ax-b||^2_2 + g(x)$ (cg)')
fig.colorbar(im, ax=ax)
ax.axis('off')
