# -*- coding: utf-8 -*-
"""
Title: 
    ventilation_floret_3d_single_channel.py

Agenda: 
    A .py script that uses compressed sensing methods to reconstruct 
non-cartesian Philips data in raw-lab-sin/data-list format. 

Author: 
    Joseph Plummer - joseph.plummer@cchmc.org 

Creation date: 
    2022-09-19
    
Modification date: 
    2026-04-30 - ASB initial commit in XIPline
    
"""

# %% Import packages

from scipy.io import savemat
from skimage.metrics import structural_similarity as ssim
import cv2
from sklearn.cluster import KMeans
import os
import nibabel as nib
from functions import tv, convexalg, util
import sigpy.mri as mr
import sigpy.plot as pl
import sigpy as sp
import ReadPhilips.readphilips as rp
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('dark_background')
#import sys
#custom_path = r"C:\Users\bda5ik\AppData\Local\anaconda3\envs\recon\Lib\site-packages"
#if custom_path not in sys.path:
#    sys.path.insert(0, custom_path)

# % Reconstruction settings
# Device settings
devnum = 0
device = sp.Device(devnum)
xp = device.xp

# Density compensation
use_dcf = True

# Subsample
subsample_pc = 100  # percentage

# Normalize data
normalize_data = False

# Multiply F by signal weighting
gamma = 1  # 1 = full HP weighting on F, 0 = no weighting on F

# Proximal gradient descent settings
if use_dcf:
    num_iters = 50
    # Regularization term for wavelet energy (1e-6~1e-2, scales logarithmically)
    lamda_1 = 8e-5
    # Regularization term for total variation energy (1e-6~1e-2, scales logarithmically)
    lamda_2 = 1e-2
    lamda_l2 = 1e-5
    # ADMM step size
    rho = 1
    ptol = 1e-1
    num_normal = 4
else:
    # More iterations required without left-preconditioning (i.e. DCF)
    num_iters = 1000
    lamda_1 = 3e-8
    lamda_2 = 3e-8
    lamda_l2 = 1e-5
    rho = 5e1
    rho = 1e3  # 186H-1044
    ptol = 1e-2
    num_normal = 11


# %% Import data and trajectories

def read_config_and_find_files(config_path=r"C:\XIPline\offline_recon\config_Xe.txt"):
    # Read the config file
    with open(config_path, 'r') as f:
        lines = f.readlines()

    # Extract datalocation and freq_offset
    datalocation = None
    freq_offset = None
    for line in lines:
        line = line.strip()
        if line.lower().startswith("datalocation"):
            datalocation = line.split('=', 1)[1].strip()
        elif line.lower().startswith("freq_offset"):
            freq_offset = float(line.split('=', 1)[1].strip())

    if datalocation is None or freq_offset is None:
        raise ValueError("Missing datalocation or freq_offset in config file.")

    # Find .data and .sin files in the datalocation
    datafile = trajfile = None
    for file in os.listdir(datalocation):
        full_path = os.path.join(datalocation, file)
        if file.endswith(".raw"):
            datafile = full_path
        elif file.endswith(".sin"):
            trajfile = full_path

    if datafile is None or trajfile is None:
        raise FileNotFoundError("Missing .data or .sin file in datalocation.")

    # Return as raw string literals
    datafile_r = fr"{datafile}"
    trajfile_r = fr"{trajfile}"
    datalocation_r = fr"{datalocation}"

    return datafile_r, trajfile_r, datalocation_r, freq_offset

config_Xe_path=r"C:\XIPline\offline_recon\config_Xe.txt"
datafile, trajfile, datalocation, freq_offset = read_config_and_find_files(config_Xe_path)
print("Data file:", datafile)
print("Trajectory file:", trajfile)
print("Data location:", datalocation)
print("Frequency offset:", freq_offset)


#%%
#location = r"D:\OneDrive - cchmc\Lab\Random Subject analysis\CPIR_FLORET_Ventilation"
# Exact filenames in that folder
#fname_coord = r"20251203_105839_CPIR_VENT_FLORET_303750_285140_13_2.sin"
#fname_data  = r"raw_211.data"
#fname_data  = r"20251203_105839_CPIR_VENT_FLORET_303750_285140_13_2.raw"

location = datalocation
fname_coord = trajfile
fname_data  = datafile

# Load data
filename_data = os.path.join(location, fname_data)
dl = rp.PhilipsData(filename_data)
dl.raw_corr = True
dl.compute()
data_load = np.array(dl.data)
data_load = np.squeeze(data_load)
data_load = np.reshape(data_load, (np.shape(data_load)[
                        0]*np.shape(data_load)[1], np.shape(data_load)[2]))
if np.size(np.shape(data_load)) != 2:
    print('Incorrect data dimensions.')

# Load traj
filename_coords = os.path.join(location, fname_coord)
traj = rp.PhilipsData(filename_coords)
traj.compute()

coords_init = traj.spparams['COORDS_EXPANDED']
coords = np.transpose(np.squeeze(coords_init))
if np.size(np.shape(coords)) == 4:
    print('Incorrect coords dimensions - will perform reshape.')
    coords = np.reshape(coords, (np.shape(
        coords)[0]*np.shape(coords)[1], np.shape(coords)[2], np.shape(coords)[3]))

print(coords.shape)

# ---- Load DCF (SDC.npy) ----
''' 
try:
    dcf_load = np.load(os.path.join(location, fname_dcf))
    dcf_load = np.squeeze(dcf_load)
    if dcf_load.ndim == 3:
        print('Incorrect dcf dimensions - will perform reshape.')
        dcf_load = np.reshape(
            dcf_load,
            (dcf_load.shape[0] * dcf_load.shape[1], dcf_load.shape[2])
        )
except FileNotFoundError:
    print("Could not load DCF.")
'''

#%% Understand data sizes and merge dimensions
N_dimensions = np.shape(coords)[2]  # dimensions
N_samp = np.shape(data_load)[1]  # samples per projection
N_proj = np.shape(data_load)[0]  # number of projections

# Plot trajectories:
plt.figure()
ax = plt.axes(projection='3d')
N_visual = 100  # Number of projections you want to show, inefficient for large n
color = iter(plt.cm.viridis(np.linspace(0, 1, N_visual)))
for i in np.linspace(0, N_proj-1, N_visual):
    i = int(i)
    c = next(color)
    ax.scatter(coords[i, :, 0], coords[i, :, 1],
               coords[i, :, 2], color=c, s=0.5, marker='.')
ax.set_zlabel('$k_z$')
ax.set_ylabel('$k_y$')
ax.set_xlabel('$k_x$')
plt.title('FLORET coordinates')
plt.show()


# Define some parameters manually
recon_resolution = 70  # 120 for ISMRM, 128 for Penn
scan_resolution = 70  # 90 for ISMRM, 100 for Penn

# %% Examine k-space

# Remove first N points if corrupted (e.g. k0 not equal to high value)
# data_load = data_load[:, 4:]
# dcf_load = dcf_load[:, 4:]
# coords = coords[:, 4:, :]
# np.save("data.npy", data_load)
# np.save("dcf.npy", dcf_load)
# np.save("coords.npy", coords)

# Visualize
fig, (ax1, ax2) = plt.subplots(2, sharex=False)
ax1.plot(abs(data_load[:, 0]), color='r')
ax1.set_ylabel('k0 intensity')
ax1.set_title("magnetization decay")
ax1.set_xlabel("projection number")
ax2.plot(abs(data_load[0, :]), color='b')
ax2.set_ylabel('k-space intensity')
ax2.set_xlabel('sample number')
ax2.set_title("fid")
fig.tight_layout(pad=1.0, w_pad=0.5, h_pad=1.0)

# Normalize data to max k0
k0 = np.abs(data_load[:, 0])

data_normalized = np.zeros_like(data_load)

for i in range(data_load.shape[0]):  # loop over rows
    data_normalized[i, :] = data_load[i, :] / k0[i]

if normalize_data:
    data_load = data_normalized

# %% Make manipulations to reconstructed volume

alpha = 1  # reconstruct onto larger volume
ImageSize = scan_resolution / alpha

resize_factor = recon_resolution/(scan_resolution)
ImageSize = int(resize_factor * ImageSize)
ImageShape = (int(ImageSize), int(ImageSize), int(ImageSize))
traj = coords / resize_factor
traj = ImageSize * traj

# Resize for visualization
# traj = traj / 0.8


# %% Subset the data

# Select only a subset of trajectories and data
subset = int(N_proj*subsample_pc/100)
data_load_subset = data_load[:subset, :]
traj_subset = traj[:subset, ...]

# Update values
N_samp = np.shape(data_load_subset)[1]  # samples per projection
N_proj = np.shape(data_load_subset)[0]  # number of projections

# Plot trajectories:
plt.figure()
ax = plt.axes(projection='3d')
N_visual = N_proj//20  # Number of projections you want to show, inefficient for large n
color = iter(plt.cm.viridis(np.linspace(0, 1, N_visual)))
for i in np.linspace(0, N_proj-1, N_visual):
    i = int(i)
    c = next(color)
    ax.scatter(traj_subset[i, :, 0], traj_subset[i, :, 1],
               traj_subset[i, :, 2], color=c, s=0.5, marker='.')
# ax.axis('off')
color_tuple = (0.0, 0.0, 0.0, 0.0)

# Remove pane background
ax.xaxis.set_pane_color(color_tuple)
ax.yaxis.set_pane_color(color_tuple)
ax.zaxis.set_pane_color(color_tuple)

# Remove axis lines
ax.xaxis.line.set_color(color_tuple)
ax.yaxis.line.set_color(color_tuple)
ax.zaxis.line.set_color(color_tuple)

# Labels
ax.set_xlabel(r'$k_x$', fontsize=18)
ax.set_ylabel(r'$k_y$', fontsize=18)
ax.set_zlabel(r'$k_z$', fontsize=18)

# Remove tick labels
ax.set_xticklabels([])
ax.set_yticklabels([])
ax.set_zticklabels([])

plt.show()

# %% Manipulate data size to deal with linear algebra operations
# Akin to adding a channel dimension for sense operations

# Add an extra dimension for channel
data_load_channel = np.reshape(data_load_subset, (1, np.shape(
    data_load_subset)[0], np.shape(data_load_subset)[1]))

# Estimate sensitivity map - assume all ones for single channel
mps = np.ones((1, ImageSize, ImageSize, ImageSize), dtype=int)

# %% Calculate density compensation

beta = 8  # default = 8
width = 4  # default = 4

try:
    dcf = abs(dcf_load)
    dcf = dcf[:subset, :]
except:
    # Density compensation
    dcf = mr.pipe_menon_dcf(traj_subset, img_shape=ImageShape,
                            beta=beta, width=width, device=device)

    # Make final dcf point smaller value
    dcf[:, -1] = np.median(dcf, axis=1)

# Visualize
fig = plt.figure()
plt.plot(sp.to_device(dcf[0, ...], -1), 'b')
plt.xlabel('Sample number')
plt.ylabel('DCF')
plt.suptitle('Sample Density Compensation Function')
plt.title(r'$ \beta $ = %.2g, width = %.2g' % (beta, width))
dcf = np.reshape(dcf, (1, np.shape(dcf)[0], np.shape(dcf)[1]))


# %% Generate a hyperpolarized weighting for the inverse NUFFT

# Simulate T2* decay
t2_star_xe = 18  # ms
readout = 15  # ms
dwell_time = readout/N_samp
relaxation = np.zeros((N_samp,))
for i in range(N_samp):
    relaxation[i] = np.exp(-i*dwell_time/t2_star_xe)

# Simulate RF decay
k = np.zeros((N_proj, N_samp))
for i in range(N_proj):
    # With T2star
    k[i, :] = k0[i] * relaxation

    # Without T2star
    # k[i, :] = k0[i]


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

# %% Initialize linear algebra operators

# Set device
with device:
    def mvd(x): return sp.to_device(x, device)
    def mvc(x): return sp.to_device(x, sp.cpu_device)

    # Move to GPU
    traj = mvd(traj_subset)
    dcf = mvd(dcf)
    mps = mvd(mps)

    # Compute linear operators
    S = sp.linop.Multiply(ImageShape, mps)
    F = sp.linop.NUFFT(mps.shape,
                       coord=traj,
                       oversamp=1.25,
                       width=4,
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

    # Define regularizing linear operators
    W1 = sp.linop.Wavelet(S.ishape, wave_name="db4")
    def g1(x): return lamda_1 * xp.linalg.norm(W1(x).ravel(), ord=1)
    prox_g1 = sp.prox.UnitaryTransform(
        sp.prox.L1Reg(W1.oshape, lamda_1), W1)

    prox_g2 = tv.ProxTV(A.ishape, lamda_2)
    def g2(x): return lamda_2 * xp.linalg.norm(prox_g2.G(x))

    lst_g = [g1, g2]
    lst_proxg = [prox_g1, prox_g2]

# %% Perform image reconstructions

# Set device
with device:
    def mvd(x): return sp.to_device(x, device)
    def mvc(x): return sp.to_device(x, sp.cpu_device)

    data = mvd(data_load_channel)

    # Compute norm
    if use_dcf:
        data_norm_dcf = xp.linalg.norm(data * dcf**0.5)
    else:
        data_norm = xp.linalg.norm(data)

    # Correct for sample density and normalize
    if use_dcf:
        b = data * dcf**0.5
        b = b/xp.linalg.norm(data_norm_dcf)
    else:
        b = data
        b = b/xp.linalg.norm(data_norm)

    # Compute a single x = A.H b operation (i.e. inverse NUFFT)
    A_dcf = D * F * S
    b_dcf = data * dcf**0.5 / xp.linalg.norm(data * dcf**0.5)
    img_nufft = mvc(abs(A_dcf.H * b_dcf))
    pl.ImagePlot(np.rot90(img_nufft[:, :, int(recon_resolution/3):int(recon_resolution/1.5):2], 2),
                 x=0, y=1, z=2,
                 title="Reconstruction using inverse NUFFT")
    
    # Make a point spread function to understand aliasing
    b_dcf_psf = xp.ones(np.shape(b)) * dcf**0.5
    img_nufft_psf = mvc(
        abs(A_dcf.H * b_dcf_psf))
    # img_nufft_psf /= np.median(abs(img_nufft_psf))

    # Solve ||Ax - b||^2 + lamda||x||^2 using conjugate gradient descent
    img_cg = mvc(abs(convexalg.cg(num_iters=200,
                                  ptol=1e-2,
                                  A=A, b=b,
                                  lamda=lamda_l2,  # l2-norm constraint for additional regularization
                                  verbose=True,
                                  draw_output=True,
                                  save=None)))
    pl.ImagePlot(np.rot90(img_cg[:, :, int(recon_resolution/3):int(recon_resolution/1.5):2], 2),
                 x=0, y=1, z=2,
                 title="Reconstruction using CG")

    # Solve ||Ax - b||^2 + lamda_1 |Wx| + lamda_2 |Gx| using ADMM
    img_cs_cg = mvc(abs(convexalg.admm(num_iters=num_iters,
                                       ptol=ptol,
                                       A=A, b=b,
                                       num_normal=num_normal,
                                       lst_proxg=lst_proxg,
                                       rho=rho,
                                       lst_g=lst_g,
                                       method="cg",
                                       verbose=True,
                                       draw_output=True)))
    pl.ImagePlot(np.rot90(img_cs_cg[:, :, int(recon_resolution/3):int(recon_resolution/1.5):2], 2),
                 x=0, y=1, z=2,
                 title='Reconstruction using CS: ADMM w/ CG')
    
# %% Manipulate images for plotting

# Translation and rotation
dx = 0
dy = 0
rot = 1

img_nufft = np.flipud(util.ps(util.ir(img_nufft, rot), dx, dy))
img_cg = np.flipud(util.ps(util.ir(img_cg, rot), dx, dy))
img_cs_cg = np.flipud(util.ps(util.ir(img_cs_cg, rot), dx, dy))


# %% Plot a 3D point spread function:

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

# Make data.
X = np.arange(recon_resolution)
Y = np.arange(recon_resolution)
X, Y = np.meshgrid(X, Y)
Z = np.log(img_nufft_psf[:, :, recon_resolution//2])
# Z = (img_nufft_psf[:, :, recon_resolution//2])

# Plot the surface.
surf = ax.plot_surface(X, Y, Z, cmap="seismic",
                       linewidth=0, antialiased=False, vmin=-3, vmax=3)
ax.set_yticklabels([])
ax.set_xticklabels([])
ax.set_zticklabels([])
ax.axis('off')
cbar = plt.colorbar(surf, location='bottom')
cbar.set_ticks([-2, 0, 2])
cbar.ax.tick_params(labelsize=22)
cbar.outline.set_color('black')

# %% Image segmentation and noise estimation (automatic)

# Perform clustering
N_clusters = 2
kmeans = KMeans(n_clusters=N_clusters, random_state=0).fit(np.reshape(
    img_cs_cg, (np.shape(np.ravel(img_cs_cg))[0], 1)))  # Reshape to (voxelcount,1) shape array
clustered = kmeans.cluster_centers_[kmeans.labels_]

# Binarize
if N_clusters == 2:
    clustered[clustered == np.percentile(clustered, 75)] = 1
    clustered[clustered == np.percentile(clustered, 25)] = 0

# Reshape clusters into images
img_cluster = np.zeros(shape=(N_clusters, np.shape(
    img_cs_cg)[0], np.shape(img_cs_cg)[1], np.shape(img_cs_cg)[2]))
labels = kmeans.labels_
for n in range(N_clusters):
    cluster_temp = []
    for i in range(len(labels)):
        if (labels[i]) == n:
            cluster_temp.append(float(clustered[i]))
        else:
            cluster_temp.append(1)
    img_cluster[n, :, :, :] = np.array(cluster_temp).reshape(img_cs_cg.shape)

# Calculate norm of each cluster*img_cs to see which one is larger (and which one is signal versus noise)
norm_cluster0 = np.linalg.norm(img_cluster[0, :, :, :]*img_cs_cg)
norm_cluster1 = np.linalg.norm(img_cluster[1, :, :, :]*img_cs_cg)

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
img_noise = cv2.erode(img_noise, kernel, iterations=5)

# Dilate the signal image (optional)
kernel = np.ones((3, 3), np.uint8)
# img_signal = cv2.dilate(img_signal, kernel, iterations=2)

# Visualize
sp.plot.ImagePlot(2*img_signal - img_noise, x=0, y=1, z=2,
                  title="Signal and noise mask, estimated using k-means",
                  vmin=0, vmax=2, colormap="YlGnBu_r")

# Calculate SNR:
snr_nufft = np.mean((img_signal*img_nufft)[img_signal != 0]
                    / np.std((img_noise*img_nufft)[img_nufft != 0]))
print("SNR for inverse NUFFT: " + str(snr_nufft))
snr_cg = np.mean((img_signal*img_cg)[img_signal != 0]
                 / np.std((img_noise*img_cg)[img_cg != 0]))
print("SNR for iterative nufft: " + str(snr_cg))
snr_cs_cg = np.mean((img_signal*img_cs_cg)[img_signal != 0]
                    / np.std((img_noise*img_cs_cg)[img_cs_cg != 0]))
print("SNR for compressed sense: " + str(snr_cs_cg))


# %% Make path to save image results

# Check whether a specified save data path exists
results_exist = os.path.exists(location + "/results")

# Create a new directory because the results path does not exist
if not results_exist:
    os.makedirs(location + "/results")
    print("A new directory inside: " + location +
          " called 'results' has been created.")

# %% Save images as Nifti files

# Build an array using matrix multiplication
scaling_affine = np.array([[1, 0, 0, 0],
                           [0, 1, 0, 0],
                           [0, 0, 1, 0],
                           [0, 0, 0, 1]])

# Rotate gamma radians about axis i
cos_gamma = np.cos(0)
sin_gamma = np.sin(0)
rotation_affine_1 = np.array([[1, 0, 0, 0],
                              [0, cos_gamma, -sin_gamma,  0],
                              [0, sin_gamma, cos_gamma, 0],
                              [0, 0, 0, 1]])
cos_gamma = np.cos(np.pi)
sin_gamma = np.sin(np.pi)
rotation_affine_2 = np.array([[cos_gamma, 0, sin_gamma, 0],
                              [0, 1, 0, 0],
                              [-sin_gamma, 0, cos_gamma, 0],
                              [0, 0, 0, 1]])
cos_gamma = np.cos(0)
sin_gamma = np.sin(0)
rotation_affine_3 = np.array([[cos_gamma, -sin_gamma, 0, 0],
                              [sin_gamma, cos_gamma, 0, 0],
                              [0, 0, 1, 0],
                              [0, 0, 0, 1]])
rotation_affine = rotation_affine_1.dot(
    rotation_affine_2.dot(rotation_affine_3))

# Apply translation
translation_affine = np.array([[1, 0, 0, 0],
                               [0, 1, 0, 0],
                               [0, 0, 1, 0],
                               [0, 0, 0, 1]])

# Multiply matrices together
aff = translation_affine.dot(rotation_affine.dot(scaling_affine))

ni_img = nib.Nifti1Image(abs(img_nufft), affine=aff)
nib.save(ni_img, location + '/results/img_nufft_' +
         str(subsample_pc) + 'pc.nii.gz')

ni_img = nib.Nifti1Image(abs(img_cs_cg), affine=aff)
nib.save(ni_img, location + '/results/img_cs_' +
         str(subsample_pc) + 'pc.nii.gz')


# Save matlab files
savemat(location + "\\img_ventilation_nufft.mat",
        mdict={'img_ventilation_nufft': img_nufft})
savemat(location + "\\img_ventilation_cs_small_rho.mat",
        mdict={'img_ventilation_cs_small_rho': img_cs_cg})
savemat(location + "\\img_ventilation_mask.mat",
        mdict={'img_ventilation_mask': img_signal})
savemat(location + "\\img_ventilation_mask_noise.mat",
        mdict={'img_ventilation_mask_noise': img_noise})


#%% Ensure output path exists
os.makedirs(location, exist_ok=True)

# Convert to float32 (recommended for NIfTI)
img_nufft   = img_nufft.astype(np.float32)
img_cs_cg   = img_cs_cg.astype(np.float32)
img_signal  = img_signal.astype(np.float32)
img_noise   = img_noise.astype(np.float32)

# Save NIfTI files
nib.save(nib.Nifti1Image(img_nufft, aff),
         os.path.join(location, 'img_ventilation_nufft.nii.gz'))

nib.save(nib.Nifti1Image(img_cs_cg, aff),
         os.path.join(location, 'img_ventilation_cs_small_rho.nii.gz'))

nib.save(nib.Nifti1Image(img_signal, aff),
         os.path.join(location, 'img_ventilation_mask.nii.gz'))

nib.save(nib.Nifti1Image(img_noise, aff),
         os.path.join(location, 'img_ventilation_mask_noise.nii.gz'))

print("All NIfTI images saved successfully")

# %% Make a beautiful plot showing magnitude and SNR

# Select slice
plot_slice = int(recon_resolution*0.45)

# Generate images
img_nufft_slice = abs(img_nufft[:, :, plot_slice])
img_nufft_slice /= np.percentile(np.ravel(img_nufft_slice), 99.9)
img_nufft_snr_slice = abs(
    img_nufft[:, :, plot_slice]/np.std((img_noise*img_nufft)[img_noise != 0]))
img_cs_slice = abs(img_cs_cg[:, :, plot_slice])
img_cs_slice /= np.percentile(np.ravel(img_cs_slice), 99.9)
img_cs_snr_slice = abs(
    img_cs_cg[:, :, plot_slice]/np.std((img_noise*img_cs_cg)[img_noise != 0]))

# Plot
plt.figure(figsize=(3, 6), dpi=100)
fig, ((ax1, ax2)) = plt.subplots(1, 2, sharey=True, gridspec_kw={'wspace': 0, 'hspace': 0},
                                 squeeze=True)
img_plot1 = ax1.imshow(img_nufft_slice, cmap="gray", vmin=0, vmax=1)
img_plot2 = ax2.imshow(img_nufft_snr_slice, cmap="jet", vmin=0, vmax=35)
ax1.axis('off')
ax2.axis('off')

plt.figure(figsize=(3, 6), dpi=100)
fig, ((ax1, ax2)) = plt.subplots(1, 2, sharey=True, gridspec_kw={'wspace': 0, 'hspace': 0},
                                 squeeze=True)
img_plot1 = ax1.imshow(img_cs_slice, cmap="gray", vmin=0, vmax=1)
img_plot2 = ax2.imshow(img_cs_snr_slice, cmap="jet", vmin=0, vmax=35)
ax1.axis('off')
ax2.axis('off')

# %% Make an image summary

# Select slice
plot_slice = int(recon_resolution*0.45)

vmax_snr = np.percentile([util.calculate_snr(img_nufft, img_nufft[img_noise == 1]),
                          util.calculate_snr(img_cg, img_cg[img_noise == 1]),
                          util.calculate_snr(
                              img_cs_cg, img_cs_cg[img_noise == 1])], 99.9)

# Set up subplot
plt.figure()
fig, axes = plt.subplots(nrows=3, ncols=2, figsize=(6, 6))
if use_dcf:
    fig.suptitle('Image reconstructions: with DCF \n max_iter=%i, \u03C1=%.1e, $p_{tol}$=%.1e \n $A^TA$=%i, \u03BB$_1$=%.1e, \u03BB$_2$=%.1e' % (
        num_iters, rho, ptol, num_normal, lamda_1, lamda_2), y=1.03)
else:
    fig.suptitle('Image reconstructions: without DCF \n max_iter=%i, \u03C1=%.1e, $p_{tol}$=%.1e \n $A^TA$=%i, \u03BB$_1$=%.1e, \u03BB$_2$=%.1e' % (
        num_iters, rho, ptol, num_normal, lamda_1, lamda_2), y=1.03)

ax = axes[0, 0]
im = ax.imshow(img_nufft[:, :, plot_slice], cmap="gray")
ax.set_title(r'$F^HD(b)$')
fig.colorbar(im, ax=ax)
ax.axis('off')
ax = axes[0, 1]
snr = ax.imshow(util.calculate_snr(
    img_nufft[:, :, plot_slice], img_nufft[img_noise == 1]), cmap="jet", vmax=vmax_snr)
ax.set_title('snr, mean = %0.4f' % snr_nufft)
fig.colorbar(snr, ax=ax)
ax.axis('off')

ax = axes[1, 0]
im = ax.imshow(img_cg[:, :, plot_slice], cmap="gray")
ax.set_title(r'$||Ax-b||^2_2$')
fig.colorbar(im, ax=ax)
ax.axis('off')
ax = axes[1, 1]
snr = ax.imshow(util.calculate_snr(
    img_cg[:, :, plot_slice], img_cg[img_noise == 1]), cmap="jet", vmax=vmax_snr)
ax.set_title('snr, mean = %0.4f' % snr_cg)
fig.colorbar(snr, ax=ax)
ax.axis('off')

ax = axes[2, 0]
im = ax.imshow(img_cs_cg[:, :, plot_slice], cmap="gray")
ax.set_title(r'$||Ax-b||^2_2 + g(x)$ (cg)')
fig.colorbar(im, ax=ax)
ax.axis('off')
ax = axes[2, 1]
snr = ax.imshow(util.calculate_snr(
    img_cs_cg[:, :, plot_slice], img_cs_cg[img_noise == 1]), cmap="jet", vmax=vmax_snr)
ax.set_title('snr, mean = %0.4f' % snr_cs_cg)
fig.colorbar(snr, ax=ax)
ax.axis('off')

# %% create .exe

'''
    To create the exe file:
        	pip install pyinstaller
        	In terminal, run:
            pyinstaller --onefile --name="Philips_3D_XeVent_FLORET_recon.exe" Philips_3D_XeVent_FLORET_recon.py
This will create a single .exe file under \dist. 

'''