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
import os
import nibabel as nib
from functions import tv, convexalg, util
import sigpy_local.mri as mr
import sigpy_local.plot as pl
import sigpy_local as sp
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

try:
    devnum = 0  # GPU
    device = sp.Device(devnum)
    xp = device.xp
    print(f"Using GPU (device {devnum})")
except Exception as e:
    print(f"GPU not available or failed: {e}")
    print("Falling back to CPU")

    devnum = -1
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

# Data folder location
def read_config_and_find_files(config_path=r"C:\XIPline\offline_recon\config_H.txt"):
    # Read the config file
    with open(config_path, 'r') as f:
        lines = f.readlines()

    # Extract datalocation and freq_offset
    datalocation = None
    for line in lines:
        line = line.strip()
        if line.lower().startswith("datalocation"):
            datalocation = line.split('=', 1)[1].strip()

    if datalocation is None:
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

    return datafile_r, trajfile_r, datalocation_r

config_H_path=r"C:\XIPline\offline_recon\config_H.txt"
datafile, trajfile, datalocation  = read_config_and_find_files(config_H_path)
print("Data file:", datafile)
print("Trajectory file:", trajfile)
print("Data location:", datalocation)


#%%
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
data_load = np.reshape(data_load, (np.shape(data_load)[0], np.shape(data_load)[1]*np.shape(data_load)[2], np.shape(data_load)[3]))
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
N_samp = np.shape(data_load)[2]  # samples per projection
N_proj = np.shape(data_load)[1]  # number of projections

# Plot trajectories:
plt.figure()
ax = plt.axes(projection='3d')
N_visual = 500  # Number of projections you want to show, inefficient for large n
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

# Extract useful information from header
scan_date = traj.header.get('sin').get('start_scan_date_time')[0]
voxel_size = float(traj.header.get('sin').get('voxel_sizes')[0][0])
slice_thickness = float(traj.header.get(
    'sin').get('slice_thickness')[0][0])
slice_spacing = traj.header.get(
    'sin').get('spacing_between_slices_arr')[0]
recon_resolution = int(traj.header.get(
    'sin').get('recon_resolutions')[0][0])
field_of_view = voxel_size * recon_resolution
scan_resolution = int(traj.header.get(
    'sin').get('scan_resolutions')[0][0])
location_center_coordinates = traj.header.get(
    'sin').get('location_center_coordinates')[0]
loc_ap_rl_fh_offcentres = traj.header.get(
    'sin').get('loc_ap_rl_fh_offcentres')[0]

# Define some parameters manually
recon_resolution = int(np.ceil(1.33 * scan_resolution))
scan_resolution  = int(np.ceil(1.33 * scan_resolution))

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
ax1.plot(abs(data_load[0,:, 0]), color='r')
ax1.set_ylabel('k0 intensity')
ax1.set_title("magnetization decay")
ax1.set_xlabel("projection number")
ax2.plot(abs(data_load[0,0, :]), color='b')
ax2.set_ylabel('k-space intensity')
ax2.set_xlabel('sample number')
ax2.set_title("fid")
fig.tight_layout(pad=1.0, w_pad=0.5, h_pad=1.0)

# Normalize data to max k0
# k0 for each channel and projection
k0 = np.abs(data_load[:, :, 0])

# Normalize
data_normalized = data_load / k0[:, :, None]

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
data_load_subset = data_load[0,:subset, :]
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
#data_load_channel = np.reshape(data_load_subset, (1, np.shape(
#    data_load_subset)[0], np.shape(data_load_subset)[1]))

data_load_channel = data_load[:,:subset, :]
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

# Generate proton weighting for inverse NUFFT
# Proton data does not need hyperpolarized RF-decay correction

k = np.ones((N_proj, N_samp), dtype=np.float32)

# Visualize
fig, ax = plt.subplots(1, 1)
ax.plot(np.abs(k[0, :]))
ax.set_ylabel('weight')
ax.set_xlabel('sample number')
ax.set_title('Proton NUFFT weighting')
fig.tight_layout()

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

    def mvd(x):
        return sp.to_device(x, device)

    def mvc(x):
        return sp.to_device(x, sp.cpu_device)

    # Move data to device
    data_all = mvd(data_load_channel)

    Nchannels = data_all.shape[0]

    img_nufft_ch = []
    #img_cg_ch = []
    #img_cs_cg_ch = []

    for ch in range(Nchannels):

        print(f"\nReconstructing channel {ch+1}/{Nchannels}")

        data = data_all[ch, :, :]

        # -----------------------------
        # Normalize data
        # -----------------------------
        if use_dcf:
            data_norm_dcf = xp.linalg.norm(data * dcf**0.5)
            b = data * dcf**0.5 / data_norm_dcf
        else:
            data_norm = xp.linalg.norm(data)
            b = data / data_norm

        # -----------------------------
        # Inverse NUFFT
        # -----------------------------
        A_dcf = D * F * S

        if use_dcf:
            b_dcf = data * dcf**0.5 / xp.linalg.norm(data * dcf**0.5)
        else:
            b_dcf = data / xp.linalg.norm(data)

        img_nufft_tmp = mvc(A_dcf.H * b_dcf)
        img_nufft_ch.append(img_nufft_tmp)

        # -----------------------------
        # PSF, only need once
        # -----------------------------
        if ch == 0:
            b_dcf_psf = xp.ones_like(b) * dcf**0.5
            img_nufft_psf = mvc(abs(A_dcf.H * b_dcf_psf))
        ''' 
        # -----------------------------
        # CG reconstruction
        # -----------------------------
        img_cg_tmp = mvc(convexalg.cg(
            num_iters=200,
            ptol=1e-2,
            A=A,
            b=b,
            lamda=lamda_l2,
            verbose=True,
            draw_output=True,
            save=None
        ))

        img_cg_ch.append(img_cg_tmp)

        # -----------------------------
        # CS reconstruction
        # -----------------------------
        img_cs_cg_tmp = mvc(convexalg.admm(
            num_iters=num_iters,
            ptol=ptol,
            A=A,
            b=b,
            num_normal=num_normal,
            lst_proxg=lst_proxg,
            rho=rho,
            lst_g=lst_g,
            method="cg",
            verbose=True,
            draw_output=True
        ))

        img_cs_cg_ch.append(img_cs_cg_tmp)
    '''
    # Convert list to arrays
    img_nufft_ch = np.asarray(img_nufft_ch)
    #img_cg_ch = np.asarray(img_cg_ch)
    #img_cs_cg_ch = np.asarray(img_cs_cg_ch)

    # Shape should be:
    # (Nchannels, Nx, Ny, Nz)

    # -----------------------------
    # Channel combination
    # Root-sum-of-squares
    # -----------------------------
    img_nufft = np.sqrt(np.sum(np.abs(img_nufft_ch)**2, axis=0))
    #img_cg = np.sqrt(np.sum(np.abs(img_cg_ch)**2, axis=0))
    #img_cs_cg = np.sqrt(np.sum(np.abs(img_cs_cg_ch)**2, axis=0))

    # -----------------------------
    # Display combined images
    # -----------------------------
    pl.ImagePlot(
        np.rot90(img_nufft[:, :, int(recon_resolution/3):int(recon_resolution/1.5):2], 2),
        x=0, y=1, z=2,
        title="Combined reconstruction using inverse NUFFT"
    )
    
    ''' 
    pl.ImagePlot(
        np.rot90(img_cg[:, :, int(recon_resolution/3):int(recon_resolution/1.5):2], 2),
        x=0, y=1, z=2,
        title="Combined reconstruction using CG"
    )

    pl.ImagePlot(
        np.rot90(img_cs_cg[:, :, int(recon_resolution/3):int(recon_resolution/1.5):2], 2),
        x=0, y=1, z=2,
        title="Combined reconstruction using CS: ADMM w/ CG"
    )    
    '''    
    
    
# %% Manipulate images for plotting

# Translation and rotation
dx = 0
dy = 0
rot = 1

img_nufft = np.flipud(util.ps(util.ir(img_nufft, rot), dx, dy))
#img_cg = np.flipud(util.ps(util.ir(img_cg, rot), dx, dy))
#img_cs_cg = np.flipud(util.ps(util.ir(img_cs_cg, rot), dx, dy))

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


# Save matlab files
savemat(location + "\\img_proton_nufft.mat",
        mdict={'img_proton_nufft': img_nufft})



#%% Ensure output path exists
os.makedirs(location, exist_ok=True)

# Convert to float32 (recommended for NIfTI)
img_nufft   = img_nufft.astype(np.float32)

# Save NIfTI files
nib.save(nib.Nifti1Image(img_nufft, aff),
         os.path.join(location, 'img_proton_nufft.nii.gz'))
''' 
nib.save(nib.Nifti1Image(img_cs_cg, aff),
         os.path.join(location, 'img_ventilation_cs_small_rho.nii.gz'))

nib.save(nib.Nifti1Image(img_signal, aff),
         os.path.join(location, 'img_ventilation_mask.nii.gz'))

nib.save(nib.Nifti1Image(img_noise, aff),
         os.path.join(location, 'img_ventilation_mask_noise.nii.gz'))
'''
print("All NIfTI images saved successfully")

# %% create .exe

'''
    To create the exe file:
        	pip install pyinstaller
        	In terminal, run:
            -------   use xe_recon env
            pyinstaller --onefile --name="Philips_3D_XeVent_FLORET_recon.exe" Philips_3D_XeVent_FLORET_recon.py
This will create a single .exe file under \dist. 
rmdir /s /q build
rmdir /s /q dist
del Philips_3D_proton_FLORET_recon.exe.spec

pyinstaller --onefile ^
  --name="Philips_3D_proton_FLORET_recon" ^
  --hidden-import=fastrlock ^
  --hidden-import=cupy ^
  --hidden-import=cupy_backends.cuda._softlink ^
  --hidden-import=cupy_backends.cuda.libs ^
  --collect-all cupy ^
  --collect-all cupy_backends ^
  --collect-all fastrlock ^
  Philips_3D_proton_FLORET_recon.py
  
'''