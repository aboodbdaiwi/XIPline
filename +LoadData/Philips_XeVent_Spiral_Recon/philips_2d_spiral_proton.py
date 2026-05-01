# -*- coding: utf-8 -*-
"""
Created on Sat Jun 20 16:47:16 2022

@author: PLUMD2
"""

# %% Import packages

import numpy as np
import matplotlib.pyplot as plt 
 
# Readphilips node - confidential
import ReadPhilips.readphilips as rp
#from readphilips.file_io import io

# Sigpy functionality
import sigpy as sp
import sigpy.mri as mr
#import sigpy.plot as pl

# Nifti functionality
import nibabel as nib
import os

# %% Reconstruction settings
debug_mod = False
devnum = -1
device = sp.Device(devnum)
xp = device.xp
def mvd(x): return sp.to_device(x, device)
def mvc(x): return sp.to_device(x, sp.cpu_device)

# Density compensation
use_dcf = True

# Iterative recon
iter_recon = False

# Reconstruct onto larger volume and crop image
crop_image = True

# Estimate a sensitivity map using JSense (not really useful unless doing CS)
estimate_mps = False

# %% Read data

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
        if file.endswith(".data"):
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

# Load data
inputfile = rp.PhilipsData(datafile)
inputfile.compute()
data_load = inputfile.data
data_load = np.squeeze(data_load)
if np.size(np.shape(data_load)) == 5:
    data_load = data_load[0, :, :, :, :]
    print('Multiple dynamics assumed on input data.')
if np.size(np.shape(data_load)) == 3:
    data_load = np.reshape(data_load, (1, np.shape(data_load)[
                           0], np.shape(data_load)[1], np.shape(data_load)[2]))
    print('Extra dimension added for channel.')
plot_slice = int(np.shape(data_load)[0]/2)
order = inputfile.header.get('list').get('ky').astype(np.int)  # Spiral order
fid = abs(data_load[0, plot_slice, 0, :])
if debug_mod:
    plt.figure()
    plt.plot(fid)
    plt.ylabel('k-space intensity')
    plt.xlabel('Sample number')
    plt.title('Free induction decay')
    plt.show()

# Load traj
inputfile = rp.PhilipsData(trajfile)
inputfile.compute()
coords = inputfile.spparams.get('COORDS')
coords_expanded = inputfile.spparams.get('COORDS_EXPANDED')
if debug_mod:
    plt.figure()
    plt.scatter(coords_expanded[0, :, :],
                coords_expanded[1, :, :],
                linewidths=0.1,
                marker=".",
                edgecolor="blue",
                s=2)
    plt.ylabel('$k_y$')
    plt.xlabel('$k_x$')
    plt.title('Spiral coordinates')
    plt.axis('square')
    plt.show()

# Understand data
N_channels = np.shape(data_load)[0]
N_slices = np.shape(data_load)[1]
N_spirals = np.shape(data_load)[2]
N_samp = np.shape(data_load)[3]
order = order[::2]  # Deal with second dynamic (only on some cases)
# First element is for noise if you look at list file
order = order[1:N_spirals+1]

# Extract useful information from header
scan_date = inputfile.header.get('sin').get('start_scan_date_time')[0]
voxel_size = float(inputfile.header.get('sin').get('voxel_sizes')[0][0])
slice_thickness = float(inputfile.header.get(
    'sin').get('slice_thickness')[0][0])
slice_spacing = inputfile.header.get(
    'sin').get('spacing_between_slices_arr')[0]
recon_resolution = int(inputfile.header.get(
    'sin').get('recon_resolutions')[0][0])
field_of_view = voxel_size * recon_resolution
scan_resolution = int(inputfile.header.get(
    'sin').get('scan_resolutions')[0][0])
location_center_coordinates = inputfile.header.get(
    'sin').get('location_center_coordinates')[0]
loc_ap_rl_fh_offcentres = inputfile.header.get(
    'sin').get('loc_ap_rl_fh_offcentres')[0]

# %% Reorder data and trajectories to match Philips ordering

coords_reordered = np.zeros(np.shape(coords_expanded))
data_reordered = np.zeros(np.shape(data_load), dtype='complex_')
for ii in range(N_spirals):
    coords_reordered[:, :, ii] = coords_expanded[:, :, order[ii]]
    data_reordered[:, :, ii, :] = data_load[:, :, order[ii], :]

# Confirm reordering
if debug_mod:
    fig, (ax1, ax2) = plt.subplots(2, sharex=True)
    ax1.plot(abs(data_load[0, plot_slice, :, 0]), color='r')
    ax1.set_ylabel('k0 intensity')
    ax1.set_title("k0 intensities before reordering")
    ax2.plot(abs(data_reordered[0, plot_slice, :, 0]), color='b')
    ax2.set_ylabel('k0 intensity')
    ax2.set_xlabel('Spiral number')
    ax2.set_title("k0 intensities after reordering")
    fig.tight_layout(pad=1.0, w_pad=0.5, h_pad=1.0)

# %% Make manipulations to reconstructed volume
if crop_image:
    alpha = 0.5  # Additional rescale factor to allow cropping
    ImageSize = scan_resolution / alpha
else:
    alpha = 1
    ImageSize = scan_resolution

resize_factor = recon_resolution/(scan_resolution)
ImageSize = int(resize_factor * ImageSize)
ImageShape = (int(ImageSize), int(ImageSize))
traj = coords_reordered[0:2, :, :] * 1 / resize_factor
traj = ImageSize * traj
traj = np.transpose(traj)

# %% Perform phase demodulation to recenter image according to x-y offset in .sin file

# Define offset from the sin file
offset_x = - float(location_center_coordinates[1])
offset_y = - float(location_center_coordinates[2])

# Convert to matrix size coordinates
offset_x = alpha * offset_x/ImageSize
offset_y = alpha * offset_y/ImageSize

# Phase modulation scale factor
translate_factor = np.exp(- 1j * 2 * np.pi *
                          (offset_x * traj[:, :, 0] + offset_y * traj[:, :, 1]))

# Plot the demodulation
if debug_mod:
    plt.figure()
    plt.plot(abs(translate_factor[0, :]), color='r')
    plt.plot(np.real(translate_factor[0, :]), color='b')
    plt.plot(np.imag(translate_factor[0, :]), color='y')
    plt.legend(['Magnitude', 'Real', 'Imaginary'])
    plt.ylabel('Scale factor')
    plt.xlabel('Sample number')
    plt.title('Phase demodulation for first projection')
    plt.show()

# Perform translation by multiplying k-space data
data_reordered = translate_factor * data_reordered

# %% Initialize image data

# Initialize proton image
# Dimensions to deal with coronal where images are in scanner (x-z) dimension
img_proton = np.zeros((ImageSize, ImageSize, N_slices))

# Initialize sensitivity map
mps_estimate = np.ones(
    (N_channels, ImageSize, ImageSize, N_slices), dtype="complex")

# %% Estimate density compensation functions

# DCF settings
beta = 8  # 8
width = 3  # 4

# Full image
dcf_full = mr.pipe_menon_dcf(
    traj, img_shape=ImageShape, beta=beta, width=width, device=device)

# %% Loop through image slices

for imslice in range(N_slices):
    data = mvd(data_reordered[:, imslice, :, :])

    # Estimate sensitivity map
    if estimate_mps:
        mps_estimate[:, :, :, imslice] = mr.app.JsenseRecon(
            data, coord=traj, weights=dcf_full,
            mps_ker_width=8, ksp_calib_width=12, lamda=0,
            img_shape=ImageShape, device=sp.to_device(0)).run()
    mps = mps_estimate[:, :, :, imslice]

    if iter_recon:

        # Iteration settings
        max_iter = 30
        lamda = 0

        if use_dcf:
            img_proton[:, :, imslice] = mvc(abs(mr.app.SenseRecon(
                data, mps, weights=dcf_full, coord=traj, lamda=lamda, max_iter=max_iter, device=device).run()))

        else:
            img_proton[:, :, imslice] = mvc(abs(mr.app.SenseRecon(
                data, mps, coord=traj, lamda=lamda, max_iter=max_iter, device=device).run()))

    else:

        S = sp.linop.Multiply(ImageShape, mps)
        F = sp.linop.NUFFT(mps.shape,
                           coord=traj,
                           oversamp=1.25,
                           width=3,
                           toeplitz=True)
        if use_dcf:
            D = sp.linop.Multiply(F.oshape, dcf_full**0.5)  # Optional
            A = D * F * S
            b = data * dcf_full**0.5
            # b = b/xp.linalg.norm(mvd(data_load) * dcf_full**0.5) # Optional: comment out for quantitative analysis

        else:
            A = F * S
            b = data
            # b = b/xp.linalg.norm(mvd(data_load)) # Optional: comment out for quantitative analysis

        # Normalize data and A operator for full image only (makes image come out at true intensity before dcf)
        LL = sp.app.MaxEig(A.N, dtype=xp.complex64, device=device,
                           max_iter=10, show_pbar=False, leave_pbar=False).run() * 1.01
        A = np.sqrt(1/LL) * A

        # Compute a single x = A.H b operation (i.e. inverse NUFFT)
        # Automatically sums along channel dimension
        img_grid = mvc(abs(A.H * b))
        img_proton[:, :, imslice] = img_grid

# %% Crop images from oversampled reconstruction
if crop_image:

    # Define center coordinates
    center_x = (field_of_view - 1)//(2*voxel_size*alpha)
    center_y = (field_of_view - 1)//(2*voxel_size*alpha)

    # Define crop boundaries
    crop_x = (np.arange(center_x - ((recon_resolution)//2),
              center_x + ((recon_resolution)//2), 1, dtype=int))
    crop_y = (np.arange(center_y - ((recon_resolution)//2),
              center_y + ((recon_resolution)//2), 1, dtype=int))

    # Make crop function
    def crp(x): return x[crop_x, :, :][:, crop_y, :]

    img_proton = crp(img_proton)

# %% Apply pixel shifts and rotations

# Translation and rotation
dx = 0
dy = 0
rot = 2

# Define functions


def ps(x): return np.roll(x, (dx, dy), axis=(0, 1))  # pixelshift
def ir(x): return np.rot90(abs(x), k=rot, axes=(0, 1))  # rotate image

# %% Visualize results

def montage3d(volume, cols=5, padding=2):
    """Create a montage from a 3D volume (shape: H x W x Slices)."""
    H, W, N = volume.shape
    rows = int(np.ceil(N / cols))

    # Normalize volume to [0, 1]
    vol = volume.astype(np.float32)
    vol -= vol.min()
    if vol.max() > 0:
        vol /= vol.max()

    # Add padding between slices
    pad_value = 0  # black
    montage = np.full((
        rows * (H + padding) - padding,
        cols * (W + padding) - padding
    ), pad_value, dtype=np.float32)

    for i in range(N):
        r, c = divmod(i, cols)
        y = r * (H + padding)
        x = c * (W + padding)
        montage[y:y+H, x:x+W] = vol[:, :, i]

    return montage

# Create montage and plot
if debug_mod:
    montage_img = montage3d(img_proton)
    plt.figure(figsize=(12, 12))
    plt.imshow(montage_img, cmap='gray')
    plt.axis('off')
    plt.title('Montage of 3D Volume')
    plt.show()
# %% Save images as Nifti files

# Build an array using matrix multiplication
scaling_affine = np.array([[voxel_size, 0, 0, 0],
                           [0, 0, -slice_thickness, 0],
                           [0, voxel_size, 0, 0],
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
translation_affine = np.array([[1, 0, 0, (field_of_view)/2 - float(location_center_coordinates[1])],
                               [0, 1, 0, float(
                                   location_center_coordinates[0])],
                               [0, 0, 1, (field_of_view)/2 -
                                float(location_center_coordinates[2])],
                               [0, 0, 0, 1]])

# Multiply matrices together
aff = translation_affine.dot(rotation_affine.dot(scaling_affine))

# Save images
ni_img = nib.Nifti1Image(img_proton, affine=aff)
nib.save(ni_img, datalocation + '/img_proton.nii.gz')
print("Image has been reconstructed successfully")

# %% create .exe

'''
    To create the exe file:
        	pip install pyinstaller
        	In terminal, run:
            pyinstaller --onefile --name="Philips_2DHSpiralRecon.exe" philips_2d_spiral_proton.py
This will create a single .exe file under \dist. 

'''
