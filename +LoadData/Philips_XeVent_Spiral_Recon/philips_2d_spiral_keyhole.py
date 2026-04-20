# -*- coding: utf-8 -*-
"""
Created on Sat Jun 20 16:47:16 2022

@author: PLUMD2
"""

# %% Import packages

#import cupy as cp
from skimage.transform import resize
from keras.models import load_model
from scipy.io import savemat
#import ants
import cv2
from scipy import ndimage
from sklearn.cluster import KMeans
import nibabel as nib
import sigpy.mri as mr
import sigpy.plot
import sigpy as sp
from readphilips.file_io import io
import readphilips as rp
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
plt.style.use('dark_background')

# %% Reconstruction settings

devnum = -1
device = sp.Device(devnum)
xp = device.xp
def mvd(x): return sp.to_device(x, device)
def mvc(x): return sp.to_device(x, sp.cpu_device)


# Density compensation
use_dcf = True

# Iterative recon
iter_recon = False

# Key radius
# This is a percentage of the maximum distance in k-space (recommended 0.025 - 0.075)
r_key = 0.04

# Offset correction
freq_offset = 150  # Hz

# Reconstruct onto larger volume and crop image
crop_image = True

# Select dynamic
dynamic = 1

# Run N4 bias correction
n4_correction = False

# %% Read data

# Data folder location
#location = io(
#    initial_dir="/mnt/c/Users/PLUMD2/Desktop/Images for Keyhole Reconstruction/").selectDirectory()
location = r"\\rds6.cchmc.org\PulMed-43\CPIR_Share\Carter\ForAbood_10Jun2025_IRC1031-140c\offline_recon_asb\gas"
# Load data
#filename = io(initial_dir=location, f_type="raw-list files",
#              ext_type="*.data").selectFile()
filename = r"\\rds6.cchmc.org\PulMed-43\CPIR_Share\Carter\ForAbood_10Jun2025_IRC1031-140c\offline_recon_asb\gas\raw_298.data"
inputfile = rp.PhilipsData(filename)
inputfile.compute()
data_load = np.array(inputfile.data)
data_load = np.squeeze(data_load)
N_dynamics = 1  # Default number of dynamics
if np.size(np.shape(data_load)) == 4:
    print('Multiple dynamics assumed on input data.')
    N_dynamics = data_load.shape[0]
    data_load = data_load[(dynamic - 1), :, :, :]  # Multiple dynamic data
elif np.size(np.shape(data_load)) == 5:
    print('Multiple dynamics and interleaving assumed on input data.')
    # Interleaved data for Abood
    N_dynamics = data_load.shape[1]
    data_load = data_load[0, (dynamic - 1), :10, :, :]
else:
    dynamic = 1
    print('Only a single dynamic was found.')
plot_slice = int(np.shape(data_load)[0]/2)
order = inputfile.header.get('list').get(
    'ky').astype(np.int)  # Spiral order
fid = abs(data_load[plot_slice, 0, :])
plt.figure()
plt.plot(fid)
plt.ylabel('k-space intensity')
plt.xlabel('Sample number')
plt.title('Free induction decay')
plt.show()

# Load traj
#filename = io(initial_dir=location, f_type="raw-lab-sin files",
#              ext_type="*.sin").selectFile()
filename = r"\\rds6.cchmc.org\PulMed-43\CPIR_Share\Carter\ForAbood_10Jun2025_IRC1031-140c\offline_recon_asb\gas\20250610_094435_CPIR_Vent_2D_VarDens_SOS_WIP.sin"
inputfile = rp.PhilipsData(filename)
inputfile.compute()
coords = inputfile.spparams.get('COORDS')
coords_expanded = inputfile.spparams.get('COORDS_EXPANDED')
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
N_slices = np.shape(data_load)[0]
N_spirals = np.shape(data_load)[1]
N_samp = np.shape(data_load)[2]
print("Received data with the following dimensions:")
print("N_slices = " + str(N_slices))
print("N_spirals = " + str(N_spirals))
print("N_samp = " + str(N_samp))
# order = order[::4] # Get every 4th element if diffusion
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
dwell_time = 1e-6*float(inputfile.header.get(
    'sin').get('sample_time_interval')[0][0])  # s
echo_time = 1e-3*float(inputfile.header.get(
    'sin').get('echo_times')[0][0])  # s
readout_duration = dwell_time * N_samp  # ms
flip_angle_applied = float(
    inputfile.header.get('sin').get('flip_angles')[0][0])

# %% Correct for frequency offset

dwell_time = readout_duration/N_samp
for ii in range(N_samp):
    ii = int(ii)
    scale_factor = np.exp(-2*np.pi*1j*(ii*dwell_time + echo_time)*freq_offset)
    data_load[:, :, ii] = scale_factor*data_load[:, :, ii]

# %% Reorder data and trajectories to match Philips ordering

coords_reordered = np.zeros_like(coords_expanded)
data_reordered = np.zeros_like(data_load, dtype='complex_')
for ii in range(N_spirals):
    coords_reordered[:, :, ii] = coords_expanded[:, :, order[ii]]
    data_reordered[:, ii, :] = data_load[:, order[ii], :]

# Confirm reordering
fig, (ax1, ax2) = plt.subplots(2, sharex=True)
ax1.plot(abs(data_load[plot_slice, :, 0]), color='r')
ax1.set_ylabel('k0 intensity')
ax1.set_title("k0 intensities before reordering")
ax2.plot(abs(data_reordered[plot_slice, :, 0]), color='b')
ax2.set_ylabel('k0 intensity')
ax2.set_xlabel('Spiral number')
ax2.set_title("k0 intensities after reordering")
fig.tight_layout(pad=1.0, w_pad=0.5, h_pad=1.0)

# %% Make manipulations to reconstructed volume
if crop_image:
    CROP_FAC = 0.5  # Additional rescale factor to allow cropping
    ImageSize = scan_resolution / CROP_FAC
else:
    CROP_FAC = 1
    ImageSize = scan_resolution

resize_factor = recon_resolution/(scan_resolution)
ImageSize = int(resize_factor * ImageSize)
ImageShape = (int(ImageSize), int(ImageSize))
traj = coords_reordered[0:2, :, :] * 1 / resize_factor
traj = ImageSize * traj  # For image size increase
traj = np.transpose(traj)

# %% Perform phase demodulation to recenter image according to x-y offset in .sin file

# Define offset from the sin file
offset_x = - float(location_center_coordinates[1])
offset_y = - float(location_center_coordinates[2])

# Convert to matrix size coordinates
offset_x = CROP_FAC * offset_x/ImageSize
offset_y = CROP_FAC * offset_y/ImageSize

# Phase modulation scale factor
translate_factor = np.exp(- 1j * 2 * np.pi *
                          (offset_x * traj[:, :, 0] + offset_y * traj[:, :, 1]))

# Plot the demodulation
plt.figure()
plt.plot(abs(translate_factor[0, :]), color='r')
plt.plot(np.real(translate_factor[0, :]), color='b')
plt.plot(np.imag(translate_factor[0, :]), color='y')
plt.legend(['Magnitude', 'Real', 'Imaginary'])
plt.ylabel('Scale factor')
plt.xlabel('Sample number')
plt.title('Phase demodulation for first projection')
plt.show()

# # Perform translation by multiplying k-space data
data_reordered = translate_factor * data_reordered

# %% Estimate sensitivity map

sens_map = np.ones((1, ImageSize, ImageSize))
# As we only have 1 channel, a sensitivity map of all 1's is assumed

# Move to GPU
mps = mvd(sens_map)

# %% Declare cut-offs for key/keyhole segregation

# Find radial k-space coordinates
traj_rad = np.sqrt(traj[:, :, 0]**2 + traj[:, :, 1]**2)

# Declare a cut off for key vs keyhole
key_rad = r_key*np.max(traj_rad[:, :])

# Define number of spiral interleaves to be allocated to each key
N_spirals_key1 = int(N_spirals/2)
N_spirals_key2 = N_spirals - N_spirals_key1
range_spirals_key1 = np.arange(0, N_spirals_key1)
range_spirals_key2 = np.arange(N_spirals_key1, N_spirals)

# Find values that fit within key radius
traj_key_bool = traj_rad < key_rad
N_key_samp = np.sum(traj_key_bool[0])

# %% Initialize image data

# Initialize
# Dimensions to deal with coronal where images are in scanner (x-z) dimension
img_full = np.zeros((ImageSize, ImageSize, N_slices))
img_key1 = np.zeros((ImageSize, ImageSize, N_slices))
img_key2 = np.zeros((ImageSize, ImageSize, N_slices))

# %% Preallocate key image trajectories

# Key 1
traj_key1 = np.zeros_like(traj)
traj_key1[range_spirals_key1, :, :] = traj[range_spirals_key1, :, :]
traj_key1[range_spirals_key2, N_key_samp:N_samp,
          :] = traj[range_spirals_key2, N_key_samp:N_samp, :]

# Key 2
traj_key2 = np.zeros_like(traj)
# Must ensure first value is not a nan
traj_key2[range_spirals_key2, :, :] = traj[range_spirals_key2, :, :]
traj_key2[range_spirals_key1, N_key_samp:N_samp,
          :] = traj[range_spirals_key1, N_key_samp:N_samp, :]

# %% Optional - remove keyhole

# TODO: Experiment with this -- does it make a difference? I think it does
# data_reordered = data_reordered[:, :, :N_key_samp]
# traj = traj[:, :N_key_samp, :]
# traj_key1 = traj_key1[:, :N_key_samp, :]
# traj_key2 = traj_key2[:, :N_key_samp, :]
# N_samp = traj.shape[1]


# %% Reshape trajectories and remove the zero (empty key) elements

# Full
traj = np.reshape(traj, (N_samp*N_spirals, 2))

# Key 1
traj_key1 = np.reshape(traj_key1, (N_samp*N_spirals, 2))
zero_row_mask1 = (traj_key1[:, 0] == 0) & (traj_key1[:, 1] == 0)
traj_key1 = traj_key1[~zero_row_mask1]

# Key 2
traj_key2 = np.reshape(traj_key2, (N_samp*N_spirals, 2))
zero_row_mask2 = (traj_key2[:, 0] == 0) & (traj_key2[:, 1] == 0)
traj_key2 = traj_key2[~zero_row_mask2]

# %% Estimate density compensation functions
if use_dcf:
    # DCF settings
    beta = 8  # 8
    width = 3  # 4

    # Full image
    dcf_full = mr.pipe_menon_dcf(
        traj, img_shape=ImageShape, beta=beta, width=width, device=device)
    dcf_full = mvd(np.reshape(dcf_full, (1, np.shape(dcf_full)[0])))

    # Modify DCF to deal with undersampled keys more efficiently
    beta = 8  # 6 # 8
    width = 3  # 7 # 4

    # Key image 1
    dcf_key1 = mr.pipe_menon_dcf(
        traj_key1, img_shape=ImageShape, beta=beta, width=width, device=device)
    dcf_key1 = mvd(np.reshape(dcf_key1, (1, np.shape(dcf_key1)[0])))

    # Key image 2
    dcf_key2 = mr.pipe_menon_dcf(
        traj_key2, img_shape=ImageShape, beta=beta, width=width, device=device)
    dcf_key2 = mvd(np.reshape(dcf_key2, (1, np.shape(dcf_key2)[0])))

# %% Loop through image slices

for imslice in range(N_slices):
    data = data_reordered[imslice, :, :]

    # Normalize k-space data to mean k0 of each key

    # Initialize variables
    data[data[:, 0] == 0] = 1e-23
    k0 = abs(data[:, 0])
    data_norm = np.zeros_like(data, dtype='complex_')
    data_norm_key1 = np.zeros_like(data, dtype='complex_')
    data_norm_key2 = np.zeros_like(data, dtype='complex_')

    # Calculate mean k0 values
    mean_k0 = abs(np.mean(k0[0:N_spirals]))
    mean_k0_key1 = abs(np.mean(k0[range_spirals_key1]))
    mean_k0_key2 = abs(np.mean(k0[range_spirals_key2]))

    for ii in range(N_spirals):
        data_norm[ii, :] = data[ii, :] * mean_k0 / k0[ii]
        data_norm_key1[ii, :] = data[ii, :] * mean_k0_key1 / k0[ii]
        data_norm_key2[ii, :] = data[ii, :] * mean_k0_key2 / k0[ii]

    # Segregate k-space data into keys and keyhole

    # Key 1
    data_key1 = np.zeros_like(data_norm_key1, dtype='complex_')
    data_key1[range_spirals_key1,
              :] = data_norm_key1[range_spirals_key1, :]
    data_key1[range_spirals_key2,
              N_key_samp:N_samp] = data_norm_key1[range_spirals_key2, N_key_samp:N_samp]

    # Key 2
    data_key2 = np.zeros_like(data_norm_key2, dtype='complex_')
    data_key2[range_spirals_key2,
              :] = data_norm_key2[range_spirals_key2, :]
    data_key2[range_spirals_key1,
              N_key_samp:N_samp] = data_norm_key2[range_spirals_key1, N_key_samp:N_samp]

    # Reshape into vectors and remove the zero elements

    # Full
    data_norm = np.reshape(data_norm, (N_samp*N_spirals))

    # Key 1
    data_key1 = np.reshape(data_key1, (N_samp*N_spirals))
    data_key1 = data_key1[~zero_row_mask1]

    # Key 2
    data_key2 = np.reshape(data_key2, (N_samp*N_spirals))
    data_key2 = data_key2[~zero_row_mask2]

    # Compute linear operators and reconstruct images

    if iter_recon:
        # Iteration settings
        max_iter = 30
        lamda = 0

        if use_dcf:
            # Full image
            data_norm = mvd(np.reshape(data_norm, (1, np.shape(data_norm)[0])))
            img_full[:, :, imslice] = mvc(abs(mr.app.SenseRecon(
                data_norm, mps, weights=dcf_full, coord=traj, lamda=lamda, max_iter=max_iter, device=device).run()))

            # Key image 1
            data_key1 = mvd(np.reshape(data_key1, (1, np.shape(data_key1)[0])))
            img_key1[:, :, imslice] = mvc(abs(mr.app.SenseRecon(
                data_key1, mps, coord=traj_key1, weights=dcf_key1, lamda=lamda, max_iter=max_iter, device=device).run()))

            # Key image 2
            data_key2 = mvd(np.reshape(data_key2, (1, np.shape(data_key2)[0])))
            img_key2[:, :, imslice] = mvc(abs(mr.app.SenseRecon(
                data_key2, mps, coord=traj_key2, weights=dcf_key2, lamda=lamda, max_iter=max_iter, device=device).run()))

        else:
            # Full image
            data_norm = mvd(np.reshape(data_norm, (1, np.shape(data_norm)[0])))
            img_full[:, :, imslice] = mvc(abs(mr.app.SenseRecon(
                data_norm, mps, coord=traj, lamda=lamda, max_iter=max_iter, device=device).run()))

            # Key image 1
            data_key1 = mvd(np.reshape(data_key1, (1, np.shape(data_key1)[0])))
            img_key1[:, :, imslice] = mvc(abs(mr.app.SenseRecon(
                data_key1, mps, coord=traj_key1, lamda=lamda, max_iter=max_iter, device=device).run()))

            # Key image 2
            data_key2 = mvd(np.reshape(data_key2, (1, np.shape(data_key2)[0])))
            img_key2[:, :, imslice] = mvc(abs(mr.app.SenseRecon(
                data_key2, mps, coord=traj_key2, lamda=lamda, max_iter=max_iter, device=device).run()))

    else:

        # Full image
        S = sp.linop.Multiply(ImageShape, mps)
        F = sp.linop.NUFFT(mps.shape,
                           coord=traj,
                           oversamp=1.5,
                           width=3,
                           toeplitz=True)
        if use_dcf:
            D = sp.linop.Multiply(F.oshape, dcf_full**0.5)  # Optional
            A = D * F * S
            b = mvd(np.reshape(data_norm, (1, np.shape(data_norm)
                    [0]))) * dcf_full**0.5
            # b = b/xp.linalg.norm(mvd(data_load) * dcf_full**0.5) # Optional: comment out for quantitative analysis
        else:
            A = F * S
            b = mvd(np.reshape(data_norm, (1, np.shape(
                data_norm)[0])))
            # b = b/xp.linalg.norm(mvd(data_load)) # Optional: comment out for quantitative analysis

        # Normalize data and A operator for full image only (makes image come out at true intensity before dcf)
        LL = sp.app.MaxEig(A.N, dtype=xp.complex64, device=device,
                           max_iter=30, show_pbar=False, leave_pbar=False).run() * 1.01
        A = np.sqrt(1/LL) * A

        # Compute a single x = A.H b operation (i.e. inverse NUFFT)
        img_full[:, :, imslice] = mvc(abs(A.H * b))

        # Key image 1
        F = sp.linop.NUFFT(mps.shape,
                           coord=traj_key1,
                           oversamp=1.5,
                           width=3,
                           toeplitz=True)
        if use_dcf:
            D = sp.linop.Multiply(F.oshape, dcf_key1**0.5)  # Optional
            A = D * F * S
            b = mvd(np.reshape(
                data_key1, (1, np.shape(data_key1)[0]))) * dcf_key1**0.5
        else:
            A = F * S
            b = mvd(np.reshape(data_key1, (1, np.shape(
                data_key1)[0])))

        # Compute a single x = A.H b operation (i.e. inverse NUFFT)
        img_key1[:, :, imslice] = mvc(abs(A.H * b))

        # Key image 2
        F = sp.linop.NUFFT(mps.shape,
                           coord=traj_key2,
                           oversamp=1.5,
                           width=3,
                           toeplitz=True)
        if use_dcf:
            D = sp.linop.Multiply(F.oshape, dcf_key2**0.5)  # Optional
            A = D * F * S
            b = mvd(np.reshape(
                data_key2, (1, np.shape(data_key2)[0]))) * dcf_key2**0.5
        else:
            A = F * S
            b = mvd(np.reshape(data_key2, (1, np.shape(
                data_key2)[0])))

        # Compute a single x = A.H b operation (i.e. inverse NUFFT)
        img_key2[:, :, imslice] = mvc(abs(A.H * b))

# %% Apply pixel shifts and rotations

# Translation and rotation
dx = 0
dy = 0
rot = 2

# Define functions


def ps(x): return np.roll(x, (dx, dy), axis=(0, 1))  # pixelshift
def ir(x): return np.rot90(abs(x), k=rot, axes=(0, 1))  # rotate image


# %% Crop images from oversampled reconstruction
if crop_image:

    # Define center coordinates
    center_x = (field_of_view - 1)//(2*voxel_size*CROP_FAC)
    center_y = (field_of_view - 1)//(2*voxel_size*CROP_FAC)

    # Define crop boundaries
    crop_x = (np.arange(center_x - ((recon_resolution)//2),
                        center_x + ((recon_resolution)//2), 1, dtype=int))
    crop_y = (np.arange(center_y - ((recon_resolution)//2),
                        center_y + ((recon_resolution)//2), 1, dtype=int))

    # Make crop function
    def crp(x): return x[crop_x, :, :][:, crop_y, :]

    img_full = crp(img_full)
    img_key1 = crp(img_key1)
    img_key2 = crp(img_key2)

# %% Export data into .npy
if dynamic == 1:
    img_dynamic1 = img_full
    np.save(location + '/img_dynamic1.npy', img_dynamic1)
    print("img_full saved as img_dynamic1.npy")
else:
    img_dynamic2 = img_full
    np.save(location + '/img_dynamic2.npy', img_dynamic2)
    print("img_full saved as img_dynamic2.npy")

# %% Visualize images

sp.plot.ImagePlot(ir(ps(img_full)), x=0, y=1, z=2,
                  colormap='gray', title="Full image")
sp.plot.ImagePlot(ir(ps(img_key1)), x=0, y=1, z=2,
                  colormap='gray', title="Key image 1")
sp.plot.ImagePlot(ir(ps(img_key2)), x=0, y=1, z=2,
                  colormap='gray', title="Key image 2")

# %% Estimate a signal mask and a noise mask using k-means

# Perform clustering
N_clusters = 2
kmeans = KMeans(n_clusters=N_clusters, random_state=0).fit(np.reshape(
    img_full, (np.shape(np.ravel(img_full))[0], 1)))  # Reshape to (voxelcount,1) shape array
clustered = kmeans.cluster_centers_[kmeans.labels_]

# Binarize
if N_clusters == 2:
    clustered[clustered == np.percentile(clustered, 75)] = 1
    clustered[clustered == np.percentile(clustered, 25)] = 0

# Reshape clusters into images
img_cluster = np.zeros(shape=(N_clusters, np.shape(img_full)[
    0], np.shape(img_full)[1], np.shape(img_full)[2]))
labels = kmeans.labels_
for n in range(N_clusters):
    cluster_temp = []
    for i in range(len(labels)):
        if (labels[i]) == n:
            cluster_temp.append(float(clustered[i]))
        else:
            cluster_temp.append(1)
    img_cluster[n, :, :, :] = np.array(
        cluster_temp).reshape(img_full.shape)

# Calculate norm of each cluster*img_full to see which one is larger (and which one is signal versus noise)
norm_cluster0 = np.linalg.norm(img_cluster[0, :, :, :]*img_full)
norm_cluster1 = np.linalg.norm(img_cluster[1, :, :, :]*img_full)

if norm_cluster0 >= norm_cluster1:
    img_signal = img_cluster[0, :, :, :]
    img_noise = img_cluster[1, :, :, :]
else:
    img_signal = img_cluster[1, :, :, :]
    img_noise = img_cluster[0, :, :, :]

# Erode the noise image
kernel = np.ones((7, 7), np.uint8)
img_noise = cv2.erode(img_noise, kernel, iterations=5)
# Optional: useful for images with lots of noise/random noise slices
# img_signal = cv2.erode(img_signal, kernel=np.ones(
#     (3, 3), np.uint8), iterations=1)
# img_signal = cv2.dilate(img_signal, kernel=np.ones(
#     (3, 3), np.uint8), iterations=2)

# Visualize
sp.plot.ImagePlot(ps(ir(2*img_signal - img_noise)), x=0, y=1, z=2,
                  title="Signal and noise mask, estimated using k-means",
                  vmin=0, vmax=2, colormap="YlGnBu_r")

# Calculate noise values
img_full_noise_std = np.std(img_full[img_noise == 1])
img_full_noise_median = np.median(img_full[img_noise == 1])

# %% Generate a mask to use for correction map

# Load a custom mask if available (ensure the name of the file matches the target)
try:
    img_mask = nib.load(location + '/img_ventilation_mask.nii.gz')
    img_mask = np.array(img_mask.get_fdata())
    img_mask = np.squeeze(img_mask)
    sp.plot.ImagePlot(ir(ps(img_mask)), x=0, y=1, z=2, title="Mask")
except:
    print("Could not load and show the img_ventilation_mask.nii from data folder.")
    img_mask = img_signal
    print("img_signal mask used for image mask.")

# %% Use the Rose Limit (SNR > 4) to threshold noise/signal regions

# Calculate rose limit threshold
rose_limit = img_full_noise_median + 4 * img_full_noise_std

# Summarize using histograms
plt.figure()
plt.hist(np.ravel(img_full[img_noise == 1]), bins=50,
         density=True,  color='blue', alpha=0.8)
plt.hist(np.ravel(img_full[img_mask == 1]), bins=50,
         density=True, color='orange', alpha=0.8)
plt.title("Histogram of img_full")
plt.legend(['Inside noise region', 'Inside mask region'])
plt.axvline(rose_limit, linestyle='dashed', color='black')
plt.xlabel('Signal intensity')
plt.ylabel('Density')
plt.show()

# Make mask correction
img_mask_corrected = np.array(img_mask)
img_mask_corrected[img_full < rose_limit] = 0

# Overwrite masking function


def mask_filter(x): return img_mask_corrected * x


# Plot corrected mask
sp.plot.ImagePlot(ir(ps(img_mask_corrected)), x=0, y=1,
                  z=2, title="Mask - Rose limit corrected")

# %% Calculate flip angle maps

# Ratio of images
ratio = img_key2/img_key1

# Correct for areas where img_key2 > img_key1 (caused by defects, aliasing, breath-hold issues, etc)
ratio[ratio > 1] = 1

# Measure flip angle
fa = np.arccos(ratio**(1/N_spirals_key2))
flip_angle_map_non_selective = np.degrees(fa)
# Force NaN regions to equal 0 (or arbitrary number)
flip_angle_map_non_selective[np.isnan(flip_angle_map_non_selective)] = 0

# Plot
sp.plot.ImagePlot(ir(ps(flip_angle_map_non_selective)), x=0, y=1, z=2,
                  colormap='viridis', vmax=25, vmin=0, title="Flip angle map - non-selective")

sp.plot.ImagePlot(ir(ps(mask_filter(flip_angle_map_non_selective))), x=0, y=1, z=2,
                  colormap='viridis', vmax=25, vmin=0, title="Masked flip angle map - non-selective")


# %% Re-calculate the flip angle map for a slice selective pulse
# plt.style.use('default')

# Load in profile
gz = np.load("pulseshape/gz_philips.npy")
gz /= np.max(gz)

# Find the peak's maximum value and index
max_index = np.argmax(gz)
peak_maximum = gz[max_index]

# Calculate half maximum value
half_maximum = peak_maximum / 2

# Find the indices where data crosses half maximum value
left_idx = np.argmin(np.abs(gz[:max_index] - half_maximum))
right_idx = max_index + np.argmin(np.abs(gz[max_index:] - half_maximum))

# Calculate FWHM, and normalize z to it
fwhm_index = right_idx-left_idx
fwhm = slice_thickness
z_res = fwhm/fwhm_index  # mm occupied by each z point
z = np.linspace(-len(gz)*z_res/2, len(gz) * z_res/2, len(gz))

# Trim the variables for better visualization
gz = gz[abs(z) < 3*slice_thickness/2]
z = z[abs(z) < 3*slice_thickness/2]

# Ideal slice selective case
gz_ideal = np.zeros_like(gz)
gz_ideal[abs(z) < slice_thickness/2] = 1

# Plot gz
fig = plt.figure(figsize=(15, 5))
ax1 = fig.add_subplot(1, 3, 1)
ax1.plot(z, gz, label="Slice-selective", color='y')
ax1.set_title('g(z) slice-select profile')
ax1.set_xlabel('z (mm)')
ax1.set_ylabel('g(z) (a.u.)')

# Signal function with slice gap
n = np.array(range(N_spirals))
n_mesh, z_mesh = np.meshgrid(n, z)
m_0 = 1
M_0 = m_0 * np.ones_like(n_mesh)


def s_n_gapped(gz, fa):

    # Convert FA to radians
    alpha = np.deg2rad(fa)

    # Signal equation
    s_n = np.zeros_like(z_mesh)
    for i in range(z.shape[0]):
        s_n[i, :] = M_0[i, :] * np.sin(alpha * gz[i])
        for j in range(N_spirals):
            s_n[i, j] *= (np.cos(alpha * gz[i]))**j
            s_n += 1e-19

    return s_n


def s_n_contiguous(gz, fa):

    # Convert FA to radians
    alpha = np.deg2rad(fa)

    # Signal equation
    s_n = np.zeros_like(z_mesh)
    for i in range(z.shape[0]):
        s_n[i, :] = M_0[i, :] * np.sin(alpha * gz[i])
        s_n[i, :] *= np.cos(alpha * gz[i-fwhm_index])**N_spirals
        for j in range(N_spirals):
            s_n[i, j] *= (np.cos(alpha * gz[i]))**j
            s_n += 1e-19

    return s_n


# Surface plot
ax2 = fig.add_subplot(1, 3, 2, projection='3d')
ax2.plot_surface(z_mesh, n_mesh, s_n_gapped(gz=gz, fa=flip_angle_applied), cmap=cm.jet,
                 linewidth=0, antialiased=False)
ax2.set_ylim([N_spirals, 0])
ax2.set_zlim([0, np.max(s_n_gapped(gz=gz, fa=flip_angle_applied))])
ax2.set_title("Slice-selective \n (gapped)", size=16)
ax2.set_xlabel("$z$", size=16)
ax2.set_ylabel("$n$", size=16)
ax2.set_zlabel(r"$M_{xy}$", size=16)

ax3 = fig.add_subplot(1, 3, 3, projection='3d')
ax3.plot_surface(z_mesh, n_mesh, s_n_contiguous(gz=gz, fa=flip_angle_applied), cmap=cm.jet,
                 linewidth=0, antialiased=False)
ax3.set_ylim([N_spirals, 0])
ax3.set_zlim([0, np.max(s_n_contiguous(gz=gz, fa=flip_angle_applied))])
ax3.set_title("Slice-selective \n (contiguous)", size=16)
ax3.set_xlabel("$z$", size=16)
ax3.set_ylabel("$n$", size=16)
ax3.set_zlabel(r"$M_{xy}$", size=16)


def s_n_mean_z(s_n):

    # Find mean over z
    s_n_mean_z = np.mean(s_n, axis=0)

    return s_n_mean_z


# Visualize
fig, (ax1) = plt.subplots(figsize=(5, 5),
                          nrows=1, ncols=1)
ax1.plot(n, s_n_mean_z(s_n_gapped(gz=gz, fa=flip_angle_applied)),
         'r', label="Gapped")
ax1.plot(n, s_n_mean_z(s_n_contiguous(gz=gz, fa=flip_angle_applied)),
         'b', label="Contiguous")
ax1.set_title("Signal measurements", size=16)
ax1.set_xlabel(r"Excitation number $n$", size=16)
ax1.set_ylabel(r"Signal $\int s_n(z) dz$", size=16)
ax1.legend()


def fa_meas(s_n, gz, fa):
    S1 = np.mean(s_n(gz=gz, fa=fa)[:, :N_spirals//2])
    S2 = np.mean(s_n(gz=gz, fa=fa)[:, N_spirals//2:])
    att = S2/S1
    fa_meas = np.arccos(att**(2/N_spirals))

    return fa_meas


# Range of flip angles
min_fa = 3
max_fa = 40
fa_range = np.linspace(min_fa, max_fa, 20)
fa_meas_gapped = np.zeros_like(fa_range)
fa_meas_contiguous = np.zeros_like(fa_range)

# Calculate the flip angle using the ideal pulse assumption equation
for i in range(fa_range.shape[0]):
    fa_meas_gapped[i] = np.rad2deg(
        fa_meas(s_n=s_n_gapped, gz=gz, fa=fa_range[i]))
    fa_meas_contiguous[i] = np.rad2deg(
        fa_meas(s_n=s_n_contiguous, gz=gz, fa=fa_range[i]))

# Visualize
fig, (ax1, ax2) = plt.subplots(figsize=(12, 5),
                               nrows=1, ncols=2)
ax1.plot(fa_range, fa_meas_gapped, 'r*--', label=r"Gapped")
ax1.plot(fa_range, fa_meas_contiguous, 'b*--', label=r"Contiguous")
ax1.plot(fa_range, fa_range, 'c*--', label="True flip angle")
ax1.set_title(r"$cos^{-1}(S_2/S_1)^{2/N_s}$ vs. True", size=16)
ax1.set_xlabel("Applied flip angle (degrees)", size=16)
ax1.set_ylabel("Measured flip angle (degrees)", size=16)
ax1.legend(fontsize=14)

# Calculate correction factor
fa_correction_gapped = fa_range/fa_meas_gapped
fa_correction_contiguous = fa_range/fa_meas_contiguous

# Fitting a 2-dimensional polynomial to the data
degree = 6  # Choose the degree of the polynomial
coeffs_gapped = np.polyfit(fa_meas_gapped, fa_correction_gapped, degree)
poly_fa_gapped = np.poly1d(coeffs_gapped)
coeffs_contiguous = np.polyfit(
    fa_meas_contiguous, fa_correction_contiguous, degree)
poly_fa_contiguous = np.poly1d(coeffs_contiguous)

# Generate the x values for the polynomial fit curve
fit_x_gapped = np.linspace(min(fa_meas_gapped), max(fa_meas_gapped), 100)
fit_y_gapped = poly_fa_gapped(fit_x_gapped)
fit_x_contiguous = np.linspace(
    min(fa_meas_contiguous), max(fa_meas_contiguous), 100)
fit_y_contiguous = poly_fa_contiguous(fit_x_contiguous)

# Plot polynomial fit
ax2.plot(fa_meas_gapped, fa_correction_gapped, 'r*',
         label="Gapped")
ax2.plot(fa_meas_contiguous, fa_correction_contiguous, 'b*',
         label="Contiguous")
ax2.plot(fit_x_gapped, fit_y_gapped, 'r', label="Poly. fit - gapped")
ax2.plot(fit_x_contiguous, fit_y_contiguous,
         'b', label="Poly. fit - contiguous")
ax2.legend(fontsize=14)
ax2.set_title(r"$FA_{meas}$ / $FA_{app}$", size=16)
ax2.set_xlabel("Measured flip angle (degrees)", size=16)
ax2.set_ylabel(r"$FA_{meas}$ / $FA_{app}$", size=16)
ax2.text(0.05, 0.5, r"$FA_{app} = \phi \times FA_{meas}$",
                    transform=ax2.transAxes,
                    bbox=dict(facecolor='white', alpha=0.5),
                    size=16)

flip_angle_map = np.zeros_like(flip_angle_map_non_selective)

# Automate correction
flip_angle_map[..., 0] = poly_fa_gapped(flip_angle_map_non_selective[..., 0]) * \
    flip_angle_map_non_selective[..., 0]

flip_angle_map[..., 1:] = poly_fa_contiguous(flip_angle_map_non_selective[..., 1:]) * \
    flip_angle_map_non_selective[..., 1:]

# Limit highest/lowest flip angle
flip_angle_map[flip_angle_map < min_fa] = min_fa
flip_angle_map[flip_angle_map > max_fa] = max_fa

sp.plot.ImagePlot(ir(ps(flip_angle_map)), x=0, y=1, z=2,
                  colormap='viridis', vmax=25, vmin=0, title="Flip angle map - slice-selective")

sp.plot.ImagePlot(ir(ps(mask_filter(flip_angle_map))), x=0, y=1, z=2,
                  colormap='viridis', vmax=25, vmin=0, title="Masked flip angle map - slice-selective")

sp.plot.ImagePlot(ir(ps(flip_angle_map - flip_angle_map_non_selective)), x=0, y=1, z=2,
                  colormap='magma', vmax=10, vmin=0, title="Difference map following FA correction")

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
ni_img = nib.Nifti1Image(img_full, affine=aff)
nib.save(ni_img, location + '/img_ventilation.nii.gz')

ni_img = nib.Nifti1Image(flip_angle_map, affine=aff)
nib.save(ni_img, location + '/flip_angle_map.nii.gz')

ni_img = nib.Nifti1Image(mask_filter(flip_angle_map), affine=aff)
nib.save(ni_img, location + '/flip_angle_map_mask.nii.gz')

# %% Calculate the correction map for a slice selective pulse

# Set up integral functions for the "average" decayed image terms


def integral_gapped(x, gz):
    integrand = np.zeros_like(gz)
    for i in range(gz.shape[0]):
        integrand[i] = np.sin(x * gz[i])
        integrand[i] *= (1 - (np.cos(x * gz[i])) ** N_spirals)
        integrand[i] /= (1 - np.cos(x * gz[i]))

    integral = np.mean(integrand)
    return integral


def integral_contiguous(x, gz):
    integrand = np.zeros_like(gz)
    for i in range(gz.shape[0]):
        integrand[i] = np.sin(x * gz[i]) * np.cos(x *
                                                  gz[i - fwhm_index]) ** N_spirals
        integrand[i] *= (1 - (np.cos(x * gz[i])) ** N_spirals)
        integrand[i] /= (1 - np.cos(x * gz[i]))

    integral = np.mean(integrand)
    return integral


def integral_non_selective(x, gz):
    gz = np.ones_like(gz)
    integrand = np.zeros_like(gz)
    for i in range(gz.shape[0]):
        integrand[i] = np.sin(x * gz[i])
        integrand[i] *= (1 - (np.cos(x * gz[i])) ** N_spirals)
        integrand[i] /= (1 - np.cos(x * gz[i]))

    integral = np.mean(integrand)
    return integral


def correction_factor(x, integral):
    correction = (N_spirals)
    correction *= (1/x)  # RF reciprocity
    xr = np.deg2rad(x)
    correction /= integral(xr, gz)

    return correction


fa_range = np.linspace(min_fa, max_fa, 200)
correction_map_gapped = np.zeros_like(fa_range)
correction_map_contiguous = np.zeros_like(fa_range)
correction_map_non_selective = np.zeros_like(fa_range)

for i in range(fa_range.shape[0]):
    correction_map_contiguous[i] = correction_factor(
        x=fa_range[i],
        integral=integral_contiguous)
    correction_map_gapped[i] = correction_factor(
        x=fa_range[i],
        integral=integral_gapped)
    correction_map_non_selective[i] = correction_factor(
        x=fa_range[i],
        integral=integral_non_selective)

correction_map_contiguous /= correction_factor(
    x=flip_angle_applied,
    integral=integral_contiguous)
correction_map_gapped /= correction_factor(
    x=flip_angle_applied,
    integral=integral_contiguous)
correction_map_non_selective /= correction_factor(
    x=flip_angle_applied,
    integral=integral_non_selective)

# Visualize
fig, (ax1, ax2) = plt.subplots(figsize=(12, 5),
                               nrows=1, ncols=2)
ax1.plot(fa_range, correction_map_contiguous, 'b', label="Contiguous")
ax1.plot(fa_range, correction_map_gapped, 'r', label="Gapped")
ax1.plot(fa_range, correction_map_non_selective, 'c', label="Non-selective")
ax1.set_title("Correction factor", size=16)
ax1.set_xlabel(r"Measured flip angle (degrees)", size=16)
ax1.set_ylabel(r"Correction factor (a.u.)", size=16)
ax1.legend()

# Fitting a 2-dimensional polynomial to the data
degree = 20  # Choose the degree of the polynomial
coeffs_gapped = np.polyfit(fa_range, correction_map_gapped, degree)
poly_cor_gapped = np.poly1d(coeffs_gapped)
coeffs_contiguous = np.polyfit(fa_range, correction_map_contiguous, degree)
poly_cor_contiguous = np.poly1d(coeffs_contiguous)
coeffs_non_selective = np.polyfit(
    fa_range, correction_map_non_selective, degree)
poly_cor_non_selective = np.poly1d(coeffs_non_selective)

# Generate the x values for the polynomial fit curve
fit_y_gapped = poly_cor_gapped(fa_range)
fit_y_contiguous = poly_cor_contiguous(fa_range)
fit_y_non_selective = poly_cor_non_selective(fa_range)

# Plot polynomial fit
ax2.plot(fa_range, correction_map_gapped, 'r*',
         label="Gapped")
ax2.plot(fa_range, correction_map_contiguous, 'b*',
         label="Contiguous")
ax2.plot(fa_range, correction_map_non_selective, 'c*',
         label="Non-selective")
ax2.plot(fa_range, fit_y_gapped, 'r', label="Poly. fit - gapped")
ax2.plot(fa_range, fit_y_contiguous,
         'b', label="Poly. fit - contiguous")
ax2.plot(fa_range, fit_y_non_selective,
         'c', label="Poly. fit - non-selective")
ax2.legend(fontsize=14)
ax2.set_title(r"Correction factor", size=16)
ax2.set_xlabel("Measured flip angle (degrees)", size=16)
ax2.set_ylabel(r"Correction factor", size=16)

# Apply corrections
correction_map = np.zeros_like(flip_angle_map)
correction_map_non_selective = np.zeros_like(
    flip_angle_map)  # Overwrites previous

# Automate correction
correction_map[..., 0] = poly_cor_gapped(flip_angle_map[..., 0])
correction_map[..., 1:] = poly_cor_contiguous(flip_angle_map[..., 1:])
correction_map_non_selective = poly_cor_non_selective(flip_angle_map)

sp.plot.ImagePlot(ir(ps(correction_map)), x=0, y=1, z=2,
                  colormap='magma', vmax=4, vmin=0, title="Correction map - slice-selective")

sp.plot.ImagePlot(ir(ps(mask_filter(correction_map))), x=0, y=1, z=2,
                  colormap='magma', vmax=4, vmin=0, title="Masked correction map - slice-selective")

sp.plot.ImagePlot(ir(ps(mask_filter(correction_map - correction_map_non_selective))), x=0, y=1, z=2,
                  colormap='magma', vmax=4, vmin=-4, title="Correction map (SS) - (NS)")


# %% Apply inpainting to correction map

# Apply inpainting
correction_map_inpainted = correction_map
for imslice in range(N_slices):
    a = mask_filter(correction_map)[:, :, imslice].astype(np.float32)
    b = np.array(1 - img_mask_corrected[:, :, imslice], dtype='uint8')
    inpainted = cv2.inpaint(a, b, 5, cv2.INPAINT_NS)
    correction_map_inpainted[:, :, imslice] = inpainted

# Visualize
sp.plot.ImagePlot(ir(ps(correction_map_inpainted)), x=0, y=1, z=2,
                  colormap='magma', title="Correction map after inpainting", vmin=0, vmax=4)

# %% Apply smoothing to correction map

# Apply smoothing
correction_map_smoothed = correction_map_inpainted

for imslice in range(N_slices):
    correction_map_smoothed[:, :, imslice] = ndimage.gaussian_filter(
        correction_map_inpainted[:, :, imslice], sigma=(4))

sp.plot.ImagePlot(ir(ps(correction_map_smoothed)), x=0, y=1, z=2,
                  colormap='magma', title="Correction map after smoothing", vmin=0, vmax=4)

# %% Apply B-spline fit

number_of_random_points = 10000
correction_map_spline = np.zeros_like(correction_map_inpainted)

for imslice in range(N_slices):
    img = ants.from_numpy(abs(correction_map_inpainted[:, :, imslice]))
    img_array = img.numpy()
    row_indices = np.random.choice(
        range(2, img_array.shape[0]), number_of_random_points)
    col_indices = np.random.choice(
        range(2, img_array.shape[1]), number_of_random_points)

    scattered_data = np.zeros((number_of_random_points, 1))
    parametric_data = np.zeros((number_of_random_points, 2))

    for i in range(number_of_random_points):
        scattered_data[i, 0] = img_array[row_indices[i], col_indices[i]]
        parametric_data[i, 0] = row_indices[i]
        parametric_data[i, 1] = col_indices[i]

    correction_map_spline[:, :, imslice] = ants.fit_bspline_object_to_scattered_data(
        scattered_data, parametric_data,
        parametric_domain_origin=[0.0, 0.0],
        parametric_domain_spacing=[1.0, 1.0],
        parametric_domain_size=img.shape,
        number_of_fitting_levels=3, mesh_size=1).numpy()

sp.plot.ImagePlot(ir(ps(correction_map_spline)), x=0, y=1, z=2,
                  colormap='magma', title="Correction map with B-spline", vmin=0, vmax=4)

# %% Apply correction map to original image

# Matrix multiply
img_corrected = img_full * correction_map_spline

# Visualize
sp.plot.ImagePlot(ir(ps(img_full)), x=0, y=1, z=2,
                  title="Uncorrected image")
sp.plot.ImagePlot(ir(ps(img_corrected)), x=0, y=1, z=2,
                  title="Keyhole corrected image")

# Test difference
img_difference = (img_corrected - img_full) * 100 / img_full
sp.plot.ImagePlot(ir(ps(img_difference)), x=0, y=1, z=2,
                  colormap='magma_r', title="(img_corrected - img_full) / img_full (%)")

# Save images
ni_img = nib.Nifti1Image(img_corrected, affine=aff)
nib.save(ni_img, location + '/img_ventilation_corrected.nii.gz')

# %% Generate an image mask using AI model
try:
    # Load AI model
    model = load_model(
        'maskmodel/Resnet_model_Xe_2D_Vent_2000epochs.hdf5', compile=False)

    # Make the same orientation as data it was trained on
    img_full_model = np.flip(np.rot90(img_full, k=3), axis=1)

    # Normalize between 0 and 1
    img_full_model = (img_full_model - np.min(img_full_model)) / \
        (np.max(img_full_model) - np.min(img_full_model))

    # Generate mask by iterating through each slice
    img_mask_model = np.zeros_like(img_full_model)
    for imslice in range(N_slices):
        img_temp = img_full_model[:, :, imslice]  # Extract slice
        # Resize to match model
        img_temp = resize(img_temp, (256, 256),
                          mode='constant', preserve_range=True)
        img_temp = np.stack((img_temp,)*3, axis=-1)  # Add 3 channels
        img_temp = np.expand_dims(img_temp, 0)  # Add point dimension
        prediction = model.predict(img_temp)
        prediction = prediction[0, :, :, 0]
        prediction = resize(prediction, (recon_resolution,
                            recon_resolution), mode='constant', preserve_range=True)
        img_mask_model[:, :, imslice] = prediction
    img_mask_model = img_mask_model > 0.9

    # Undo the transformation to work with this data
    img_mask_model = np.rot90(np.flip(img_mask_model, axis=1), k=-3)

    # Map Booleans to integers
    img_mask_model = img_mask_model.astype(np.int16)

    # Plot the resnet mask
    sp.plot.ImagePlot(ir(ps(img_mask_model)), x=0, y=1, z=2,
                      title="ResNet mask - Machine Learning")

    # Save images
    ni_img = nib.Nifti1Image(img_mask_model, affine=aff, dtype=np.int16)
    nib.save(ni_img, location + '/img_ventilation_mask_resnet.nii.gz')

except:
    print("Could not save ResNet mask. Check for bugs.")
    img_mask_model = np.zeros_like(img_full_model)

# %% Perform N4 bias field correction

# Reload a mask in case one is already saved and you don't want the 'auto' mask
try:
    img_mask = nib.load(location + '/img_ventilation_mask.nii.gz')
    img_mask = np.array(img_mask.get_fdata())
    img_mask = np.squeeze(img_mask)
    sp.plot.ImagePlot(ir(ps(img_mask)), x=0, y=1, z=2, title="Mask")
except:
    print("Could not load and show the img_ventilation_mask.nii from data folder.")
    img_mask = img_signal
    print("img_signal mask used for image mask.")

# N4 correction
if n4_correction:
    img_n4 = ants.n4_bias_field_correction(image=ants.from_numpy(abs(img_full)),
                                           mask=ants.from_numpy(
        abs(img_mask)),
        shrink_factor=1,
        convergence={
        'iters': [100, 50, 50, 50],
        'tol': 1e-7},
        spline_param=200,
        verbose=False,
        return_bias_field=False)

    # Revert back to numpy
    img_n4 = img_n4.numpy()

    # Visualize
    sp.plot.ImagePlot(ir(ps(img_n4)), x=0, y=1, z=2,
                      title="N4 corrected image")

    # Save images
    ni_img = nib.Nifti1Image(img_n4, affine=aff)
    nib.save(ni_img, location + '/img_ventilation_N4.nii.gz')
else:
    print("N4 bias correction not performed.")
    img_n4 = np.zeros_like(img_mask)

# %% Calculate flip angle uncertainty map - applies to non-selective only

# TODO: Re-calculate/derive using slice selective model

# M. Costa, et. al., ISMRM Proceedings 2021
ratio = img_key2/img_key1
ratio[ratio >= 1] = 1 - 1e-10
flip_angle_map_uncertainty = (ratio**(1/N_spirals_key1)) * np.std(img_key1[img_noise == 1]) / (N_spirals_key1 * img_key1) \
    * np.sqrt((1 + (ratio**(- N_spirals))) / (1 - (ratio**2)))

# Make equal to a percentage of applied flip angle
flip_angle_map_uncertainty = (
    flip_angle_map_uncertainty / flip_angle_applied) * 100 / 5.2
# Note the additional 5.2 factor placed according to Mariah

# Deal with infinites/NaNs and cap off at 100% uncertainty
flip_angle_map_uncertainty[np.isnan(flip_angle_map_uncertainty)] = 100
flip_angle_map_uncertainty[np.isinf(flip_angle_map_uncertainty)] = 100

# Visualize
sp.plot.ImagePlot(ps(ir(mask_filter(flip_angle_map_uncertainty))), x=0, y=1, z=2, colormap='magma',
                  vmax=50, vmin=0, title='Flip angle uncertainty (%) (up to 20%)')

# %% Dummy code for dual dynamic flip angle map calculation
try:
    img_dynamic1 = np.load(location + '/img_dynamic1.npy')
    img_dynamic2 = np.load(location + '/img_dynamic2.npy')

    # Ratio of images
    ratio = img_dynamic2/img_dynamic1

    # Correct for areas where img_key2 > img_key1 (caused by defects, aliasing, breath-hold issues, etc)
    ratio[ratio > 1] = 1

    # Measure flip angle
    fa_paired_image_non_selective = np.rad2deg(
        np.arccos(ratio**(1/(N_spirals))))
    # Force NaN regions to equal 0 (or arbitrary number)
    fa_paired_image_non_selective[np.isnan(fa_paired_image_non_selective)] = 0

    # Plot
    sp.plot.ImagePlot(ir(ps(mask_filter(fa_paired_image_non_selective))), x=0, y=1, z=2,
                      colormap='viridis', vmax=25, vmin=0, title="Flip angle map - dual dynamic - non-selective")


except:
    print('Paired image flip angle map not calculated.')

# Make flip angle corrections for slice selective profile for dual dynamic
# TODO: Clean up as this currently overwrites the keyhole correction variables

# TODO: if possible, can you please write a basic docstring. RH

try:
    # Signal function with slice gap
    n_paired_image = np.array(range(N_spirals*N_dynamics))
    n_mesh_paired_image, z_mesh = np.meshgrid(n_paired_image, z)
    m_0 = 1
    M_0 = m_0 * np.ones_like(n_mesh_paired_image)

    # Update functions for dual dynamic number of spirals

    def s_n_gapped_paired_image(gz, fa):

        # Convert FA to radians
        alpha = np.deg2rad(fa)

        # Signal equation
        s_n = np.zeros_like(z_mesh)
        for i in range(z.shape[0]):
            s_n[i, :] = M_0[i, :] * np.sin(alpha * gz[i])
            for j in range(N_spirals*N_dynamics):
                s_n[i, j] *= (np.cos(alpha * gz[i]))**j
                s_n += 1e-19

        return s_n

    def s_n_contiguous_paired_image(gz, fa):

        # Convert FA to radians
        alpha = np.deg2rad(fa)

        # Signal equation
        s_n = np.zeros_like(z_mesh)
        for i in range(z.shape[0]):
            s_n[i, :] = M_0[i, :] * np.sin(alpha * gz[i])
            s_n[i, :] *= np.cos(alpha * gz[i-fwhm_index]
                                )**(N_spirals*N_dynamics)
            for j in range(N_spirals*N_dynamics):
                s_n[i, j] *= (np.cos(alpha * gz[i]))**j
                s_n += 1e-19

        return s_n

    def fa_meas_paired_image(s_n, gz, fa):
        S1 = np.mean(s_n(gz=gz, fa=fa)[:, :N_spirals])
        S2 = np.mean(s_n(gz=gz, fa=fa)[:, N_spirals:])
        att = S2/S1
        fa_meas = np.arccos(att**(1/N_spirals))

        return fa_meas

    # Range of flip angles
    min_fa = 3
    max_fa = 30
    fa_range = np.linspace(min_fa, max_fa, 30)
    fa_meas_gapped = np.zeros_like(fa_range)
    fa_meas_contiguous = np.zeros_like(fa_range)

    # Calculate the flip angle using the ideal pulse assumption equation
    for i in range(fa_range.shape[0]):
        fa_meas_gapped[i] = np.rad2deg(
            fa_meas_paired_image(s_n=s_n_gapped_paired_image, gz=gz, fa=fa_range[i]))
        fa_meas_contiguous[i] = np.rad2deg(
            fa_meas_paired_image(s_n=s_n_contiguous_paired_image, gz=gz, fa=fa_range[i]))

    # Plot gz
    fig = plt.figure(figsize=(15, 5))
    plt.suptitle("Paired-image signal dynamics", fontsize=16)
    ax1 = fig.add_subplot(1, 3, 1)
    ax1.plot(z, gz, label="Slice-selective", color='y')
    ax1.set_title('g(z) slice-select profile')
    ax1.set_xlabel('z (mm)')
    ax1.set_ylabel('g(z) (a.u.)')

    # Surface plot
    ax2 = fig.add_subplot(1, 3, 2, projection='3d')
    ax2.plot_surface(z_mesh, n_mesh_paired_image, s_n_gapped_paired_image(gz=gz, fa=flip_angle_applied), cmap=cm.jet,
                     linewidth=0, antialiased=False)
    ax2.set_ylim([N_spirals*N_dynamics, 0])
    ax2.set_zlim(
        [0, np.max(s_n_gapped_paired_image(gz=gz, fa=flip_angle_applied))])
    ax2.set_title("Slice-selective \n (gapped)", size=16)
    ax2.set_xlabel("$z$", size=16)
    ax2.set_ylabel("$n$", size=16)
    ax2.set_zlabel(r"$M_{xy}$", size=16)

    ax3 = fig.add_subplot(1, 3, 3, projection='3d')
    ax3.plot_surface(z_mesh, n_mesh_paired_image, s_n_contiguous_paired_image(gz=gz, fa=flip_angle_applied), cmap=cm.jet,
                     linewidth=0, antialiased=False)
    ax3.set_ylim([N_spirals*N_dynamics, 0])
    ax3.set_zlim(
        [0, np.max(s_n_contiguous_paired_image(gz=gz, fa=flip_angle_applied))])
    ax3.set_title("Slice-selective \n (contiguous)", size=16)
    ax3.set_xlabel("$z$", size=16)
    ax3.set_ylabel("$n$", size=16)
    ax3.set_zlabel(r"$M_{xy}$", size=16)

    # Visualize
    fig, (ax1, ax2) = plt.subplots(figsize=(12, 5),
                                   nrows=1, ncols=2)
    plt.suptitle("Paired-image flip angle measurements", fontsize=16)
    ax1.plot(fa_range, fa_meas_gapped, 'r*--', label=r"Gapped")
    ax1.plot(fa_range, fa_meas_contiguous, 'b*--', label=r"Contiguous")
    ax1.plot(fa_range, fa_range, 'c*--', label="True flip angle")
    ax1.set_title(r"$cos^{-1}(S_2/S_1)^{1/N_s}$ vs. True", size=16)
    ax1.set_xlabel("Applied flip angle (degrees)", size=16)
    ax1.set_ylabel("Measured flip angle (degrees)", size=16)
    ax1.legend(fontsize=14)

    # Calculate correction factor
    fa_correction_gapped = fa_range/fa_meas_gapped
    fa_correction_contiguous = fa_range/fa_meas_contiguous

    # Fitting a 2-dimensional polynomial to the data
    degree = 10  # Choose the degree of the polynomial
    coeffs_gapped = np.polyfit(fa_meas_gapped, fa_correction_gapped, degree)
    poly_fa_gapped = np.poly1d(coeffs_gapped)
    coeffs_contiguous = np.polyfit(
        fa_meas_contiguous, fa_correction_contiguous, degree)
    poly_fa_contiguous = np.poly1d(coeffs_contiguous)

    # Generate the x values for the polynomial fit curve
    fit_x_gapped = np.linspace(min(fa_meas_gapped), max(fa_meas_gapped), 100)
    fit_y_gapped = poly_fa_gapped(fit_x_gapped)
    fit_x_contiguous = np.linspace(
        min(fa_meas_contiguous), max(fa_meas_contiguous), 100)
    fit_y_contiguous = poly_fa_contiguous(fit_x_contiguous)

    # Plot polynomial fit
    ax2.plot(fa_meas_gapped, fa_correction_gapped, 'r*',
             label="Gapped")
    ax2.plot(fa_meas_contiguous, fa_correction_contiguous, 'b*',
             label="Contiguous")
    ax2.plot(fit_x_gapped, fit_y_gapped, 'r', label="Poly. fit - gapped")
    ax2.plot(fit_x_contiguous, fit_y_contiguous,
             'b', label="Poly. fit - contiguous")
    ax2.legend(fontsize=14)
    ax2.set_title(r"$FA_{meas}$ / $FA_{app}$", size=16)
    ax2.set_xlabel("Measured flip angle (degrees)", size=16)
    ax2.set_ylabel(r"$FA_{meas}$ / $FA_{app}$", size=16)
    ax2.text(0.05, 0.5, r"$FA_{app} = \phi \times FA_{meas}$",
                        transform=ax2.transAxes,
                        bbox=dict(facecolor='white', alpha=0.5),
                        size=16)

    # Apply flip angle corrections for slice-selective mode
    flip_angle_map_paired_image = np.zeros_like(fa_paired_image_non_selective)

    # Automate correction
    flip_angle_map_paired_image[..., 0] = poly_fa_gapped(fa_paired_image_non_selective[..., 0]) * \
        fa_paired_image_non_selective[..., 0]

    flip_angle_map_paired_image[..., 1:] = poly_fa_contiguous(fa_paired_image_non_selective[..., 1:]) * \
        fa_paired_image_non_selective[..., 1:]

    sp.plot.ImagePlot(ir(ps(mask_filter(flip_angle_map_paired_image))), x=0, y=1, z=2,
                      colormap='viridis', vmax=25, vmin=0, title="Flip angle map - dual dynamic - slice-selective")

    # Limit highest/lowest flip angle
    flip_angle_map_paired_image[flip_angle_map_paired_image < min_fa] = min_fa
    flip_angle_map_paired_image[flip_angle_map_paired_image > max_fa] = max_fa

except:
    if N_dynamics != 1:
        print(
            "Could not perform the slice-selective flip angle correction for paired image.")
    else:
        print("Single dynamic so slice-selective correction not applied to dual dynamic.")


# %% Make some beautiful montages for figure generation


def create_montage(input, slices=None):
    """
    A function that creates a 2D N_pixel x (N_slice * N_pixel) array from a 3D N_pixel x N_pixel x N_slice array. 

    Args:
        input (ndarray): 3D array, in form N_pixel x N_pixel x N_slice.
        slices (list of ints): slices to make into montage. Defaults to middle 5 slices.
    """
    if np.size(np.shape(input)) != 3:
        raise ValueError("Incorrect data dimensions.")

    if slices == None:
        median = np.shape(input)[2]//2
        slices = [median-2, median-1, median, median+1, median+2]

    output = input[:, :, slices[0]]
    for i in slices[1:]:
        temp = input[:, :, i]
        output = np.hstack((output, temp))

    return output


def reorient(input):
    """
    A function that reorients an image array into a plotable object using plt.imshow().

    Args:
        input (ndarray): 3D array, in form N_pixel x N_pixel x N_slice.
    """
    output = np.flip(np.rot90(input), 0)
    return output

# Use custom mask


def mask_filter(x): return img_mask_corrected * x  # Segmented and rose limit
def mask_filter(x): return img_mask_model * x  # AI
# def mask_filter(x): return img_mask * x  # Custom
# def mask_filter(x): return img_mask_model * img_mask * x  # Custom * AI


# Select range of slices to montage

slices = [5, 7, 9, 10, 11, 13, 14]  # Big
slices = [4, 5, 6, 7, 8, 9, 10, 11]  # Medium

# slices = [2, 3, 4, 5, 6, 7]  # Small
# slices = range(3, N_slices-3)  # Auto
slices = range(4, N_slices-4)  # Auto
# slices = (6,)  # Median

img_array_size = (2*len(slices) + 1, 2)

# Perform relevant normalizations
# img_full_norm = img_full / \
#     np.percentile(np.ravel(mask_filter(img_full)), 99.9)
# img_corrected_norm = img_corrected / \
#     np.percentile(np.ravel(mask_filter(img_corrected)), 99.9)
# img_n4_norm = img_n4 / \
#     np.percentile(np.ravel(mask_filter(img_n4)), 99.9)

img_full_norm = img_full / \
    np.percentile(np.ravel(mask_filter(img_corrected)), 99.9)
img_corrected_norm = img_corrected / \
    np.percentile(np.ravel(mask_filter(img_corrected)), 99.9)
img_n4_norm = img_n4 / \
    np.percentile(np.ravel(mask_filter(img_corrected)), 99.9)

# Uncorrected images
img_array = reorient(img_full_norm)
montage = create_montage((img_array), slices)

plt.figure(figsize=img_array_size, dpi=100)
plt.imshow(montage, cmap="gray", vmin=0, vmax=1)
plt.axis('off')
cbar = plt.colorbar(location='right')
cbar.set_ticks(ticks=[0, 0.5, 1])
cbar.ax.set_yticklabels(['0.0', '0.5', '1.0'])
cbar.ax.tick_params(labelsize=16)
plt.show()

# Corrected images
img_array = reorient(img_corrected_norm)
montage = create_montage((img_array), slices)

plt.figure(figsize=img_array_size, dpi=100)
plt.imshow(montage, cmap="gray", vmin=0, vmax=1)
plt.axis('off')
cbar = plt.colorbar(location='right')
cbar.set_ticks(ticks=[0, 0.5, 1])
cbar.ax.set_yticklabels(['0.0', '0.5', '1.0'])
cbar.ax.tick_params(labelsize=16)
plt.show()

# N4 images
if n4_correction:
    img_array = reorient(img_n4_norm)
    montage = create_montage((img_array), slices)

    plt.figure(figsize=img_array_size, dpi=100)
    plt.imshow(montage, cmap="gray", vmin=0, vmax=1)
    plt.axis('off')
    cbar = plt.colorbar(location='right')
    cbar.set_ticks(ticks=[0, 0.5, 1])
    cbar.ax.set_yticklabels(['0.0', '0.5', '1.0'])
    cbar.ax.tick_params(labelsize=16)
    plt.show()
else:
    print("N4 correction was not performed in this reconstruction.")

# Flip angle map
sigma = 4  # Gaussian kernel size for smoothing
flip_angle_map_smoothed = np.ndarray(np.shape(flip_angle_map))
for imslice in range(N_slices):
    flip_angle_map_smoothed[:, :, imslice] = ndimage.gaussian_filter(
        flip_angle_map[:, :, imslice], sigma=(sigma))
img_array = reorient(mask_filter(flip_angle_map_smoothed))
montage = create_montage((img_array), slices)

plt.figure(figsize=img_array_size, dpi=100)
cmap = plt.cm.get_cmap('viridis').copy()
cmap.set_under('black')
fa_vmax = int(np.percentile(
    np.round(flip_angle_map_smoothed[img_mask == 1] / 2) * 2, 95))
fa_vmin = int(np.percentile(
    np.round(flip_angle_map_smoothed[img_mask == 1] / 2) * 2, 0.1))
plt.imshow(montage, cmap=cmap, vmin=fa_vmin, vmax=fa_vmax)
plt.axis('off')
cbar = plt.colorbar(location='right')
cbar.set_ticks(ticks=[fa_vmin, (fa_vmin+fa_vmax)/2, fa_vmax])
cbar.ax.set_yticklabels(
    [str(fa_vmin) + ' ', str((fa_vmin+fa_vmax)//2) + ' ', str(fa_vmax) + ' '])
cbar.ax.tick_params(labelsize=16)
plt.show()

try:
    # Flip angle map paired image
    flip_angle_map_paired_image_smoothed = np.ndarray(np.shape(flip_angle_map))
    for imslice in range(N_slices):
        flip_angle_map_paired_image_smoothed[:, :, imslice] = ndimage.gaussian_filter(
            flip_angle_map_paired_image[:, :, imslice], sigma=(sigma))

    img_array = reorient(mask_filter(flip_angle_map_paired_image_smoothed))
    montage = create_montage((img_array), slices)

    plt.figure(figsize=img_array_size, dpi=100)
    cmap = plt.cm.get_cmap('viridis').copy()
    cmap.set_under('black')
    plt.imshow(montage, cmap=cmap, vmin=fa_vmin, vmax=fa_vmax)
    plt.axis('off')
    cbar = plt.colorbar(location='right')
    cbar.set_ticks(ticks=[fa_vmin, (fa_vmin+fa_vmax)/2, fa_vmax])
    cbar.ax.set_yticklabels(
        [str(fa_vmin) + ' ', str((fa_vmin+fa_vmax)//2) + ' ', str(fa_vmax) + ' '])
    cbar.ax.tick_params(labelsize=16)
    plt.show()

    # Difference flip angle map
    img_array = reorient(mask_filter(
        flip_angle_map_paired_image_smoothed - flip_angle_map))
    # Make noise mask areas very negative for prettier plotting
    img_array_noise = reorient(mask_filter(np.ones(np.shape(flip_angle_map))))
    img_array_noise[img_array_noise == 0] = -100
    img_array_noise[img_array_noise == 1] = 0
    montage = create_montage((img_array+img_array_noise), slices)

    plt.figure(figsize=img_array_size, dpi=100)
    cmap = plt.cm.get_cmap('bwr').copy()
    cmap.set_under('black')
    plt.imshow(montage, cmap=cmap, vmin=-5, vmax=5)
    plt.axis('off')
    cbar = plt.colorbar(location='right')
    cbar.set_ticks(ticks=[-5, 0, 5])
    cbar.ax.set_yticklabels(['-5 ', '0  ', '+5 '])
    cbar.ax.tick_params(labelsize=16)
    plt.show()

    # SSIM difference
    from skimage.metrics import structural_similarity as ssim
    (ssim_score, ssim_diff) = ssim(mask_filter(
        flip_angle_map_paired_image_smoothed), mask_filter(flip_angle_map_smoothed), full=True)

    img_array = reorient(mask_filter(ssim_diff))
    montage = create_montage((img_array), slices)
    print("SSIM paired-image and keyhole: ", np.round(ssim_score, 3))

    plt.figure(figsize=img_array_size, dpi=100)
    cmap = plt.cm.get_cmap('jet').copy()
    cmap.set_under('black')
    plt.imshow(montage, cmap=cmap, vmin=1e-10, vmax=1)
    plt.axis('off')
    cbar = plt.colorbar(location='right')
    cbar.set_ticks(ticks=[0, 0.5, 1])
    cbar.ax.set_yticklabels(['0.0', '0.5', '1.0'])
    cbar.ax.tick_params(labelsize=16)
    plt.show()

    # L2-norm difference
    l2_diff = np.linalg.norm(np.ravel(mask_filter(
        flip_angle_map_paired_image_smoothed - flip_angle_map)), 2)
    print("L2-norm difference paired-image and keyhole: ", np.round(l2_diff, 3))

    # RMSE
    rmse = (np.mean(np.ravel(mask_filter(
        flip_angle_map_paired_image_smoothed - flip_angle_map))**2))**0.5
    print("RMSE paired-image and keyhole: ", np.round(rmse, 3))

except:
    print("Paired image approach was not used in this reconstruction.")

# Flip angle map uncertainty
img_array = reorient(mask_filter(flip_angle_map_uncertainty))
montage = create_montage((img_array), slices)

plt.figure(figsize=img_array_size, dpi=100)
cmap = plt.cm.get_cmap('gnuplot2').copy()
cmap.set_under('black')
plt.imshow(montage, cmap=cmap, vmin=1e-10, vmax=20)
plt.axis('off')
cbar = plt.colorbar(location='right')
cbar.set_ticks(ticks=[0, 10, 20])
cbar.ax.set_yticklabels(['0 ', '10 ', '20 '])
cbar.ax.tick_params(labelsize=16)
plt.show()

# Correction map
img_array = reorient(mask_filter(correction_map))
montage = create_montage((img_array), slices)

plt.figure(figsize=img_array_size, dpi=100)
cmap = plt.cm.get_cmap('magma').copy()
cmap.set_under('black')
cor_vmax = 3.0
# plt.imshow(montage, cmap=cmap, vmin=1/cor_vmax, vmax=cor_vmax)
# plt.axis('off')
# cbar = plt.colorbar(location='right')
# cbar.set_ticks(ticks=[1/cor_vmax, cor_vmax/2, cor_vmax])
# cbar.ax.set_yticklabels([str(1/cor_vmax), str(cor_vmax/2), str(cor_vmax)])
# cbar.ax.tick_params(labelsize=16)
# plt.show()
plt.figure(figsize=img_array_size, dpi=100)
plt.imshow(montage, cmap=cmap, vmin=0, vmax=cor_vmax)
plt.axis('off')
cbar = plt.colorbar(location='right')
cbar.set_ticks(ticks=[0, cor_vmax/2, cor_vmax])
cbar.ax.set_yticklabels([str(0.0), str(cor_vmax/2), str(cor_vmax)])
cbar.ax.tick_params(labelsize=16)
plt.show()

# Correction map spline
img_array = reorient(correction_map_spline)
montage = create_montage((img_array), slices)

plt.figure(figsize=img_array_size, dpi=100)
cmap = plt.cm.get_cmap('magma').copy()
cmap.set_under('black')
# plt.imshow(montage, cmap=cmap, vmin=1/cor_vmax, vmax=cor_vmax)
# plt.axis('off')
# cbar = plt.colorbar(location='right')
# cbar.set_ticks(ticks=[1/cor_vmax, cor_vmax/2, cor_vmax])
# cbar.ax.set_yticklabels([str(1/cor_vmax), str(cor_vmax/2), str(cor_vmax)])
# cbar.ax.tick_params(labelsize=16)
# plt.show()
plt.figure(figsize=img_array_size, dpi=100)
plt.imshow(montage, cmap=cmap, vmin=0, vmax=cor_vmax)
plt.axis('off')
cbar = plt.colorbar(location='right')
cbar.set_ticks(ticks=[0, cor_vmax/2, cor_vmax])
cbar.ax.set_yticklabels([str(0.0), str(cor_vmax/2), str(cor_vmax)])
cbar.ax.tick_params(labelsize=16)
plt.show()

# Correction map with outline
try:
    G = sp.linop.FiniteDifference(ishape=img_array.shape, axes=(0, 1))
    outline = G.N * img_mask_model
    outline_thresh = 0.01
    outline[outline < outline_thresh] = 0
    outline[outline >= outline_thresh] = 1e3
    img_array = reorient(outline)
    montage = create_montage((img_array), slices)

    plt.figure(figsize=img_array_size, dpi=100)
    plt.imshow(montage, cmap="gray", vmin=0, vmax=1)
    plt.axis('off')
    cbar = plt.colorbar(location='right')
    cbar.set_ticks(ticks=[0, 0.5, 1])
    cbar.ax.set_yticklabels(['0.0', '0.5', '1.0'])
    cbar.ax.tick_params(labelsize=16)
    plt.show()

    correction_map_outline = correction_map_spline - outline*correction_map_spline
    img_array = reorient(correction_map_outline)
    montage = create_montage((img_array), slices)

    plt.figure(figsize=img_array_size, dpi=100)
    cmap = plt.cm.get_cmap('magma').copy()
    cmap.set_under('black')
    # plt.imshow(montage, cmap=cmap, vmin=1/cor_vmax, vmax=cor_vmax)
    # plt.axis('off')
    # cbar = plt.colorbar(location='right')
    # cbar.set_ticks(ticks=[1/cor_vmax, cor_vmax/2, cor_vmax])
    # cbar.ax.set_yticklabels([str(1/cor_vmax), str(cor_vmax/2), str(cor_vmax)])
    # cbar.ax.tick_params(labelsize=16)
    # plt.show()
    plt.figure(figsize=img_array_size, dpi=100)
    plt.imshow(montage, cmap=cmap, vmin=0, vmax=cor_vmax)
    plt.axis('off')
    cbar = plt.colorbar(location='right')
    cbar.set_ticks(ticks=[0, cor_vmax/2, cor_vmax])
    cbar.ax.set_yticklabels([str(0.0), str(cor_vmax/2), str(cor_vmax)])
    cbar.ax.tick_params(labelsize=16)
    plt.show()
except:
    print("Could not draw gradient of mask for the correction spline.")

# ResNet mask
try:
    img_array = reorient(img_mask_model)
    montage = create_montage((img_array), slices)

    plt.figure(figsize=img_array_size, dpi=100)
    plt.imshow(montage, cmap="gray", vmin=0, vmax=1)
    plt.axis('off')
    cbar = plt.colorbar(location='right')
    cbar.set_ticks(ticks=[0, 0.5, 1])
    cbar.ax.set_yticklabels(['0.0', '0.5', '1.0'])
    cbar.ax.tick_params(labelsize=16)
    plt.show()
except:
    print("ResNet mask not plotted.")
# %%
