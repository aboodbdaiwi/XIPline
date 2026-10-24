# -*- coding: utf-8 -*-
"""
Title: 
    spiral_2d_single_channel.py

Agenda: 
    A .py script that uses compressed sensing methods to reconstruct 
non-cartesian Philips data in raw-lab-sin/data-list format. 

Author: 
    Joseph Plummer - joseph.plummer@cchmc.org 

Creation date: 
    2022-07-06
    
Modification date: 
    2022-07-06 - JWP initial commit
    
"""

# %% Import packages
from scipy.io import savemat
import cv2
from sklearn.cluster import KMeans
if __name__ == "__main__":
    # Base
    import numpy as np
    import matplotlib.pyplot as plt
    plt.style.use('dark_background')

    # Readphilips node - confidential
    import readphilips as rp
    from readphilips.file_io import io

    # Sigpy functionality
    import sigpy as sp
    import sigpy.plot as pl
    import sigpy.mri as mr

    # Compressed sensing
    from functions import tv, convexalg, util

    # Nifti functionality
    import nibabel as nib

    # Save results
    import os

    # %% Reconstruction settings

    # Device settings
    devnum = 0
    device = sp.Device(devnum)
    xp = device.xp

    # Density compensation
    use_dcf = False

    # Reconstruct onto larger volume and crop image
    crop_image = True

    # Offset correction
    freq_offset = 0  # Hz

    # Consensus settings
    if use_dcf:
        num_iters = 10
        lamda_1 = 5e-5
        lamda_2 = 10e-1  # 5e-2
        rho = 1
        ptol = 1e-1
        num_normal = 4
    else:
        num_iters = 30
        lamda_1 = 1e-7
        lamda_2 = 1e-5
        rho = 5e1
        ptol = 1e-2
        num_normal = 4

    # %% Import data and trajectories

    # Load data
    location = io(
        initial_dir="C:\\Users\\PLUMD2\\Desktop\\Diffusion Images and Analysis").selectDirectory()
    # location = io(initial_dir = "I:\\Woods_CPIR_Images\\IRC803_DUKE_ILD").selectDirectory()
    filename = io(initial_dir=location, f_type="raw-list files",
                  ext_type="*.data").selectFile()
    inputfile = rp.PhilipsData(filename)
    inputfile.compute()
    data_load = np.array(inputfile.data)
    data_load = np.squeeze(data_load)
    order = inputfile.header.get('list').get('ky').astype(int)
    plot_slice = int(np.shape(data_load)[1]/2)
    fid = abs(data_load[0, plot_slice, 0, :])
    plt.figure()
    plt.plot(fid)
    plt.ylabel('k-space intensity')
    plt.xlabel('Sample number')
    plt.title('Free induction decay')
    plt.show()

    # Load traj
    filename = io(initial_dir=location, f_type="raw-lab-sin files",
                  ext_type="*.sin").selectFile()
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
    plt.title('Spiral coordinates for single slice')
    plt.axis('square')
    plt.show()

    # Understand data
    N_bvalue = np.shape(data_load)[0]
    N_slices = np.shape(data_load)[1]
    N_spirals = np.shape(data_load)[2]
    N_samp = np.shape(data_load)[3]
    order = order[::N_bvalue]  # Deal with repeated ordering for b-value
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
    recon_resolution = int(1.5 * recon_resolution)
    field_of_view = voxel_size * recon_resolution
    scan_resolution = int(inputfile.header.get(
        'sin').get('scan_resolutions')[0][0])
    location_center_coordinates = inputfile.header.get(
        'sin').get('location_center_coordinates')[0]
    loc_ap_rl_fh_offcentres = inputfile.header.get(
        'sin').get('loc_ap_rl_fh_offcentres')[0]
    readout_duration = float(inputfile.header.get(
        'sin').get('sample_time_interval')[0][0])
    flip_angle_applied = float(
        inputfile.header.get('sin').get('flip_angles')[0][0])

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

    # %% Correct for frequency offset
    dwell_time = readout_duration*1e-3/N_samp
    for ii in range(N_samp):
        ii = int(ii)
        scale_factor = np.exp(-2*np.pi*1j*ii*freq_offset*dwell_time)
        data_load[:, :, :, ii] = scale_factor*data_load[:, :, :, ii]

    # %% Reorder data and trajectories to match Philips ordering

    coords_reordered = np.zeros(np.shape(coords_expanded))
    data_reordered = np.zeros(np.shape(data_load), dtype='complex_')
    for jj in range(N_bvalue):
        for ii in range(N_spirals):
            coords_reordered[:, :, ii] = coords_expanded[:, :, order[ii]]
            data_reordered[jj, :, ii, :] = data_load[jj, :, order[ii], :]

    # Confirm reordering
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
    for ii in range(N_bvalue):
        data_reordered[ii, :, :, :] = translate_factor * \
            data_reordered[ii, :, :, :]

    # %% Analyze and manipulate the data

    # Select only a subset of trajectories and data
    # subset = 9
    subset = N_spirals
    data_subset = data_reordered[:, :, :subset, :]
    traj_subset = traj[:subset, ...]

    # Plot trajectories:
    plt.figure()
    N_visual = np.shape(traj_subset)[0]
    color = iter(plt.cm.viridis(np.linspace(0, 1, N_visual)))
    for i in np.linspace(0, N_visual-1, N_visual):
        i = int(i)
        c = next(color)
        plt.scatter(traj_subset[i, :, 0],
                    traj_subset[i, :, 1], color=c, s=1, marker='.')
    plt.ylabel('$k_y$', fontsize=18)
    plt.xlabel('$k_x$', fontsize=18)
    # plt.title('Spiral coordinates for single slice')
    plt.xticks([])
    plt.yticks([])
    plt.axis('square')
    plt.show()

    # Density compensation
    beta = 8  # default = 8
    width = 3  # default = 4
    dcf = mr.pipe_menon_dcf(traj_subset, img_shape=ImageShape,
                            beta=beta, width=width, device=device)
    dcf[..., -N_spirals:] = np.median(dcf)
    fig = plt.figure()
    plt.plot(sp.to_device(dcf[0, ...], -1))
    plt.xlabel('Sample number')
    plt.ylabel('DCF')
    plt.suptitle('Sample Density Compensation Function')
    plt.title(r'$ \beta $ = %.2g, width = %.2g' % (beta, width))
    dcf = np.reshape(dcf, (1, np.shape(dcf)[0], np.shape(dcf)[1]))

    # Initialize b-value image
    img_bvalue_nufft = np.zeros((int(resize_factor * scan_resolution), int(
        resize_factor * scan_resolution), N_slices, N_bvalue), dtype='complex64')
    img_bvalue_cg = np.zeros((int(resize_factor * scan_resolution), int(
        resize_factor * scan_resolution), N_slices, N_bvalue), dtype='complex64')
    img_bvalue_cs = np.zeros((int(resize_factor * scan_resolution), int(
        resize_factor * scan_resolution), N_slices, N_bvalue), dtype='complex64')

    # Loop through each b-value
    for bvalue in range(N_bvalue):
        data_bvalue = data_subset[bvalue, :, :, :]

        # Add an extra dimension for channel
        data_multi_channel = np.reshape(data_bvalue, (1, np.shape(
            data_bvalue)[0], np.shape(data_bvalue)[1], np.shape(data_bvalue)[2]))

        # Estimate sensitivity map - assume all ones for single channel
        mps = np.ones((1, ImageSize, ImageSize), dtype=int)

        # Set device
        with device:
            def mvd(x): return sp.to_device(x, device)
            def mvc(x): return sp.to_device(x, sp.cpu_device)

            # Move to GPU
            data = mvd(data_multi_channel)
            data_subset = mvd(data_subset)
            traj = mvd(traj_subset)
            dcf = mvd(dcf)
            mps = mvd(mps)

            # Compute norm
            data_norm_dcf = xp.linalg.norm(data_subset * dcf**0.5)
            data_norm = xp.linalg.norm(data_subset)

            # Compute linear operators
            S = sp.linop.Multiply(ImageShape, mps)
            F = sp.linop.NUFFT(mps.shape,
                               coord=traj,
                               oversamp=1.25,
                               width=4,
                               toeplitz=True)
            D = sp.linop.Multiply(F.oshape, dcf**0.5)  # Optional

            if use_dcf:
                A = D * F * S
            else:
                A = F * S

            # Calculate maximum eigenvalue
            LL = sp.app.MaxEig(A.N, dtype=xp.complex64,
                               device=device).run() * 1.01

            # Divide by svd eigenvalue to ensure maximum eigenvalue falls between [0,1]
            A = np.sqrt(1/LL) * A

            # Calculate a DCF A for the img_nufft
            A_dcf = D * F * S
            LL_dcf = sp.app.MaxEig(A.N, dtype=xp.complex64,
                                   device=device).run() * 1.01
            A_dcf = np.sqrt(1/LL_dcf) * A_dcf

            # Define regularizing linear operators
            W1 = sp.linop.Wavelet(S.ishape, wave_name="db4")
            def g1(x): return lamda_1 * xp.linalg.norm(W1(x).ravel(), ord=1)
            prox_g1 = sp.prox.UnitaryTransform(
                sp.prox.L1Reg(W1.oshape, lamda_1), W1)

            prox_g2 = tv.ProxTV(A.ishape, lamda_2)
            def g2(x): return lamda_2 * xp.linalg.norm(prox_g2.G(x))

            lst_g = [g1, g2]
            lst_proxg = [prox_g1, prox_g2]

            # Initialize data storage

            img_nufft = np.zeros((ImageSize, ImageSize, N_slices))
            img_cg = np.zeros((ImageSize, ImageSize, N_slices))
            img_cs_cg = np.zeros((ImageSize, ImageSize, N_slices))

            # Loop through slices and reconstruct

            for imslice in range(N_slices):
                data_slice = data[:, imslice, :, :]

                # Correct for sample density and normalize
                if use_dcf:
                    b = data_slice * dcf**0.5
                    b = b/xp.linalg.norm(data_norm_dcf)
                else:
                    b = data_slice
                    b = b/xp.linalg.norm(data_norm)

                # Compute a single x = A.H b operation (i.e. inverse NUFFT)
                b_dcf = data_slice * dcf**0.5 / \
                    xp.linalg.norm(data_norm_dcf * dcf**0.5)
                img_nufft[:, :, imslice] = mvc(abs(A_dcf.H * b_dcf))

                # Solve ||Ax - b||^2 + lamda||x||^2 using conjugate gradient descent
                img_cg[:, :, imslice] = mvc(abs(convexalg.cg(num_iters=40,
                                            ptol=5e-1,
                                            A=A, b=b,
                                            lamda=0,  # l2-norm constraint for additional regularization
                                            verbose=True,
                                            draw_output=False,
                                            save=None)))

                # Solve ||Ax - b||^2 + lamda_1 |Wx| + lamda_2 |Gx| using ADMM
                img_cs_cg[:, :, imslice] = mvc(abs(convexalg.admm(num_iters=num_iters,
                                                                  ptol=ptol,
                                                                  A=A, b=b,
                                                                  num_normal=num_normal,
                                                                  lst_proxg=lst_proxg,
                                                                  rho=rho,
                                                                  lst_g=lst_g,
                                                                  method="cg",
                                                                  verbose=True,
                                                                  draw_output=False)))
            #  Manipulate images for plotting

            # Translation and rotation
            dx = 0
            dy = 0
            rot = 1

            img_nufft = np.flipud(util.ps(util.ir(img_nufft, rot), dx, dy))
            img_cg = np.flipud(util.ps(util.ir(img_cg, rot), dx, dy))
            img_cs_cg = np.flipud(util.ps(util.ir(img_cs_cg, rot), dx, dy))

            # Crop images from oversampled reconstruction
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

                # Crop
                img_nufft = crp(img_nufft)
                img_cg = crp(img_cg)
                img_cs_cg = crp(img_cs_cg)

            # Store an image of your choice for the bvalue image
            img_bvalue_nufft[:, :, :, bvalue] = img_nufft
            img_bvalue_cg[:, :, :, bvalue] = img_cg
            img_bvalue_cs[:, :, :, bvalue] = img_cs_cg

# %% Plot images

fig = plt.figure()
visual_slice = N_slices//2
pl.ImagePlot(abs(img_bvalue_nufft[:, :, visual_slice, [
             2, 3, 0, 1]]), x=0, y=1, z=2,  hide_axes=True)
pl.ImagePlot(abs(img_bvalue_cg[:, :, visual_slice, [
             2, 3, 0, 1]]), x=0, y=1, z=2, hide_axes=True)
pl.ImagePlot(abs(img_bvalue_cs[:, :, visual_slice, [
             2, 3, 0, 1]]), x=0, y=1, z=2, hide_axes=True)

# %% Apply global flip angle correction

# Calculate flip angle correction for each b-value
flip_angle_applied_rad = np.radians(flip_angle_applied)

# Intialize
correction_factor = np.ones((N_bvalue, 1)) / flip_angle_applied_rad
for b_value in range(N_bvalue):
    correction_factor[b_value] *= N_spirals
    correction_factor[b_value] /= np.sin(flip_angle_applied_rad)
    correction_factor[b_value] /= np.cos(flip_angle_applied_rad)**(b_value)
    correction_factor[b_value] *= (1 -
                                   np.cos(flip_angle_applied_rad)**(N_bvalue))
    correction_factor[b_value] /= (1 -
                                   np.cos(flip_angle_applied_rad)**(N_bvalue*N_spirals))


# %% Image segmentation (automatic)

# Extract first b-value (where signal is highest)
img_bvalue_0th = abs(img_bvalue_cs[:, :, :, 0])

# Perform clustering
N_clusters = 2
kmeans = KMeans(n_clusters=N_clusters, random_state=0).fit(np.reshape(
    img_bvalue_0th, (np.shape(np.ravel(img_bvalue_0th))[0], 1)))  # Reshape to (voxelcount,1) shape array
clustered = kmeans.cluster_centers_[kmeans.labels_]

# Binarize
if N_clusters == 2:
    clustered[clustered == np.percentile(clustered, 75)] = 1
    clustered[clustered == np.percentile(clustered, 25)] = 0

# Reshape clusters into images
img_cluster = np.zeros(shape=(N_clusters, np.shape(img_bvalue_0th)[
                       0], np.shape(img_bvalue_0th)[1], np.shape(img_bvalue_0th)[2]))
labels = kmeans.labels_
for n in range(N_clusters):
    cluster_temp = []
    for i in range(len(labels)):
        if (labels[i]) == n:
            cluster_temp.append(float(clustered[i]))
        else:
            cluster_temp.append(1)
    img_cluster[n, :, :, :] = np.array(
        cluster_temp).reshape(img_bvalue_0th.shape)

# Calculate norm of each cluster*img_full to see which one is larger (and which one is signal versus noise)
norm_cluster0 = np.linalg.norm(img_cluster[0, :, :, :]*img_bvalue_0th)
norm_cluster1 = np.linalg.norm(img_cluster[1, :, :, :]*img_bvalue_0th)

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
img_signal = cv2.dilate(img_signal, kernel, iterations=2)

# Visualize
sp.plot.ImagePlot(2*img_signal - img_noise, x=0, y=1, z=2,
                  title="Signal and noise mask, estimated using k-means",
                  vmin=0, vmax=2, colormap="YlGnBu_r")

# %% Re-run the same segmentation without the airway

img_bvalue_3rd = abs(img_bvalue_cs[:, :, :, 3])

# Perform clusting
N_clusters = 2
kmeans = KMeans(n_clusters=N_clusters, random_state=0).fit(np.reshape(
    img_bvalue_3rd, (np.shape(np.ravel(img_bvalue_3rd))[0], 1)))  # Reshape to (voxelcount,1) shape array
clustered = kmeans.cluster_centers_[kmeans.labels_]

# Binarize
if N_clusters == 2:
    clustered[clustered == np.percentile(clustered, 75)] = 1
    clustered[clustered == np.percentile(clustered, 25)] = 0

# Reshape clusters into images
img_cluster = np.zeros(shape=(N_clusters, np.shape(img_bvalue_3rd)[
                       0], np.shape(img_bvalue_3rd)[1], np.shape(img_bvalue_3rd)[2]))
labels = kmeans.labels_
for n in range(N_clusters):
    cluster_temp = []
    for i in range(len(labels)):
        if (labels[i]) == n:
            cluster_temp.append(float(clustered[i]))
        else:
            cluster_temp.append(1)
    img_cluster[n, :, :, :] = np.array(
        cluster_temp).reshape(img_bvalue_3rd.shape)

# Calculate norm of each cluster*img_full to see which one is larger (and which one is signal versus noise)
norm_cluster0 = np.linalg.norm(img_cluster[0, :, :, :]*img_bvalue_3rd)
norm_cluster1 = np.linalg.norm(img_cluster[1, :, :, :]*img_bvalue_3rd)

if norm_cluster0 >= norm_cluster1:
    img_signal_no_airway = img_cluster[0, :, :, :]
    img_noise_no_airway = img_cluster[1, :, :, :]
else:
    img_signal_no_airway = img_cluster[1, :, :, :]
    img_noise_no_airway = img_cluster[0, :, :, :]

# Binarize noise mask
img_noise_no_airway[img_noise_no_airway != 1] = 0

# Erode the noise image
kernel = np.ones((7, 7), np.uint8)
img_noise_no_airway = cv2.erode(img_noise_no_airway, kernel, iterations=5)

# Open the signal image (optional)
kernel = np.ones((3, 3), np.uint8)
img_signal_no_airway = cv2.erode(img_signal_no_airway, kernel, iterations=1)
img_signal_no_airway = cv2.dilate(img_signal_no_airway, kernel, iterations=2)

# Visualize
sp.plot.ImagePlot(2*img_signal_no_airway - img_noise_no_airway, x=0, y=1, z=2,
                  title="Signal and noise mask, estimated using k-means",
                  vmin=0, vmax=2, colormap="YlGnBu_r")

# Alternatively, load a mask
try:
    img_signal = nib.load(location + '/img_mask_manual.nii.gz')
    img_signal = np.array(img_signal.get_data())
    # drop singleton dim in case temporal dim exists
    img_signal = np.squeeze(img_signal)
    # img_signal = ir(img_signal)
    fig = plt.figure()
    pl.ImagePlot(abs(img_signal), x=0, y=1, z=2, title="Manual image mask")
except:
    print("Could not load and show the img_mask_manual.nii from results folder.")

# %% Calculate ADC

# Initialize settings
ADC_map_nufft = np.zeros((recon_resolution, recon_resolution, N_slices))
ADC_map_cs = np.zeros((recon_resolution, recon_resolution, N_slices))
b_values = np.array([0, 10, 20, 30])

# Iterate through image data to calculate ADC
for i in range(recon_resolution):
    for j in range(recon_resolution):
        for k in range(N_slices):
            signal = abs(img_bvalue_nufft[i, j, k, :])/correction_factor[:, 0]
            log_signal = np.log(signal)
            model = np.polyfit(b_values, log_signal, 1)
            slope = model[0]
            ADC_map_nufft[i, j, k] = - slope
            signal = abs(img_bvalue_cs[i, j, k, :])/correction_factor[:, 0]
            log_signal = np.log(signal)
            model = np.polyfit(b_values, log_signal, 1)
            slope = model[0]
            ADC_map_cs[i, j, k] = - slope

# Cleanup
ADC_map_nufft[ADC_map_nufft < 0] = 0
ADC_map_nufft[ADC_map_nufft > 0.14] = 0
ADC_map_cs[ADC_map_cs < 0] = 0
ADC_map_cs[ADC_map_cs > 0.14] = 0

# Statistics
ADC_mean_nufft = np.nanmean(
    ADC_map_nufft[img_signal_no_airway*ADC_map_nufft != 0])
ADC_mean_cs = np.nanmean(ADC_map_cs[img_signal_no_airway*ADC_map_cs != 0])

# Mask for plotting
ADC_map_nufft[img_signal*ADC_map_nufft == 0] = np.nan
ADC_map_cs[img_signal*ADC_map_cs == 0] = np.nan
# ADC_map_nufft[img_signal_no_airway*ADC_map_nufft == 0] = np.nan
# ADC_map_cs[img_signal_no_airway*ADC_map_cs == 0] = np.nan

# Plot
pl.ImagePlot(abs(ADC_map_nufft[:, :, 2:8]), x=0, y=1, z=2, title="NUFFT: mean ADC = %3.4f" % (
    ADC_mean_nufft), colormap="jet", vmax=0.06, vmin=0)
pl.ImagePlot(abs(ADC_map_cs[:, :, 2:8]), x=0, y=1, z=2, title="CS: mean ADC = %3.4f" % (
    ADC_mean_cs), colormap="jet", vmax=0.06, vmin=0)
pl.ImagePlot(abs(ADC_map_nufft[:, :, 2:8]), x=0, y=1, z=2, title="NUFFT: mean ADC = %3.4f" % (
    ADC_mean_nufft), colormap="jet", vmax=0.06, vmin=0, hide_axes=True)
pl.ImagePlot(abs(ADC_map_cs[:, :, 2:8]), x=0, y=1, z=2, title="CS: mean ADC = %3.4f" % (
    ADC_mean_cs), colormap="jet", vmax=0.06, vmin=0, hide_axes=True)

# Summarize using histograms
plt.figure()
plt.hist(np.ravel(ADC_map_nufft), bins=100,
         density=True,  color='blue', alpha=0.5)
plt.hist(np.ravel(ADC_map_cs), bins=100, density=True,  color='red', alpha=0.5)
plt.title("ADC histograms")
plt.legend(['nufft', 'cs'])
plt.xlabel('ADC value')
plt.ylabel('Density')
plt.show()

# %% Make beautiful figure

# Select slice
plot_slice = int(N_slices*0.5)

# Generate images
ADC_map_nufft_slice = abs(ADC_map_nufft[:, :, plot_slice])
ADC_map_cs_slice = abs(ADC_map_cs[:, :, plot_slice])

# Plot b-values
fig, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, sharey=True,
                                         gridspec_kw={'wspace': 0, 'hspace': 0}, squeeze=True)
img_bvalue_nufft_norm = img_bvalue_nufft / \
    np.percentile(
        abs(np.ravel(img_bvalue_nufft[:, :, int(plot_slice), 0])), 99.5)
img_plot1 = ax1.imshow(abs(np.rot90(img_bvalue_nufft_norm[:, :, int(
    plot_slice), 0], 3)), cmap="gray", vmin=0, vmax=1)
img_plot2 = ax2.imshow(abs(np.rot90(img_bvalue_nufft_norm[:, :, int(
    plot_slice), 1], 3)), cmap="gray", vmin=0, vmax=1)
img_plot3 = ax3.imshow(abs(np.rot90(img_bvalue_nufft_norm[:, :, int(
    plot_slice), 2], 3)), cmap="gray", vmin=0, vmax=1)
img_plot4 = ax4.imshow(abs(np.rot90(img_bvalue_nufft_norm[:, :, int(
    plot_slice), 3], 3)), cmap="gray", vmin=0, vmax=1)
ax1.axis('off')
ax2.axis('off')
ax3.axis('off')
ax4.axis('off')
fig.set_size_inches(11, 4, forward=True)
fig.set_dpi(100)
cax = fig.add_axes([ax4.get_position().x1+0.01,
                   ax4.get_position().y0, 0.02, ax4.get_position().height])
cb = plt.colorbar(img_plot4, cax=cax)
cb.outline.set_edgecolor('black')
cb.ax.tick_params(labelsize=16)
cb.set_ticks([0, 0.5, 1])
cb.set_ticklabels([0.0, 0.5, 1.0])

# Plot b-values
fig, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4, sharey=True,
                                         gridspec_kw={'wspace': 0, 'hspace': 0}, squeeze=True)
img_bvalue_cs_norm = img_bvalue_cs / \
    np.percentile(abs(np.ravel(img_bvalue_cs[:, :, int(plot_slice), 0])), 99.5)
img_plot1 = ax1.imshow(abs(np.rot90(
    img_bvalue_cs_norm[:, :, int(plot_slice), 0], 3)), cmap="gray", vmin=0, vmax=1)
img_plot2 = ax2.imshow(abs(np.rot90(
    img_bvalue_cs_norm[:, :, int(plot_slice), 1], 3)), cmap="gray", vmin=0, vmax=1)
img_plot3 = ax3.imshow(abs(np.rot90(
    img_bvalue_cs_norm[:, :, int(plot_slice), 2], 3)), cmap="gray", vmin=0, vmax=1)
img_plot4 = ax4.imshow(abs(np.rot90(
    img_bvalue_cs_norm[:, :, int(plot_slice), 3], 3)), cmap="gray", vmin=0, vmax=1)
ax1.axis('off')
ax2.axis('off')
ax3.axis('off')
ax4.axis('off')
fig.set_size_inches(11, 4, forward=True)
fig.set_dpi(100)
cax = fig.add_axes([ax4.get_position().x1+0.01,
                   ax4.get_position().y0, 0.02, ax4.get_position().height])
cb = plt.colorbar(img_plot4, cax=cax)
cb.outline.set_edgecolor('black')
cb.ax.tick_params(labelsize=16)
cb.set_ticks([0, 0.5, 1])
cb.set_ticklabels([0.0, 0.5, 1.0])

# Plot ADC
fig, (ax1) = plt.subplots(1, 1, sharey=True, gridspec_kw={
    'wspace': 0, 'hspace': 0}, squeeze=True)
img_plot1 = ax1.imshow(
    np.rot90(ADC_map_nufft_slice[:, :], 3), cmap="jet", vmax=0.06)
ax1.axis('off')
fig.set_size_inches(4, 4, forward=True)
fig.set_dpi(100)
cax = fig.add_axes([ax1.get_position().x1+0.01,
                   ax1.get_position().y0, 0.02, ax1.get_position().height])
cb = plt.colorbar(img_plot1, cax=cax)
cb.outline.set_edgecolor('black')
cb.ax.tick_params(labelsize=16)
cb.set_ticks([0, 0.03, 0.06])
cb.set_ticklabels([0.0, 0.03, 0.06])

# Plot ADC
fig, (ax1) = plt.subplots(1, 1, sharey=True, gridspec_kw={
    'wspace': 0, 'hspace': 0}, squeeze=True)
img_plot1 = ax1.imshow(
    np.rot90(ADC_map_cs_slice[:, :], 3), cmap="jet", vmax=0.06)
ax1.axis('off')
fig.set_size_inches(4, 4, forward=True)
fig.set_dpi(100)
cax = fig.add_axes([ax1.get_position().x1+0.01,
                   ax1.get_position().y0, 0.02, ax1.get_position().height])
cb = plt.colorbar(img_plot1, cax=cax)
cb.outline.set_edgecolor('black')
cb.ax.tick_params(labelsize=16)
cb.set_ticks([0, 0.03, 0.06])
cb.set_ticklabels([0.0, 0.03, 0.06])

# Generate images
ADC_map_cs_slice = abs(ADC_map_cs[:, :, plot_slice])

# Plot
plt.figure()
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, sharey=True, gridspec_kw={'wspace': 0, 'hspace': 0},
                                    squeeze=True)
img_plot1 = ax1.imshow(
    np.rot90(ADC_map_cs[:, :, int(N_slices*0.4)], 3), cmap="jet", vmax=0.06)
img_plot2 = ax2.imshow(
    np.rot90(ADC_map_cs[:, :, int(N_slices*0.5)], 3), cmap="jet", vmax=0.06)
img_plot3 = ax3.imshow(
    np.rot90(ADC_map_cs[:, :, int(N_slices*0.6)], 3), cmap="jet", vmax=0.06)
ax1.axis('off')
ax2.axis('off')
ax3.axis('off')

# %% Save images as Nifti files


# Save matlab files
savemat(location + "\\img_diffusion_nufft.mat",
        mdict={'img_diffusion_nufft': img_bvalue_nufft})
savemat(location + "\\img_diffusion_cs.mat",
        mdict={'img_diffusion_cs': img_bvalue_cs})
savemat(location + "\\img_diffusion_mask_auto.mat",
        mdict={'img_diffusion_mask_auto': img_signal})
savemat(location + "\\img_diffusion_mask_airway_rm_auto.mat",
        mdict={'img_diffusion_mask_airway_rm_auto': img_signal_no_airway})
savemat(location + "\\img_diffusion_noise_mask_auto.mat",
        mdict={'img_diffusion_noise_mask_auto': img_noise})
savemat(location + "\\k_space_ordered.mat",
        mdict={'k_space_ordered': mvc(data_subset)})
savemat(location + "\\traj_ordered.mat",
        mdict={'traj_ordered': mvc(traj_subset)})

# Save the result as a Nifti
ni_img = nib.Nifti1Image(abs(img_bvalue_nufft), affine=aff)
nib.save(ni_img, location + '\\img_bvalue_nufft')

ni_img = nib.Nifti1Image(abs(img_bvalue_cs), affine=aff)
nib.save(ni_img, location + '\\img_bvalue_cs')
# %%
