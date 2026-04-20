

# %% Import packages
''' '''
#import cupy as cp
import nibabel as nib
import readphilips as rp
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('dark_background')
import sigpy.mri as mr
import sigpy as sp
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

# Offset correction
#freq_offset = 150  # Hz

# Reconstruct onto larger volume and crop image
crop_image = True

# Select dynamic
dynamic = 1

# %% Read data


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

    return datafile_r, trajfile_r, datalocation_r, freq_offset

config_Xe_path=r"C:\XIPline\offline_recon\config_Xe.txt"
datafile, trajfile, datalocation, freq_offset = read_config_and_find_files(config_Xe_path)
print("Data file:", datafile)
print("Trajectory file:", trajfile)
print("Data location:", datalocation)
print("Frequency offset:", freq_offset)

inputfile = rp.PhilipsData(datafile)
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
if debug_mod:
    plt.figure()
    plt.plot(fid)
    plt.ylabel('k-space intensity')
    plt.xlabel('Sample number')
    plt.title('Free induction decay')
    plt.show()

# Load traj
#filename = r"\\rds6.cchmc.org\PulMed-43\CPIR_Share\Carter\ForAbood_10Jun2025_IRC1031-140c\offline_recon_asb\gas\20250610_094435_CPIR_Vent_2D_VarDens_SOS_WIP.sin"
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
if debug_mod:
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

# # Perform translation by multiplying k-space data
data_reordered = translate_factor * data_reordered

# %% Estimate sensitivity map
sens_map = np.ones((1, ImageSize, ImageSize))
# As we only have 1 channel, a sensitivity map of all 1's is assumed

# Move to GPU
mps = mvd(sens_map)

# Find radial k-space coordinates
traj_rad = np.sqrt(traj[:, :, 0]**2 + traj[:, :, 1]**2)

# Initialize
# Dimensions to deal with coronal where images are in scanner (x-z) dimension
img_full = np.zeros((ImageSize, ImageSize, N_slices))

# Reshape trajectories 
traj = np.reshape(traj, (N_samp*N_spirals, 2))

# %% Estimate density compensation functions
if use_dcf:
    # DCF settings
    beta = 8  # 8
    width = 3  # 4

    # Full image
    dcf_full = mr.pipe_menon_dcf(
        traj, img_shape=ImageShape, beta=beta, width=width, device=device)
    dcf_full = mvd(np.reshape(dcf_full, (1, np.shape(dcf_full)[0])))

# %% Loop through image slices

for imslice in range(N_slices):
    data = data_reordered[imslice, :, :]

    # Normalize k-space data to mean k0 of each key
    # Initialize variables
    data[data[:, 0] == 0] = 1e-23
    k0 = abs(data[:, 0])
    data_norm = np.zeros_like(data, dtype='complex_')

    # Calculate mean k0 values
    mean_k0 = abs(np.mean(k0[0:N_spirals]))

    for ii in range(N_spirals):
        data_norm[ii, :] = data[ii, :] * mean_k0 / k0[ii]

    # Reshape into vectors and remove the zero elements
    data_norm = np.reshape(data_norm, (N_samp*N_spirals))

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
        else:
            # Full image
            data_norm = mvd(np.reshape(data_norm, (1, np.shape(data_norm)[0])))
            img_full[:, :, imslice] = mvc(abs(mr.app.SenseRecon(
                data_norm, mps, coord=traj, lamda=lamda, max_iter=max_iter, device=device).run()))
    else:
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

# %% Visualize images
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
    montage_img = montage3d(img_full)
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
ni_img = nib.Nifti1Image(img_full, affine=aff)
nib.save(ni_img, datalocation + '/img_ventilation.nii.gz')
print("Image has been reconstructed successfully")

# %% create .exe

'''
    To create the exe file:
        	pip install pyinstaller
        	In terminal, run:
            pyinstaller --onefile --name="Philips_2DXeSpiralRecon.exe" Philips_2D_VentSpiral_recon.py
This will create a single .exe file under \dist. 

'''


