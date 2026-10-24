# -*- coding: utf-8 -*-
"""Philips 3D hyperpolarized 129Xe multi-b-value FLORET reconstruction.

Adapted from Philips_3D_XeVent_FLORET_recon.py. Keeps its single-channel
3D NUFFT, CG, ADMM wavelet+TV, DCF and HP weighting reconstruction model.

IMPORTANT: Set B_VALUES, B_ORDER and acquisition-specific parameters below.
Raw projection-to-b-value association MUST be validated against the .list
file / pulse sequence before interpreting ADC maps.

Dependencies: the SAME environment and local modules as the ventilation script:
  sigpy_local, ReadPhilips.readphilips, functions.tv, functions.convexalg,
  functions.util, numpy, scipy, matplotlib, nibabel, scikit-learn, opencv.
"""
import os
import sys
import json
import warnings
from pathlib import Path

if getattr(sys, 'frozen', False):
    base_path = sys._MEIPASS
    os.environ['PATH'] = base_path + os.pathsep + os.environ.get('PATH', '')

import numpy as np
import matplotlib.pyplot as plt
import nibabel as nib
from scipy.io import savemat
from sklearn.cluster import KMeans
import cv2
import sigpy_local as sp
import sigpy_local.mri as mr
import ReadPhilips.readphilips as rp
from functions import tv, convexalg, util

plt.style.use('dark_background')

# ======================= USER CONFIGURATION ================================
CONFIG_PATH = r'C:\XIPline\offline_recon\config_Xe.txt'
# Override data directory here, or leave None to use config_Xe.txt.
DATA_DIRECTORY = None
FREQUENCY_OFFSET_HZ = None  # None = use config value (currently not applied; see note)

# EDIT TO MATCH THE ACQUIRED SEQUENCE. Example values from the OLD 2D script;
# these are NOT verified for the new 3D FLORET diffusion acquisition.
B_VALUES = [10, 0, 0]  # s/cm^2 if using the conventional Xe diffusion scale
# 'interleaved': b0,b1,...,bN,b0,b1,...; 'blocked': all b0, then all b1...
# 'index_file': one zero-based b-value index per raw projection in CSV/NPY.
B_ORDER = 'auto'  # 'auto' uses the explicit b-value axis for 4D raw data
B_INDEX_FILE = None
# Optional explicit projection indices per b value, overrides B_ORDER:
# e.g. {0:[0,4,8], 1:[1,5,9], ...}. Leave None for automatic grouping.
PROJECTION_INDICES = None

USE_DCF = True
SUBSAMPLE_PC = 100
NORMALIZE_DATA = False  # preserves ventilation-script default
GAMMA = 1.0            # 1 = full HP decay weighting; 0 = none
T2STAR_XE_MS = 18.0
READOUT_MS = 15.0
USE_MEASURED_K0 = True  # same weighting strategy as ventilation script

# Same defaults as ventilation reconstruction (with DCF enabled).
NUM_ITERS = 50
LAMDA_1 = 8e-5
LAMDA_2 = 1e-2
LAMDA_L2 = 1e-5
RHO = 1.0
PTOL = 1e-1
NUM_NORMAL = 4
CG_ITERS = 200
CG_PTOL = 1e-2

# As in ventilation code, 1.33x nominal scan resolution.
RECON_SCALE = 1.5
# Preserve original image manipulation: util.ir(...,1) -> util.ps -> flipud.
APPLY_VENTILATION_ORIENTATION = True
# Optional ADC map; fitting assumes correctly assigned b values and signal
# correction for RF depletion across b-values (not performed automatically).
COMPUTE_ADC = True  # separate optional post-reconstruction analysis; default OFF
MATCHED_TE_B0_GROUP = 1  # set 1 or 2 only after TE validation
ADC_RANGE = (0.0, 0.14)  # cm^2/s, only for display/summary masking
MIN_SIGNAL_FRACTION = 0.05
PLOT_RESULTS = True
PLOT_TRAJECTORIES = True
SAVE_TRAJECTORY_PLOTS = True
TRAJECTORY_MAX_PROJECTIONS = 165
TRAJECTORY_POINTS_PER_PROJECTION = 1455
TRAJECTORY_LINE_WIDTH = 0.95
TRAJECTORY_ALPHA = 0.80
TRAJECTORY_PLOT_MODE = 'line'  # 'line' or 'scatter'
TRAJECTORY_MARKER_SIZE = 2.5
# Display orientation and slice controls (does not alter reconstructed NIfTI)
DISPLAY_AXIS = 1  # 0, 1, or 2; try 1 for coronal-like view
DISPLAY_SLICE = 'auto'  # 'auto' chooses brightest plane, or set integer
DISPLAY_PERMUTE = (0, 1, 2)  # e.g. (2, 1, 0)
DISPLAY_FLIP_AXES = ()  # e.g. (0, 2); applies only to previews
DISPLAY_ROT90 = 0  # counterclockwise 90-degree rotations in preview
# Optional transform to saved volumes (leave identity unless validated)
OUTPUT_PERMUTE = (0, 1, 2)
OUTPUT_FLIP_AXES = ()
RUN_CG = False  # CG implementation retained below but disabled
SAVE_MAT = True
SAVE_PER_B = True
# ==========================================================================


def load_config():
    data_dir = DATA_DIRECTORY
    freq = FREQUENCY_OFFSET_HZ
    if data_dir is None:
        if not os.path.isfile(CONFIG_PATH):
            raise FileNotFoundError('Configuration not found: ' + CONFIG_PATH)
        with open(CONFIG_PATH, encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if not line or '=' not in line:
                    continue
                key, value = line.split('=', 1)
                if key.strip().lower() == 'datalocation':
                    data_dir = value.strip().strip('"')
                elif key.strip().lower() == 'freq_offset' and freq is None:
                    freq = float(value)
    if not data_dir:
        raise ValueError('Set DATA_DIRECTORY or datalocation in config_Xe.txt')
    root = Path(data_dir)
    if not root.is_dir():
        raise NotADirectoryError(str(root))
    raw = sorted(root.glob('*.raw'))
    if not raw:
        raw = sorted(root.glob('*.data'))
    sin = sorted(root.glob('*.sin'))
    if len(raw) != 1 or len(sin) != 1:
        raise ValueError('Expected exactly one .raw/.data and one .sin in %s; '
                         'found %d raw and %d sin. Isolate one scan per folder.'
                         % (root, len(raw), len(sin)))
    if freq not in (None, 0):
        warnings.warn('freq_offset=%s Hz was read, but no dwell-time-based '
                      'frequency demodulation is applied. Verify correction '
                      'before quantitative use.' % freq)
    return root, raw[0], sin[0], freq


def group_indices(nproj, nb, root):
    if PROJECTION_INDICES is not None:
        groups = [np.asarray(PROJECTION_INDICES[i], dtype=int) for i in range(nb)]
    elif B_ORDER == 'index_file':
        if not B_INDEX_FILE:
            raise ValueError('B_INDEX_FILE required for B_ORDER=index_file')
        p = Path(B_INDEX_FILE)
        if not p.is_absolute():
            p = root / p
        labels = np.load(p) if p.suffix.lower() == '.npy' else np.loadtxt(p, delimiter=',', dtype=int)
        labels = np.asarray(labels).reshape(-1)
        if len(labels) != nproj:
            raise ValueError('Index file has %d entries; expected %d' % (len(labels), nproj))
        groups = [np.flatnonzero(labels == i) for i in range(nb)]
    else:
        if nproj % nb:
            raise ValueError('%d projections cannot be divided evenly into %d b-values. '
                             'Supply explicit PROJECTION_INDICES or B_INDEX_FILE.' % (nproj, nb))
        if B_ORDER == 'interleaved':
            groups = [np.arange(i, nproj, nb) for i in range(nb)]
        elif B_ORDER == 'blocked':
            groups = [np.arange(i * (nproj // nb), (i + 1) * (nproj // nb)) for i in range(nb)]
        else:
            raise ValueError('Unknown B_ORDER: ' + str(B_ORDER))
    if any(len(g) == 0 for g in groups):
        raise ValueError('At least one b-value has zero projections')
    flat = np.concatenate(groups)
    if np.any(flat < 0) or np.any(flat >= nproj) or len(np.unique(flat)) != len(flat):
        raise ValueError('Projection indices out of range or duplicated')
    if len(flat) != nproj:
        warnings.warn('Only %d of %d projections assigned to b-values' % (len(flat), nproj))
    return groups


def load_philips(raw_path, sin_path):
    dl = rp.PhilipsData(str(raw_path))
    dl.raw_corr = True
    dl.compute()
    raw = np.asarray(dl.data)
    data = np.squeeze(raw)
    original_shape = data.shape
    # The first dimension of 4D diffusion data is provisionally the b-value axis.
    # C-order flattening keeps all projections for each b-value contiguous.
    if data.ndim >= 3:
        data = data.reshape(-1, data.shape[-1])
    if data.ndim != 2:
        raise ValueError('Raw data must flatten to (projections,samples); got %s' % (data.shape,))
    tr = rp.PhilipsData(str(sin_path))
    tr.compute()
    coords = np.transpose(np.squeeze(np.asarray(tr.spparams['COORDS_EXPANDED'])))
    if coords.ndim == 4:
        coords = coords.reshape(coords.shape[0]*coords.shape[1], coords.shape[2], coords.shape[3])
    if coords.ndim != 3 or coords.shape[-1] != 3:
        raise ValueError('Expected 3D FLORET coordinates (projections,samples,3); got %s' % (coords.shape,))
    # A .sin may provide a trajectory for every b-value, or just one
    # common trajectory set repeated for each b-value.
    if coords.shape[1] != data.shape[1]:
        raise ValueError('Sample mismatch: data %s vs trajectory %s' % (data.shape, coords.shape))
    if coords.shape[0] != data.shape[0]:
        if len(original_shape) == 4 and coords.shape[0] == int(np.prod(original_shape[1:-1])):
            print('Trajectory has one common set per b-value; repeating across first raw axis.')
            coords = np.tile(coords, (original_shape[0], 1, 1))
        else:
            raise ValueError('Projection mismatch: raw shape %s -> %s, trajectory %s. '
                             'Cannot safely infer b-value trajectory mapping.' %
                             (original_shape, data.shape, coords.shape))
    if len(original_shape) == 4:
        print('Detected 4D raw data:', original_shape)
        print('Assuming first axis indexes b-values; confirm with Philips LIST metadata.')
    h = tr.header.get('sin')
    def scalar(key, default):
        try:
            return float(np.asarray(h.get(key)).reshape(-1)[0])
        except (TypeError, ValueError, IndexError):
            return default
    # Some ReadPhilips builds omit valid SIN fields from the parsed header.
    # Read their literal values from the supplied ASCII .sin as fallback.
    import re
    sin_text = Path(sin_path).read_text(encoding='utf-8', errors='replace')
    def sin_values(key):
        match = re.search(r'^\s*\d+\s+\d+\s+\d+\s*:\s*' + re.escape(key) +
                          r'\s*:\s*([^\r\n]*)', sin_text, re.MULTILINE | re.IGNORECASE)
        if not match:
            return []
        return [float(v) for v in re.findall(r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[Ee][-+]?\d+)?', match.group(1))]
    def header_or_sin(key):
        value = scalar(key, np.nan)
        if np.isfinite(value):
            return value
        values = sin_values(key)
        if values:
            print('Read %s directly from ASCII SIN: %s' % (key, values))
            return values[0]
        return np.nan
    meta = {'voxel_size_header': header_or_sin('voxel_sizes'),
            'scan_resolution_header': header_or_sin('scan_resolutions'),
            'recon_resolution_header': header_or_sin('recon_resolutions'),
            'slice_thickness_header': header_or_sin('slice_thickness'),
            'diffusion_b_factors_sin': sin_values('diffusion_b_factors'),
            'raw_shape': list(original_shape)}
    return data, coords, meta


def plot_floret_trajectories(coords, groups, bvalues, image_size, out):
    """Plot ONE common FLORET trajectory (first group) with visible curves."""
    if not PLOT_TRAJECTORIES:
        return
    from matplotlib import cm
    indices = np.asarray(groups[0])
    count = max(1, int(np.floor(len(indices) * SUBSAMPLE_PC / 100)))
    indices = indices[:count]
    n_show = min(len(indices), max(1, TRAJECTORY_MAX_PROJECTIONS))
    chosen = indices[np.unique(np.linspace(0, len(indices) - 1, n_show, dtype=int))]
    fig = plt.figure(figsize=(9, 8))
    ax = fig.add_subplot(111, projection='3d')
    colors = cm.viridis(np.linspace(0.1, 1.0, len(chosen)))
    for color, idx in zip(colors, chosen):
        curve = coords[idx, :TRAJECTORY_POINTS_PER_PROJECTION, :] * image_size
        if TRAJECTORY_PLOT_MODE == 'scatter':
            ax.scatter(curve[:, 0], curve[:, 1], curve[:, 2],
                       color=color, s=TRAJECTORY_MARKER_SIZE,
                       alpha=TRAJECTORY_ALPHA, linewidths=0)
        else:
            ax.plot(curve[:, 0], curve[:, 1], curve[:, 2],
                    color=color, linewidth=TRAJECTORY_LINE_WIDTH,
                    alpha=TRAJECTORY_ALPHA)
    ax.set(xlabel='kx', ylabel='ky', zlabel='kz',
           title='3D FLORET trajectory (common to all diffusion contrasts)')
    ax.view_init(elev=24, azim=-55)
    fig.tight_layout()
    if SAVE_TRAJECTORY_PLOTS:
        fig.savefig(out / 'trajectory_FLORET_common.png', dpi=220,
                    bbox_inches='tight')
    plt.show(block=False)


def orient_volume(vol):
    """User-configurable saved-volume axis permutation/flips."""
    vol = np.transpose(vol, OUTPUT_PERMUTE)
    for ax in OUTPUT_FLIP_AXES:
        vol = np.flip(vol, axis=ax)
    return np.ascontiguousarray(vol)


def display_volume(vol):
    """Apply preview-only orientation; does not affect saved NIfTI."""
    vol = np.transpose(vol, DISPLAY_PERMUTE)
    for ax in DISPLAY_FLIP_AXES:
        vol = np.flip(vol, axis=ax)
    return vol


def select_display_slice(vol, axis=DISPLAY_AXIS):
    """Pick high-signal slice instead of assuming central slice is informative."""
    if DISPLAY_SLICE != 'auto':
        z = int(DISPLAY_SLICE)
        if not 0 <= z < vol.shape[axis]:
            raise ValueError('DISPLAY_SLICE outside volume dimension')
        return z
    # Smooth summed magnitude profile to avoid selecting isolated noisy plane.
    other_axes = tuple(j for j in range(3) if j != axis)
    profile = np.sum(np.abs(vol), axis=other_axes)
    from scipy.ndimage import gaussian_filter1d
    profile = gaussian_filter1d(profile, sigma=1.5)
    return int(np.argmax(profile))


def extract_display_slice(vol, index, axis=DISPLAY_AXIS):
    plane = np.take(vol, index, axis=axis)
    return np.rot90(plane, k=DISPLAY_ROT90)


def plot_recon_preview(volumes, bvalues, out):
    """Coronal-style slice preview with independent, robust image scaling."""
    nb = len(bvalues)
    ref = display_volume(volumes['cs'][..., 1 if nb > 1 else 0])
    z = select_display_slice(ref)
    print('Preview: axis=%d, slice=%d (of %d)' % (DISPLAY_AXIS, z, ref.shape[DISPLAY_AXIS]))
    fig, axes = plt.subplots(2, nb, figsize=(4.3*nb, 8), squeeze=False)
    for row, name in enumerate(('nufft', 'cs')):
        for i in range(nb):
            img = extract_display_slice(display_volume(volumes[name][..., i]), z)
            positive = img[np.isfinite(img) & (img > 0)]
            vmax = float(np.percentile(positive, 99.5)) if positive.size else 1.0
            axes[row, i].imshow(img, cmap='gray', vmin=0, vmax=max(vmax, 1e-12),
                                interpolation='nearest')
            axes[row, i].set_title('%s group %d (b=%g), slice %d' %
                                   (name.upper(), i, bvalues[i], z))
            axes[row, i].axis('off')
            print('  %s group %d: min %.3g max %.3g p99.5 %.3g' %
                  (name, i, np.nanmin(img), np.nanmax(img), vmax))
    fig.tight_layout()
    fig.savefig(out / 'bvalue_reconstruction_comparison.png', dpi=180)
    plt.show(block=False)


def safe_snr(image, signal_mask, noise_mask):
    a = image[signal_mask]
    b = image[noise_mask]
    if not len(a) or len(b) < 2 or np.std(b) == 0:
        return float('nan')
    return float(np.mean(a)/np.std(b))


def reconstruct_one(data_np, coords_np, image_shape, device):
    """Ventilation script's F, D, K, CG and CS-ADMM applied to ONE b-value."""
    nproj, nsamp = data_np.shape
    xp = device.xp
    # DCF per b-value, because sampling geometry may differ between b groups.
    with device:
        traj = sp.to_device(coords_np, device)
        data = sp.to_device(data_np.reshape(1, nproj, nsamp), device)
        mps = sp.to_device(np.ones((1,) + image_shape, dtype=np.int8), device)
        dcf = mr.pipe_menon_dcf(traj, img_shape=image_shape,
                                beta=8, width=4, device=device)
        dcf[:, -1] = xp.median(dcf, axis=1)
        dcf = dcf.reshape(1, nproj, nsamp)

        # Same measured k0 times assumed T2* relaxation model as vent script.
        k0 = xp.abs(data[0, :, 0])
        if NORMALIZE_DATA:
            if xp.any(k0 <= 0):
                raise ValueError('Cannot normalize projection with zero k0')
            data = data / k0[None, :, None]
        relaxation = xp.exp(-xp.arange(nsamp) * (READOUT_MS / nsamp) / T2STAR_XE_MS)
        k = k0[:, None] * relaxation[None, :] if USE_MEASURED_K0 else xp.broadcast_to(relaxation, (nproj, nsamp))
        if GAMMA != 0 and bool(sp.to_device(xp.any(k <= 0), sp.cpu_device)):
            raise ValueError('Zero k0 detected: HP weighting cannot be used safely')

        S = sp.linop.Multiply(image_shape, mps)
        F = sp.linop.NUFFT(mps.shape, coord=traj, oversamp=1.25, width=4, toeplitz=True)
        D = sp.linop.Multiply(F.oshape, dcf**0.5)
        K = sp.linop.Multiply(F.oshape, k[None, ...]**GAMMA)
        A = K * D * F * S if USE_DCF else K * F * S
        ll = sp.app.MaxEig(A.N, dtype=xp.complex64, device=device).run() * 1.01
        A = np.sqrt(1.0 / float(sp.to_device(ll, sp.cpu_device))) * A

        W1 = sp.linop.Wavelet(S.ishape, wave_name='db4')
        def g1(x):
            return LAMDA_1 * xp.linalg.norm(W1(x).ravel(), ord=1)
        prox_g1 = sp.prox.UnitaryTransform(sp.prox.L1Reg(W1.oshape, LAMDA_1), W1)
        prox_g2 = tv.ProxTV(A.ishape, LAMDA_2)
        def g2(x):
            return LAMDA_2 * xp.linalg.norm(prox_g2.G(x))

        data_norm = xp.linalg.norm(data * dcf**0.5) if USE_DCF else xp.linalg.norm(data)
        b = (data * dcf**0.5 if USE_DCF else data) / data_norm
        A_dcf = D * F * S
        b_dcf = data * dcf**0.5 / xp.linalg.norm(data * dcf**0.5)
        nufft = sp.to_device(xp.abs(A_dcf.H * b_dcf), sp.cpu_device)
        # CG intentionally disabled to focus on NUFFT and CS.
        # To restore it, uncomment the following two lines:
        # cg = convexalg.cg(num_iters=CG_ITERS, ptol=CG_PTOL, A=A, b=b,
        #                   lamda=LAMDA_L2, verbose=True, draw_output=False, save=None)
        # cg = sp.to_device(xp.abs(cg), sp.cpu_device)
        cs = convexalg.admm(num_iters=NUM_ITERS, ptol=PTOL, A=A, b=b,
                            num_normal=NUM_NORMAL, lst_proxg=[prox_g1, prox_g2],
                            rho=RHO, lst_g=[g1, g2], method='cg',
                            verbose=True, draw_output=False)
        cs = sp.to_device(xp.abs(cs), sp.cpu_device)
    nufft = np.asarray(nufft, dtype=np.float32)
    cs = np.asarray(cs, dtype=np.float32)
    if APPLY_VENTILATION_ORIENTATION:
        nufft, cs = [np.flipud(util.ps(util.ir(x, 1), 0, 0)) for x in (nufft, cs)]
    return orient_volume(nufft), orient_volume(cs)


def adc_fit(images, bvalues, mask):
    """Log-linear monoexponential ADC (cm^2/s for b in s/cm^2).

    NOTE: no RF flip-angle correction is applied. If different b values
    have unequal magnetization depletion, ADC is biased. This output is
    provisional until acquisition-specific correction is validated.
    """
    b = np.asarray(bvalues, dtype=float)
    if len(np.unique(b)) < 2:
        raise ValueError('ADC requires >=2 distinct b-values')
    ref = images[..., np.argmin(b)]
    good = mask & (ref > 0)
    good &= np.all(np.isfinite(images), axis=-1)
    good &= np.all(images > MIN_SIGNAL_FRACTION * ref[..., None], axis=-1)
    log_s = np.log(np.maximum(images, np.finfo(np.float32).tiny))
    bc = b - b.mean()
    slope = np.tensordot(log_s, bc, axes=([-1], [0])) / np.sum(bc**2)
    adc = (-slope).astype(np.float32)
    adc[~good] = np.nan
    return adc


def main():
    bvalues = np.asarray(B_VALUES, dtype=float)
    if bvalues.ndim != 1 or len(bvalues) < 2:
        raise ValueError('Set B_VALUES to the actual 2+ acquired diffusion b-values')
    root, raw_path, sin_path, freq = load_config()
    out = root / 'results_diffusion_floret'
    out.mkdir(exist_ok=True)
    print('Raw:', raw_path, '\nTrajectory:', sin_path)
    data, coords, meta = load_philips(raw_path, sin_path)
    sin_b = meta.get('diffusion_b_factors_sin', [])
    if sin_b:
        print('Diffusion b-factors recorded in SIN:', sin_b)
        if len(sin_b) == len(bvalues) and not np.allclose(sin_b, bvalues):
            warnings.warn('Configured B_VALUES %s differ from SIN diffusion_b_factors %s. '
                          'Verify ordering before ADC fitting.' % (bvalues.tolist(), sin_b))
    nproj, nsamp = data.shape
    if B_ORDER == 'auto':
        raw_shape = meta['raw_shape']
        if len(raw_shape) != 4:
            raise ValueError('B_ORDER=auto requires 4D raw data. Set interleaved/blocked/index_file explicitly.')
        if len(bvalues) != raw_shape[0]:
            raise ValueError('Raw first axis contains %d groups, but B_VALUES has %d entries (%s). '
                             'Set B_VALUES to the actual acquired values in raw-axis order.' %
                             (raw_shape[0], len(bvalues), bvalues.tolist()))
        n_per_b = int(np.prod(raw_shape[1:-1]))
        groups = [np.arange(i*n_per_b, (i+1)*n_per_b) for i in range(len(bvalues))]
    else:
        groups = group_indices(nproj, len(bvalues), root)
    print('Raw projections:', nproj, 'samples:', nsamp)
    for i, g in enumerate(groups):
        print('b=%g: %d projections, first indices %s' % (bvalues[i], len(g), g[:12]))
    if not np.isfinite(meta['scan_resolution_header']):
        raise ValueError('Could not read scan_resolutions from .sin')
    scan_resolution = int(np.ceil(RECON_SCALE * meta['scan_resolution_header']))
    image_shape = (scan_resolution,) * 3
    print('3D reconstruction shape:', image_shape)
    plot_floret_trajectories(coords, groups, bvalues, scan_resolution, out)
    try:
        device = sp.Device(0)
        with device:
            _ = sp.to_device(np.ones(1), device)
        print('Using GPU 0')
    except Exception as exc:
        print('GPU unavailable; using CPU:', exc)
        device = sp.Device(-1)

    # Same orientation matrix as ventilation; intentionally no invented
    # scanner-space voxel scaling or position (source vent affine is unit).
    aff = np.diag([-1., 1., -1., 1.])
    nb = len(bvalues)
    volumes = {name: np.zeros(image_shape + (nb,), dtype=np.float32)
               for name in ('nufft', 'cs')}
    stats = []
    for i, indices in enumerate(groups):
        count = max(1, int(np.floor(len(indices) * SUBSAMPLE_PC / 100)))
        indices = indices[:count]
        print('\n=== Reconstructing b=%g (%d projections) ===' % (bvalues[i], len(indices)))
        # Important: each b-value gets its own matching k-space coordinates.
        tr = coords[indices].astype(np.float32) * scan_resolution
        d = data[indices].astype(np.complex64)
        nufft, cs = reconstruct_one(d, tr, image_shape, device)
        for name, vol in zip(('nufft', 'cs'), (nufft, cs)):
            if vol.shape != image_shape:
                raise ValueError('%s reconstructed shape %s != %s' % (name, vol.shape, image_shape))
            volumes[name][..., i] = vol
            if SAVE_PER_B:
                nib.save(nib.Nifti1Image(vol, aff), str(out / ('img_%s_group%02d_b%g.nii.gz' % (name, i, bvalues[i]))))
        # Same 2-class K-means masking approach as ventilation script.
        km = KMeans(n_clusters=2, random_state=0, n_init=10).fit(cs.reshape(-1, 1))
        means = km.cluster_centers_.ravel()
        signal = (km.labels_.reshape(image_shape) == np.argmax(means))
        noise = ~signal
        # The source's cv2.erode is 2D; apply slice by slice explicitly.
        noise = np.stack([cv2.erode(noise[..., z].astype(np.uint8), np.ones((7, 7), np.uint8), iterations=5)
                          for z in range(image_shape[2])], axis=-1).astype(bool)
        row = {'b_value': float(bvalues[i]), 'n_projections': int(len(indices))}
        for name in volumes:
            row['snr_' + name] = safe_snr(volumes[name][..., i], signal, noise)
        stats.append(row)
        if i == 1:  # first b=0 reference group; verify against TE metadata
            reference_mask, reference_noise = signal, noise
        print('SNR:', row)

    for name, vol in volumes.items():
        nib.save(nib.Nifti1Image(vol, aff), str(out / ('img_diffusion_%s_4d.nii.gz' % name)))
        if SAVE_MAT:
            savemat(str(out / ('img_diffusion_%s_4d.mat' % name)),
                    {'img_diffusion_' + name: vol, 'b_values': bvalues}, do_compression=True)
    nib.save(nib.Nifti1Image(reference_mask.astype(np.uint8), aff), str(out / 'img_diffusion_mask.nii.gz'))
    nib.save(nib.Nifti1Image(reference_noise.astype(np.uint8), aff), str(out / 'img_diffusion_noise_mask.nii.gz'))

    # ---------------- OPTIONAL ADC POST-PROCESSING ----------------
    # b=10,0,0 includes two distinct TE contrasts. ADC requires the b=0
    # acquisition with TE matching b=10. Set MATCHED_TE_B0_GROUP after
    # checking Philips metadata. Disabled by default.
    adc_statistics = {}
    if COMPUTE_ADC:
        if MATCHED_TE_B0_GROUP is None:
            raise ValueError('Set MATCHED_TE_B0_GROUP (1 or 2) before ADC calculation')
        hi = int(np.argmax(bvalues))
        lo = int(MATCHED_TE_B0_GROUP)
        if hi == lo or bvalues[hi] <= bvalues[lo]:
            raise ValueError('Invalid matched-TE ADC group assignment')
        for name in ('nufft', 'cs'):
            high = volumes[name][..., hi]
            low = volumes[name][..., lo]
            valid = reference_mask & (high > 0) & (low > 0)
            adc = np.full(image_shape, np.nan, np.float32)
            adc[valid] = (-np.log(high[valid] / low[valid]) /
                          (bvalues[hi] - bvalues[lo]))
            nib.save(nib.Nifti1Image(adc, aff), str(out / ('ADC_%s_UNCORRECTED.nii.gz' % name)))
            if SAVE_MAT:
                savemat(str(out / ('ADC_%s_UNCORRECTED.mat' % name)), {'ADC': adc})
            # Descriptive statistics inside the reference lung mask only.
            # No physiologic clipping: inspect outliers before interpretation.
            adc_values = adc[valid & np.isfinite(adc)]
            count = int(adc_values.size)
            mean_adc = float(np.mean(adc_values)) if count else None
            sd_adc = float(np.std(adc_values, ddof=1)) if count > 1 else None
            adc_statistics[name] = {'mean_adc_cm2_s': mean_adc,
                                    'sd_adc_cm2_s': sd_adc, 'n_voxels': count}
            if count > 1:
                print(f'{name.upper()} ADC: {mean_adc:.5f} ± {sd_adc:.5f} cm²/s '
                      f'(mean ± sample SD; n={count} voxels; uncorrected)')
            else:
                print(f'{name.upper()} ADC: insufficient valid voxels (n={count})')
    # -------------------------------------------------------------

    metadata = dict(meta, b_values=bvalues.tolist(), b_order=B_ORDER,
                    projection_counts=[len(g) for g in groups], frequency_offset_hz=freq,
                    reconstruction_shape=list(image_shape), use_dcf=USE_DCF,
                    gamma=GAMMA, t2star_ms=T2STAR_XE_MS, readout_ms=READOUT_MS,
                    lamda_1=LAMDA_1, lamda_2=LAMDA_2, rho=RHO,
                    num_iters=NUM_ITERS, snr=stats, adc_statistics=adc_statistics,
                    display_axis=DISPLAY_AXIS, display_permute=list(DISPLAY_PERMUTE),
                    output_permute=list(OUTPUT_PERMUTE), run_cg=False,
                    warnings=['Verify raw projection ordering from acquisition metadata',
                              'ADC has no diffusion-specific RF depletion correction',
                              'NIfTI affine copied from ventilation script, not spatially calibrated',
                              'freq_offset not applied'])
    with open(out / 'reconstruction_summary.json', 'w', encoding='utf-8') as f:
        json.dump(metadata, f, indent=2)

    if PLOT_RESULTS:
        plot_recon_preview(volumes, bvalues, out)
    if COMPUTE_ADC:
        adc = nib.load(str(out / 'ADC_cs_UNCORRECTED.nii.gz')).get_fdata()
        ref = display_volume(volumes['cs'][..., 1])
        z = select_display_slice(ref)
        fig2, ax = plt.subplots(figsize=(6, 5))
        im = ax.imshow(extract_display_slice(display_volume(adc), z),
                       cmap='jet', vmin=ADC_RANGE[0], vmax=ADC_RANGE[1])
        ax.set_title('Provisional ADC (uncorrected)')
        ax.axis('off')
        fig2.colorbar(im, ax=ax, label='ADC (cm²/s)')
        fig2.tight_layout()
        fig2.savefig(out / 'ADC_cs_UNCORRECTED_preview.png', dpi=180)
    if PLOT_RESULTS or PLOT_TRAJECTORIES or COMPUTE_ADC:
        plt.show()
    print('\nFinished. Results saved in:', out)


if __name__ == '__main__':
    main()
