# -*- coding: utf-8 -*-
"""Quantitative Philips 3D hyperpolarized 129Xe multi-b-value FLORET reconstruction.

Adapted from Philips_3D_XeVent_FLORET_recon.py. Keeps its single-channel 3D NUFFT and wavelet+TV ADMM reconstruction,
while applying SHARED scaling across all acquisitions for quantitative ADC.
Uses an explicitly normalized forward-adjoint power iteration to avoid unstable
MaxEig(A_raw.N) evaluation; never normalizes b-value volumes independently.
ADDITIONAL VALIDATION REQUIRED: acquisition index/order, complex raw data gains,
DCF consistency, and RF/HP-decay correction.

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
NORMALIZE_DATA = False  # MUST remain False for quantitative inter-contrast fitting
GAMMA = 1.0            # 1 = full HP decay weighting; 0 = none
T2STAR_XE_MS = 18.0
READOUT_MS = 15.0
USE_MEASURED_K0 = True  # k0 is RELATIVE within each group; NOT global amplitude

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
APPLY_VENTILATION_ORIENTATION = False  # avoid unverified util.ps / util.ir intensity changes
# Optional ADC map; fitting assumes correctly assigned b values and signal
# correction for RF depletion across b-values (not performed automatically).
COMPUTE_ADC = True
MATCHED_TE_B0_GROUP = 1  # set 1 or 2 only after TE validation
ADC_RANGE = (0.0, 0.14)  # cm^2/s, only for display/summary masking
MIN_SIGNAL_FRACTION = 0.05
# Do not force the observed ADC to match any expected physiological range.
# This is a QC-only range for secondary descriptive statistics.
ADC_QC_RANGE = (0.0, 0.14)
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
    # CAUTION: raw first axis is ASSUMED to index diffusion contrasts for B_ORDER='auto'.
    # Check Philips LIST/sequence metadata! Do not trust dimensionality alone.
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


def _to_float(value):
    """Move a scalar from CPU/GPU to a plain Python float."""
    return float(np.asarray(sp.to_device(value, sp.cpu_device)).reshape(-1)[0])


def _stable_device_norm(value, xp):
    """Compute norm in float64 precision to reduce float32 overflow risk."""
    return _to_float(xp.linalg.norm(value.astype(xp.complex128, copy=False)))


def _shared_operator_scale(A_raw, device, n_iter=20):
    """Estimate ONE safe operator scale from the matched-TE b=0 reference.

    Instead of MaxEig(A_raw.N), which may use an unstable optimized Toeplitz
    normal operator, explicitly apply A_raw.H(A_raw(x)) with normalized probes.
    A preliminary forward-probe normalization reduces overflow risk.

    The returned scalar MUST be reused for all other diffusion contrasts.
    This is an operator-norm preconditioner, not a per-b-value image rescaling.
    """
    xp = device.xp
    rng = np.random.default_rng(20261110)
    shape = tuple(A_raw.ishape)
    probe_cpu = (rng.standard_normal(shape).astype(np.float32) +
                 1j*rng.standard_normal(shape).astype(np.float32))
    probe_cpu /= np.linalg.norm(probe_cpu.astype(np.complex128))

    with device:
        x = sp.to_device(probe_cpu.astype(np.complex64), device)
        # First check forward operator without computing its normal operator.
        try:
            probe_response = A_raw * x
            probe_gain = _stable_device_norm(probe_response, xp)
        except Exception as exc:
            raise RuntimeError('Failed forward-operator probe; inspect NUFFT coordinates, '
                               'DCF and k0 weighting') from exc
        if not np.isfinite(probe_gain) or probe_gain <= 0:
            raise ValueError('Non-finite or zero forward-operator norm (%.6g). '
                             'Check DCF, trajectory units and k0 weighting.' % probe_gain)

        pre_scale = float(1.0 / probe_gain)
        A_pre = pre_scale * A_raw
        for it in range(n_iter):
            # Explicit composition, avoiding the NUFFT/Toeplitz .N shortcut.
            y = A_pre.H * (A_pre * x)
            ynorm = _stable_device_norm(y, xp)
            if not np.isfinite(ynorm) or ynorm <= 0:
                raise ValueError('Operator power iteration %d returned %.6g. '
                                 'Cannot safely run CS with this acquisition.' % (it+1, ynorm))
            x = y / ynorm

        # For unit-norm x, ||A_pre x||^2 is the positive Rayleigh quotient.
        lam = _stable_device_norm(A_pre * x, xp) ** 2
        if not np.isfinite(lam) or lam <= 0:
            raise ValueError('Power-iteration eigenvalue invalid (%.6g)' % lam)
        # Like the original code, add 1 percent margin.
        scale = float(pre_scale / np.sqrt(1.01 * lam))
        if not np.isfinite(scale) or scale <= 0:
            raise ValueError('Final shared operator scale invalid (%.6g)' % scale)

    print('  [Operator power iteration] forward probe gain %.6g | '
          'preconditioned eigenvalue %.6g | SHARED operator scale %.6g' %
          (probe_gain, lam, scale))
    return scale


def reconstruct_one(data_np, coords_np, image_shape, device, *,
                    common_data_norm=None, common_A_scale=None, group_label=''):
    """Reconstruct one contrast using SHARED data and operator scaling.

    The first call (matched-TE b=0) establishes BOTH normalization constants.
    ALL other contrasts receive those identical constants. The measured k0
    weighting is normalized *within each group* so that k0 amplitude is NOT
    accidentally treated as an extrinsic gain and divided out of image contrast.
    """
    if NORMALIZE_DATA:
        raise ValueError('NORMALIZE_DATA must be False for quantitative ADC/T2*')
    nproj, nsamp = data_np.shape
    xp = device.xp
    with device:
        traj = sp.to_device(coords_np, device)
        data = sp.to_device(data_np.reshape(1, nproj, nsamp), device)
        mps = sp.to_device(np.ones((1,) + image_shape, dtype=np.int8), device)
        dcf = mr.pipe_menon_dcf(traj, img_shape=image_shape,
                                beta=8, width=4, device=device)
        dcf[:, -1] = xp.median(dcf, axis=1)
        dcf = dcf.reshape(1, nproj, nsamp)
        # An invalid or highly ill-scaled DCF may cause a zero/NaN MaxEig.
        # Do NOT silently replace it or normalize each group independently.
        dcf_ok = bool(sp.to_device(xp.all(xp.isfinite(dcf) & (dcf > 0)), sp.cpu_device))
        if not dcf_ok:
            raise ValueError('DCF has zero, negative or non-finite entries. '
                             'Check trajectory units and pipe_menon_dcf output.')
        dcf_min = _to_float(xp.min(dcf))
        dcf_med = _to_float(xp.median(dcf))
        dcf_max = _to_float(xp.max(dcf))
        print('  [DCF QC] min %.6g | median %.6g | max %.6g' %
              (dcf_min, dcf_med, dcf_max))

        # k0 modulation tracks *relative temporal decay* within this contrast,
        # not its absolute image amplitude. Using raw k0 as K would absorb b-value
        # attenuation into the forward operator, making high-b and low-b images
        # artificially similar, even with shared b normalization.
        k0 = xp.abs(data[0, :, 0])
        positive = k0 > 0
        if not bool(sp.to_device(xp.all(positive), sp.cpu_device)) and GAMMA != 0 and USE_MEASURED_K0:
            raise ValueError('Zero k0 sample: cannot safely create relative HP weights')
        median_k0 = _to_float(xp.median(k0[positive]))
        if median_k0 <= 0 or not np.isfinite(median_k0):
            raise ValueError('Invalid k0 median')
        relaxation = xp.exp(-xp.arange(nsamp, dtype=xp.float32)
                            * (READOUT_MS / nsamp) / T2STAR_XE_MS)
        if USE_MEASURED_K0:
            k = (k0 / median_k0)[:, None] * relaxation[None, :]
        else:
            k = xp.broadcast_to(relaxation, (nproj, nsamp))

        S = sp.linop.Multiply(image_shape, mps)
        F = sp.linop.NUFFT(mps.shape, coord=traj,
                           oversamp=1.25, width=4, toeplitz=True)
        D = sp.linop.Multiply(F.oshape, dcf**0.5)
        K = sp.linop.Multiply(F.oshape, k[None, ...]**GAMMA)
        A_raw = K * D * F * S if USE_DCF else K * F * S

        # Estimate the operator scale once from the matched-TE b=0 group.
        # A_raw.N may invoke a custom Toeplitz normal approximation that
        # becomes zero/NaN; explicit forward+adjoint power iteration is safer.
        # Crucially, all other groups use the SAME common_A_scale.
        if common_A_scale is None:
            common_A_scale = _shared_operator_scale(A_raw, device)
        A = common_A_scale * A_raw

        W1 = sp.linop.Wavelet(S.ishape, wave_name='db4')
        def g1(x):
            return LAMDA_1 * xp.linalg.norm(W1(x).ravel(), ord=1)
        prox_g1 = sp.prox.UnitaryTransform(sp.prox.L1Reg(W1.oshape, LAMDA_1), W1)
        prox_g2 = tv.ProxTV(A.ishape, LAMDA_2)
        def g2(x):
            return LAMDA_2 * xp.linalg.norm(prox_g2.G(x))

        weighted_data = data * dcf**0.5 if USE_DCF else data
        if common_data_norm is None:
            common_data_norm = _stable_device_norm(weighted_data, xp)
            if common_data_norm <= 0 or not np.isfinite(common_data_norm):
                raise ValueError('Invalid shared data norm')
        b = weighted_data / common_data_norm
        # NUFFT adjoint: the same data normalization across all b/TE groups.
        # Use the same DCF forward/adjoint convention for every group.
        A_dcf = D * F * S
        b_dcf = data * dcf**0.5 / common_data_norm
        nufft = sp.to_device(xp.abs(A_dcf.H * b_dcf), sp.cpu_device)

        # CS image uses the shared operator scale, data normalization and
        # regularization weights. Unlike NUFFT, CS is nonlinear and should
        # be validated against an amplitude-known phantom / synthetic data.
        cs = convexalg.admm(num_iters=NUM_ITERS, ptol=PTOL, A=A, b=b,
                            num_normal=NUM_NORMAL, lst_proxg=[prox_g1, prox_g2],
                            rho=RHO, lst_g=[g1, g2], method='cg',
                            verbose=True, draw_output=False)
        cs = sp.to_device(xp.abs(cs), sp.cpu_device)

        dc = sp.to_device(dcf, sp.cpu_device)
        report = dict(group=group_label, median_abs_raw_k0=median_k0,
                      median_dcf=float(np.median(dc)),
                      common_data_norm=common_data_norm,
                      common_A_scale=common_A_scale)
        print('  [Quant scales] group=%s raw k0 median %.6g | DCF median %.6g | '
              'common_data_norm %.6g | common_A_scale %.6g' %
              (group_label, median_k0, report['median_dcf'],
               common_data_norm, common_A_scale))

    nufft = np.asarray(nufft, dtype=np.float32)
    cs = np.asarray(cs, dtype=np.float32)
    if APPLY_VENTILATION_ORIENTATION:
        warnings.warn('APPLY_VENTILATION_ORIENTATION=True calls util.ir/util.ps; '
                      'confirm they preserve intensity before fitting ADC')
        nufft, cs = [np.flipud(util.ps(util.ir(x, 1), 0, 0))
                     for x in (nufft, cs)]
    return (orient_volume(nufft), orient_volume(cs),
            common_data_norm, common_A_scale, report)


def adc_fit_polyfit(image_b0, image_high, b_high, mask):
    """Voxelwise log-linear fitting with np.polyfit, matching MATLAB polyfit.

    np.polyfit can fit MANY voxels at once with y=(2,nVoxels): each column
    receives its own independent linear fit. Equal-TE images are required.
    Returns raw signed ADC (NaN outside valid signal voxels) and fit mask.
    """
    ref = np.asarray(image_b0, dtype=np.float64)
    high = np.asarray(image_high, dtype=np.float64)
    valid = (np.asarray(mask, dtype=bool) &
             np.isfinite(ref) & np.isfinite(high) &
             (ref > 0) & (high > 0) &
             (high >= MIN_SIGNAL_FRACTION * ref))
    adc = np.full(ref.shape, np.nan, dtype=np.float32)
    if np.any(valid):
        log_signals = np.vstack((np.log(ref[valid]), np.log(high[valid])))
        # Each voxel is independently fit: log(S) = intercept - b*ADC.
        slopes = np.polyfit(np.array([0., float(b_high)]), log_signals, 1)[0]
        adc[valid] = (-slopes).astype(np.float32)
    return adc, valid


def print_adc_qc(name, low, high, adc, valid):
    """Diagnostics: raw amplitude ratio, raw ADC, and QC-only range summary."""
    n = int(np.count_nonzero(valid))
    if n == 0:
        print('%s: no valid ADC voxels' % name)
        return {'n': 0}
    v = adc[valid].astype(float)
    ratios = (high[valid] / low[valid]).astype(float)
    within = np.isfinite(v) & (v > ADC_QC_RANGE[0]) & (v < ADC_QC_RANGE[1])
    report = dict(n=n, mean_adc_raw=float(np.mean(v)),
                  median_adc_raw=float(np.median(v)),
                  sd_adc_raw=float(np.std(v, ddof=1)) if n > 1 else None,
                  median_high_over_low=float(np.median(ratios)),
                  fraction_adc_negative=float(np.mean(v < 0)),
                  fraction_in_qc_range=float(np.mean(within)))
    if np.any(within):
        report['mean_adc_qc'] = float(np.mean(v[within]))
        report['sd_adc_qc'] = (float(np.std(v[within], ddof=1))
                               if int(np.count_nonzero(within)) > 1 else None)
    print('  %s: median Shigh/Sb0=%.4f | ADC raw mean=%.5f median=%.5f | '
          'negative=%.1f%% | ADC QC 0..0.14 fraction=%.1f%% | QC mean=%s' %
          (name, report['median_high_over_low'], report['mean_adc_raw'],
           report['median_adc_raw'], 100*report['fraction_adc_negative'],
           100*report['fraction_in_qc_range'],
           ('%.5f' % report['mean_adc_qc']) if 'mean_adc_qc' in report else 'N/A'))
    return report


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
    # Raw k0 amplitudes are untouched by reconstruction. Diagnose before CS.
    print('\n=== RAW k-SPACE CENTER CONTRAST DIAGNOSTICS ===')
    raw_k0_medians = []
    for i, group in enumerate(groups):
        med = float(np.median(np.abs(data[group, 0])))
        raw_k0_medians.append(med)
        print('  raw group %d, b=%g: median abs(first sample)=%.6g' %
              (i, bvalues[i], med))
    if len(groups) >= 2:
        ratio = raw_k0_medians[0] / max(raw_k0_medians[1], 1e-20)
        print('  RAW k0 b10/b0-long ratio: %.4f (if first sample is k0)' % ratio)
        print('  For ADC=0.055 and delta-b=10, ideal equal-gain ratio: %.4f' %
              np.exp(-0.055*10))
        if len(groups[0]) == len(groups[1]):
            cr = np.max(np.abs(coords[groups[0]]-coords[groups[1]]))
            print('  Max absolute trajectory difference b10 vs b0 long: %.6g' % cr)
    if len(groups) >= 3 and len(groups[1]) == len(groups[2]):
        cr = np.max(np.abs(coords[groups[1]]-coords[groups[2]]))
        print('  Max absolute trajectory difference b0 long vs b0 short: %.6g' % cr)
    # Keep raw group-level signal contrasts: NO independent mean/max/RMS normalization.

    try:
        device = sp.Device(0)
        with device:
            _ = sp.to_device(np.ones(1), device)
        print('Using GPU 0')
    except Exception as exc:
        print('GPU unavailable; using CPU:', exc)
        device = sp.Device(-1)

    # Validate b group selection before reconstruction. Do not silently mix TE.
    if MATCHED_TE_B0_GROUP != 1 or list(bvalues) != [10.0, 0.0, 0.0]:
        warnings.warn('This script was audited for [b10-long, b0-long, b0-short] '
                      'in groups [0,1,2]. Verify metadata manually!')
    # Same orientation matrix as ventilation; intentionally no invented
    # scanner-space voxel scaling or position (source vent affine is unit).
    aff = np.diag([-1., 1., -1., 1.])
    nb = len(bvalues)
    volumes = {name: np.zeros(image_shape + (nb,), dtype=np.float32)
               for name in ('nufft', 'cs')}
    stats = [None] * nb
    reconstruction_scales = [None] * nb
    common_data_norm = None
    common_A_scale = None
    # Reconstruct matched-TE b=0 FIRST to define normalization/operator scales.
    # This does NOT alter group ordering in the saved 4D image.
    reconstruction_order = [int(MATCHED_TE_B0_GROUP)] + [
        i for i in range(nb) if i != int(MATCHED_TE_B0_GROUP)]
    print('Reconstruction processing order (b0 reference first):', reconstruction_order)
        
    # ============================================================
    # RAW b=0 LONG-TE vs SHORT-TE SIGNAL DIAGNOSTIC
    # ============================================================

    idx_long = groups[1]
    idx_short = groups[2]

    assert len(idx_long) == len(idx_short), \
        "Long and short TE groups have different projection counts"

    # Assumes the first readout sample is at k-space center
    S_long = np.abs(data[idx_long, 0]).astype(np.float64)
    S_short = np.abs(data[idx_short, 0]).astype(np.float64)

    valid = (
        np.isfinite(S_long) &
        np.isfinite(S_short) &
        (S_long > 0) &
        (S_short > 0)
    )

    ratio = S_short[valid] / S_long[valid]

    print("\n===== RAW T2* SIGNAL DIAGNOSTICS =====")
    print("Median long TE k0:", np.median(S_long[valid]))
    print("Median short TE k0:", np.median(S_short[valid]))
    print("Median short/long ratio:", np.median(ratio))
    print("Ratio > 1:", 100*np.mean(ratio > 1), "%")

    print("Expected ratio T2*=18 ms:", np.exp(10/18))
    print("Expected ratio T2*=24 ms:", np.exp(10/24))

    # Verify that the trajectories are paired
    traj_difference = np.max(
        np.abs(coords[idx_long] - coords[idx_short])
    )
    print("Maximum trajectory difference:", traj_difference)

    # Check that the first samples are near k-space center
    print("Median initial |k|, long TE:",
        np.median(np.linalg.norm(coords[idx_long, 0, :], axis=1)))
    print("Median initial |k|, short TE:",
        np.median(np.linalg.norm(coords[idx_short, 0, :], axis=1)))

    # Plot paired raw k0 amplitudes
    plt.figure(figsize=(10, 4))

    plt.subplot(1, 2, 1)
    plt.plot(S_long, label="b=0 Long TE")
    plt.plot(S_short, label="b=0 Short TE")
    plt.xlabel("Projection")
    plt.ylabel("Raw k0 magnitude")
    plt.title("Raw signal comparison")
    plt.legend()
    plt.grid(True)

    plt.subplot(1, 2, 2)
    plt.hist(ratio, bins=50)
    plt.axvline(np.exp(10/18), color='r',
                linestyle='--', label='Expected 18 ms')
    plt.axvline(np.exp(10/24), color='b',
                linestyle='--', label='Expected 24 ms')
    plt.xlabel("Short TE / Long TE")
    plt.ylabel("Projection count")
    plt.title("Raw signal ratio")
    plt.legend()
    plt.grid(True)

    plt.tight_layout()
    plt.show()

    # Three acquisitions:
    # 0 = b10 long TE
    # 1 = b0 long TE
    # 2 = b0 short TE (nominal)

    S1 = np.abs(data[groups[0], 0]).astype(float)
    S2 = np.abs(data[groups[1], 0]).astype(float)
    S3 = np.abs(data[groups[2], 0]).astype(float)

    valid = (
        np.isfinite(S1) & np.isfinite(S2) & np.isfinite(S3) &
        (S1 > 0) & (S2 > 0) & (S3 > 0)
    )

    S1, S2, S3 = S1[valid], S2[valid], S3[valid]

    r12 = S1 / S2
    r32 = S3 / S2
    r31 = S3 / S1

    print('\n=== THREE-CONTRAST RAW SIGNAL ANALYSIS ===')
    print(f'S1/S2 (b10 long / b0 long): {np.median(r12):.4f}')
    print(f'S3/S2 (short / b0 long):  {np.median(r32):.4f}')
    print(f'S3/S1 (short / b10 long):{np.median(r31):.4f}')

    print('\nExpected S3/S1 if both have b10 and')
    print('their TE difference is 10 ms:', np.exp(10/18))

    # Hypothetical T2* if S1 and S3 have matched b-values
    if np.median(r31) > 1:
        apparent_t2 = 10 / np.log(np.median(r31))
        print(f'Hypothetical S1/S3 T2*: {apparent_t2:.2f} ms')

    plt.figure(figsize=(10, 5))
    plt.plot(r12, label='S1/S2: b10 / b0 long')
    plt.plot(r32, label='S3/S2: short / long')
    plt.plot(r31, label='S3/S1: short / high b long')

    plt.axhline(np.exp(10/18), color='r', linestyle='--',
                label='Expected 18 ms TE ratio')
    plt.xlabel('Projection index')
    plt.ylabel('Raw k0 magnitude ratio')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()
    
    for i in reconstruction_order:
        indices = groups[i]
        count = max(1, int(np.floor(len(indices) * SUBSAMPLE_PC / 100)))
        indices = indices[:count]
        print('\n=== Reconstructing group=%d b=%g (%d projections) ===' %
              (i, bvalues[i], len(indices)))
        # Important: each b-value uses its own matching k-space coordinates.
        tr = coords[indices].astype(np.float32) * scan_resolution
        d = data[indices].astype(np.complex64)
        nufft, cs, common_data_norm, common_A_scale, scale_report = reconstruct_one(
            d, tr, image_shape, device,
            common_data_norm=common_data_norm, common_A_scale=common_A_scale,
            group_label='group%02d_b%g' % (i, bvalues[i]))
        reconstruction_scales[i] = scale_report
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
        stats[i] = row
        if i == int(MATCHED_TE_B0_GROUP):  # matched long-TE b=0 reference
            reference_mask, reference_noise = signal, noise
        print('SNR:', row)

    for name, vol in volumes.items():
        nib.save(nib.Nifti1Image(vol, aff), str(out / ('img_diffusion_%s_4d.nii.gz' % name)))
        if SAVE_MAT:
            savemat(str(out / ('img_diffusion_%s_4d.mat' % name)),
                    {'img_diffusion_' + name: vol, 'b_values': bvalues}, do_compression=True)
    nib.save(nib.Nifti1Image(reference_mask.astype(np.uint8), aff), str(out / 'img_diffusion_mask.nii.gz'))
    nib.save(nib.Nifti1Image(reference_noise.astype(np.uint8), aff), str(out / 'img_diffusion_noise_mask.nii.gz'))

    # ---------------- ADC: voxelwise log-linear np.polyfit ----------------
    # Pair exactly the SAME-TE contrasts: group0=b10 long TE, group1=b0 long TE.
    # Group2=b0 short TE is excluded from ADC fitting.
    adc_statistics = {}
    if COMPUTE_ADC:
        hi = 0
        lo = int(MATCHED_TE_B0_GROUP)
        if hi == lo or bvalues[hi] <= bvalues[lo]:
            raise ValueError('Invalid matched-TE ADC group assignment')
        for name in ('nufft', 'cs'):
            high = volumes[name][..., hi]
            low = volumes[name][..., lo]
            adc, valid = adc_fit_polyfit(low, high, bvalues[hi]-bvalues[lo],
                                         reference_mask)
            nib.save(nib.Nifti1Image(adc, aff), str(out / ('ADC_%s_UNCORRECTED.nii.gz' % name)))
            if SAVE_MAT:
                savemat(str(out / ('ADC_%s_UNCORRECTED.mat' % name)),
                        {'ADC': adc, 'valid_mask': valid.astype(np.uint8),
                         'b_values': np.array([bvalues[lo], bvalues[hi]])})
            adc_statistics[name] = print_adc_qc(name.upper(), low, high, adc, valid)
    # ----------------------------------------------------------------------

    metadata = dict(meta, b_values=bvalues.tolist(), b_order=B_ORDER,
                    projection_counts=[len(g) for g in groups], frequency_offset_hz=freq,
                    reconstruction_shape=list(image_shape), use_dcf=USE_DCF,
                    gamma=GAMMA, t2star_ms=T2STAR_XE_MS, readout_ms=READOUT_MS,
                    lamda_1=LAMDA_1, lamda_2=LAMDA_2, rho=RHO,
                    num_iters=NUM_ITERS, snr=stats, adc_statistics=adc_statistics,
                    shared_data_norm=common_data_norm, shared_operator_scale=common_A_scale,
                    group_scaling_diagnostics=reconstruction_scales,
                    display_axis=DISPLAY_AXIS, display_permute=list(DISPLAY_PERMUTE),
                    output_permute=list(OUTPUT_PERMUTE), run_cg=False,
                    warnings=['Verify raw projection ordering from acquisition metadata',
                              'ADC still requires independent validation: raw gains, RF/HP, DCF, CS nonlinearity',
                              'Unverified util.ps and util.ir intentionally disabled',
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
