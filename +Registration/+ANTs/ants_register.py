# RH: rewritten to match Windows exactly (previous version used ants.registration()'s
# Python convenience wrapper, whose internal defaults silently differ from the schedule
# XIPline uses on Windows: gradient step is hardcoded to 0.25 there vs Affine[0.1] here,
# sampling rate defaults to 0.2 vs 1 (100%) here, iterations/shrink-factors/smoothing-sigmas
# all differ, and there is no way to request BSpline interpolation through that wrapper).
"""ANTsPy registration wrapper used by AntsRegistration.m on macOS/Linux, where no
native antsRegistration CLI binary is available.

This calls the same compiled antsRegistration library function that both the
native Windows .exe and ANTsPy itself wrap, with the exact argument list
AntsRegistration.m uses on Windows, so registration behaves identically on
both platforms (same algorithm, schedule, sampling, and interpolation).

Usage:
    python ants_register.py <fixed.nii> <moving.nii> <output_moving.nii.gz> <output_prefix> <transform_spec>

    transform_spec is the ANTs --transform argument, e.g. "Affine[0.1]"
"""
import sys

from ants.internal import get_lib_fn, process_arguments


def main():
    if len(sys.argv) != 6:
        print('Usage: ants_register.py <fixed> <moving> <output_moving> <output_prefix> <transform_spec>')
        sys.exit(2)

    fixed_path, moving_path, output_moving, output_prefix, transform_spec = sys.argv[1:6]

    args = [
        '--dimensionality', '3',
        '--float', '0',
        '--interpolation', 'BSpline',
        '--metric', 'MI[%s,%s,1,32,Regular,1]' % (fixed_path, moving_path),
        '--transform', transform_spec,
        '--convergence', '[20x20x20,1e-6,20]',
        '--shrink-factors', '4x2x1',
        '--smoothing-sigmas', '0x0x0',
        '--output', '[%s,%s]' % (output_prefix, output_moving),
        '--verbose', '1',
    ]
    processed = process_arguments(args)
    libfn = get_lib_fn('antsRegistration')
    ret = libfn(processed)
    if ret != 0:
        raise RuntimeError('antsRegistration failed with exit code %d' % ret)


if __name__ == '__main__':
    main()
