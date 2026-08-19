# RH: rewritten to match Windows exactly, same reasoning as ants_register.py —
# calls the antsApplyTransforms library function directly with the same argument
# list AntsRegistration.m uses on Windows, instead of ants.apply_transforms()'s
# Python wrapper.
"""ANTsPy apply-transform wrapper used by AntsRegistration.m on macOS/Linux,
where no native antsApplyTransforms CLI binary is available.

Usage:
    python ants_apply_transform.py <fixed.nii> <moving.nii> <output.nii.gz> <transform.mat>
"""
import sys

from ants.internal import get_lib_fn, process_arguments


def main():
    if len(sys.argv) != 5:
        print('Usage: ants_apply_transform.py <fixed> <moving> <output> <transform_mat>')
        sys.exit(2)

    fixed_path, moving_path, output_path, transform_path = sys.argv[1:5]

    args = [
        '-d', '3',
        '-e', '0',
        '-i', moving_path,
        '-r', fixed_path,
        '-o', output_path,
        '-t', transform_path,
    ]
    processed = process_arguments(args)
    libfn = get_lib_fn('antsApplyTransforms')
    ret = libfn(processed)
    if ret != 0:
        raise RuntimeError('antsApplyTransforms failed with exit code %d' % ret)


if __name__ == '__main__':
    main()
