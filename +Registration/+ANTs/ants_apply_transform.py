# RH: new file
"""ANTsPy apply-transform wrapper used by AntsRegistration.m on macOS/Linux,
where no native antsApplyTransforms CLI binary is available.

Usage:
    python ants_apply_transform.py <fixed.nii> <moving.nii> <output.nii.gz> <transform.mat> <interpolator>
"""
import sys

import ants


def main():
    if len(sys.argv) != 6:
        print('Usage: ants_apply_transform.py <fixed> <moving> <output> <transform_mat> <interpolator>')
        sys.exit(2)

    fixed_path, moving_path, output_path, transform_path, interpolator = sys.argv[1:6]

    fixed = ants.image_read(fixed_path)
    moving = ants.image_read(moving_path)

    warped = ants.apply_transforms(fixed=fixed, moving=moving,
                                    transformlist=[transform_path],
                                    interpolator=interpolator)
    ants.image_write(warped, output_path)


if __name__ == '__main__':
    main()
