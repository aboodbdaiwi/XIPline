# RH: new file
"""ANTsPy registration wrapper used by AntsRegistration.m on macOS/Linux,
where no native antsRegistration CLI binary is available.

Usage:
    python ants_register.py <fixed.nii> <moving.nii> <output_moving.nii.gz> <output_transform.mat> <type_of_transform>
"""
import shutil
import sys

import ants


def main():
    if len(sys.argv) != 6:
        print('Usage: ants_register.py <fixed> <moving> <output_moving> <output_transform_mat> <type_of_transform>')
        sys.exit(2)

    fixed_path, moving_path, output_moving, output_transform, type_of_transform = sys.argv[1:6]

    fixed = ants.image_read(fixed_path)
    moving = ants.image_read(moving_path)

    result = ants.registration(fixed=fixed, moving=moving, type_of_transform=type_of_transform)

    ants.image_write(result['warpedmovout'], output_moving)
    shutil.copyfile(result['fwdtransforms'][0], output_transform)


if __name__ == '__main__':
    main()
