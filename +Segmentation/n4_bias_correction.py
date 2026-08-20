# RH: added — ANTsPy N4 wrapper used by N4_bias_correction.m on macOS/Linux,
# where no native N4BiasFieldCorrection CLI binary is available.
"""
Usage:
    python n4_bias_correction.py <image.nii> <weight.nii> <shrink_factor> <output_corrected.nii> <output_bias.nii>
"""
import sys

import ants


def main():
    if len(sys.argv) != 6:
        print('Usage: n4_bias_correction.py <image> <weight> <shrink_factor> <output_corrected> <output_bias>')
        sys.exit(2)

    image_path, weight_path, shrink_factor, output_corrected, output_bias = sys.argv[1:6]

    image = ants.image_read(image_path)
    weight = ants.image_read(weight_path)
    shrink_factor = int(shrink_factor)

    corrected = ants.n4_bias_field_correction(image, weight_mask=weight, shrink_factor=shrink_factor)
    bias = ants.n4_bias_field_correction(image, weight_mask=weight, shrink_factor=shrink_factor,
                                          return_bias_field=True)

    ants.image_write(corrected, output_corrected)
    ants.image_write(bias, output_bias)


if __name__ == '__main__':
    main()
