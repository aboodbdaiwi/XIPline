# RH: new file — used by VentilationFunctions/N4_bias_correction.m on macOS/Linux,
# where no native N4BiasFieldCorrection CLI binary is available. Calls ANTs' own
# compiled N4BiasFieldCorrection library function directly (same approach as
# Registration.ANTs's ants_register.py) with the exact same argument list
# N4_bias_correction.m uses on Windows, instead of ants.n4_bias_field_correction()'s
# Python wrapper, which does not expose mask + weight + convergence + spline +
# histogram settings all together the way the CLI does.
"""
Usage:
    python n4_bias_correction.py <image> <mask> <weight> <shrink_factor> <convergence> <bspline> <histogram> <output_corrected> <output_bias>

    convergence, bspline, and histogram are the ANTs -c/-b/-t argument values,
    e.g. "[25,0]", "[1x1x14,3]", "[0.75,0.01,100]"
"""
import sys

from ants.internal import get_lib_fn, process_arguments


def main():
    if len(sys.argv) != 10:
        print('Usage: n4_bias_correction.py <image> <mask> <weight> <shrink_factor> '
              '<convergence> <bspline> <histogram> <output_corrected> <output_bias>')
        sys.exit(2)

    image, mask, weight, shrink, convergence, bspline, histogram, output_corrected, output_bias = sys.argv[1:10]

    args = [
        '-d', '3',
        '-i', image,
        '-s', shrink,
        '-x', mask,
        '-w', weight,
        '-c', convergence,
        '-b', bspline,
        '-t', histogram,
        '-o', '[%s,%s]' % (output_corrected, output_bias),
        '-v', '1',  # RH: verbose, so MATLAB's system() call actually shows N4 progress
    ]
    processed = process_arguments(args)
    libfn = get_lib_fn('N4BiasFieldCorrection')
    ret = libfn(processed)
    if ret != 0:
        raise RuntimeError('N4BiasFieldCorrection failed with exit code %d' % ret)


if __name__ == '__main__':
    main()
