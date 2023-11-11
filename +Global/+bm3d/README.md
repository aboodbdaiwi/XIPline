# MATLAB wrapper for BM3D denoising - from Tampere with love

MATLAB wrapper for BM3D for stationary correlated noise (including white noise) for color,
grayscale and multichannel images and deblurring.

BM3D is an algorithm for attenuation of additive spatially correlated
stationary (aka colored) Gaussian noise. This package provides a wrapper
for the BM3D binaries for use for grayscale, color and other multichannel images
for denoising and deblurring.

This implementation is based on Y. Mäkinen, L. Azzari, A. Foi, 2020,
"Collaborative Filtering of Correlated Noise: Exact Transform-Domain Variance for Improved Shrinkage and Patch Matching",
in IEEE Transactions on Image Processing, vol. 29, pp. 8339-8354.

The package contains the BM3D binaries compiled for:
- Windows (R2018b)
- Linux (R2017b)
- Mac OSX (R2018a)

The binaries are available for non-commercial use only. For details, see LICENSE.

Authors:
	Ymir Mäkinen   <ymir.makinen@tuni.fi>
	Lucio Azzari
	Alessandro Foi



