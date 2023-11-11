% Multichannel (for color images, see bm3d_demo_rgb.m) BM3D denoising demo file, based on 
% Y. Mäkinen, L. Azzari, A. Foi, 2020,
% "Collaborative Filtering of Correlated Noise: Exact Transform-Domain Variance for Improved Shrinkage and Patch Matching",
% in IEEE Transactions on Image Processing, vol. 29, pp. 8339-8354.

% For multichannel images, block matching is performed on the first channel.
% ---
% The location of the BM3D files -- this folder only contains demo data
addpath('bm3d');

% Experiment specifications   

% Load noise-free image

% The multichannel example data is acquired from: http://www.bic.mni.mcgill.ca/brainweb/
% (C.A. Cocosco, V. Kollokian, R.K.-S. Kwan, A.C. Evans,
%  "BrainWeb: Online Interface to a 3D MRI Simulated Brain Database"
% NeuroImage, vol.5, no.4, part 2/4, S425, 1997
% -- Proceedings of 3rd International Conference on Functional Mapping of the Human Brain, Copenhagen, May 1997.

data_name = 'brainslice.mat';
y = load(data_name, 'slice_sample');
y = y.slice_sample;

% Possible noise types to be generated 'gw', 'g1', 'g2', 'g3', 'g4', 'g1w',
% 'g2w', 'g3w', 'g4w'.
noise_type =  'g2';
noise_var = 0.02; % Noise variance
seed = 0; % seed for pseudorandom noise realization
% Generate noise with given PSD
[noise, PSD, kernel] = getExperimentNoise(noise_type, noise_var, seed, size(y));
% N.B.: For the sake of simulating a more realistic acquisition scenario,
% the generated noise is *not* circulant. Therefore there is a slight
% discrepancy between PSD and the actual PSD computed from infinitely many
% realizations of this noise with different seeds.

% Generate noisy image corrupted by additive spatially correlated noise
% with noise power spectrum PSD
z = y + noise;

% Call BM3D With the default settings.
y_est = BM3D(z, PSD);

% If different channels have different PSDs, they can be passed
% concatenated in the 3rd dimension (for all settings):
% y_est = BM3D(z, cat(3, PSD1, PSD2));

% To include refiltering:
% y_est = BM3D(z, PSD, 'refilter');

% For other settings, use BM3DProfile.
% profile = BM3DProfile(); % equivalent to profile = BM3DProfile('np');
% profile.gamma = 6;  % redefine value of gamma parameter
% y_est = BM3D(z, PSD, profile);

% Note: For white noise, you may instead of the PSD 
% also pass a standard deviation 
% y_est = BM3D(z, sqrt(noise_var));
% or multiple standard deviations
% y_est = BM3D(z, [sigma1, sigma2]);


psnr = getPSNR(y, y_est)

% PSNR ignoring 16-pixel wide borders (as used in the paper), due to refiltering potentially leaving artifacts
% on the pixels near the boundary of the image when noise is not circulant
psnr_cropped = getCroppedPSNR(y, y_est, [16, 16])

figure,
subplot(1, 3, 1);
imshow(y);
title('y');
subplot(1, 3, 2);
imshow(z);
title('z');
subplot(1, 3, 3);
imshow(y_est);
title('y_{est}');

