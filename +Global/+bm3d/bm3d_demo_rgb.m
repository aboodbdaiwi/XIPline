% RGB image BM3D denoising demo file, based on 
% Y. Mäkinen, L. Azzari, A. Foi, 2020,
% "Collaborative Filtering of Correlated Noise: Exact Transform-Domain Variance for Improved Shrinkage and Patch Matching",
% in IEEE Transactions on Image Processing, vol. 29, pp. 8339-8354.

% For color images, block matching is performed on the luminance channel.
% ---
% The location of the BM3D files -- this folder only contains demo data
addpath('bm3d');

% Experiment specifications   
imagename = 'image_Lena512rgb.png';

% Load noise-free image
y = im2double(imread(imagename));

% Possible noise types to be generated 'gw', 'g1', 'g2', 'g3', 'g4', 'g1w',
% 'g2w', 'g3w', 'g4w'.
noise_type =  'g3';
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

% Call (color) BM3D With the default settings.
y_est = CBM3D(z, PSD);

% To include refiltering:
%y_est = CBM3D(z, PSD, 'refilter');

% For other settings, use BM3DProfile.
% profile = BM3DProfile(); % equivalent to profile = BM3DProfile('np');
% profile.gamma = 6;  % redefine value of gamma parameter
% y_est = CBM3D(z, PSD, profile);

% Note: For white noise, you may instead of the PSD 
% also pass a standard deviation 
% y_est = CBM3D(z, sqrt(noise_var));


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
