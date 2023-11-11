% Deblurring BM3D denoising demo file, based on
% Y. Mäkinen, L. Azzari, A. Foi, 2020,
% "Collaborative Filtering of Correlated Noise: Exact Transform-Domain Variance for Improved Shrinkage and Patch Matching",
% in IEEE Transactions on Image Processing, vol. 29, pp. 8339-8354.
% ---

% The location of the BM3D files -- this folder only contains demo data
addpath('bm3d');

% Experiment specifications   
imagename = 'cameraman256.png';

% Load noise-free image
y = im2double(imread(imagename));

% Generate blurry + noisy image
experiment_number = 1;

if experiment_number==1
    sigma=sqrt(2)/255;
    for x1=-7:7; for x2=-7:7; v(x1+8,x2+8)=1/(x1^2+x2^2+1); end, end; v=v./sum(v(:));
end
if experiment_number==2
    sigma=sqrt(8)/255;
    s1=0; for a1=-7:7; s1=s1+1; s2=0; for a2=-7:7; s2=s2+1; v(s1,s2)=1/(a1^2+a2^2+1); end, end;  v=v./sum(v(:));
end
if experiment_number==3
    BSNR=40; sigma=-1; % if "sigma=-1", then the value of sigma depends on the BSNR
    v=ones(9); v=v./sum(v(:));
end
if experiment_number==4
    sigma=7/255;
    v=[1 4 6 4 1]'*[1 4 6 4 1]; v=v./sum(v(:));  % PSF
end
if experiment_number==5
    sigma=2/255;
    v=fspecial('gaussian', 25, 1.6);
end
if experiment_number==6
    sigma=8/255;
    v=fspecial('gaussian', 25, .4);
end

y_blur = imfilter(y, v(end:-1:1,end:-1:1), 'circular'); % performs blurring (by circular convolution)


if sigma == -1;   %% check whether to use BSNR in order to define value of sigma
    sigma=sqrt(norm(y_blur(:)-mean(y_blur(:)),2)^2 /(size(y_blur, 1)*size(y_blur, 2)*10^(BSNR/10))); % compute sigma from the desired BSNR
end


z = y_blur + sigma*randn(size(y_blur));



% Call BM3D With the default settings.
y_est = BM3DDEB(z, sigma, v);

% To include refiltering:
%y_est = BM3DDEB(z, sigma, v, 'refilter');

% For other settings, use BM3DProfile.
% profile = BM3DProfile(); % equivalent to profile = BM3DProfile('np');
% profile.gamma = 6;  % redefine value of gamma parameter
% y_est = BM3DDEB(z, sigma, v, profile);

% Note: Although all of the examples here use white noise, you could also
% add correlated noise to the blurred image, in which case:
% y_est = BM3D(z, PSD, v, sqrt(noise_var));


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
