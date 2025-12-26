function denoised_images = Xe_GLHOSVD_4D(noisy_images, noise_mask, kglobal, klocal, patchsize, sw, step)
%% Function to denoise diffusion weighted HP Xe MRI images
% (supports 3D or 4D input)
%
% Inputs:
%   noisy_images: 3D (X×Y×Z) or 4D (X×Y×Z×N)
%   noise_mask  : background noise mask (same XYZ size)
%   kglobal, klocal, patchsize, sw, step: denoise params (optional)
%
% Output:
%   denoised_images: same dimensionality as input (3D or 4D)

%% If denoising parameters not specified, use these defaults
if ~exist('kglobal','var') || isempty(kglobal)
    kglobal = 0.4;
end
if ~exist('klocal','var') || isempty(klocal)
    klocal = 0.8;
end
if ~exist('patchsize','var') || isempty(patchsize)
    patchsize = 6;
end
if ~exist('sw','var') || isempty(sw)
    sw = 15;
end
if ~exist('step','var') || isempty(step)
    step = 2;
end

%% Handle 3D vs 4D
in_sz = size(noisy_images);
nd = ndims(noisy_images);

if nd == 3
    is3D = true;
    noisy4D = reshape(noisy_images, in_sz(1), in_sz(2), in_sz(3), 1); % make it 4D
elseif nd == 4
    is3D = false;
    noisy4D = noisy_images;
else
    error('Xe_GLHOSVD_4D: noisy_images must be 3D or 4D.');
end

% Basic mask sanity
if ~isequal(size(noise_mask,1), size(noisy4D,1)) || ~isequal(size(noise_mask,2), size(noisy4D,2)) || ~isequal(size(noise_mask,3), size(noisy4D,3))
    error('Xe_GLHOSVD_4D: noise_mask must match the first 3 dims of noisy_images.');
end

nVol = size(noisy4D,4);

%% Calculate base_noise_std
% If 3D, do this ONCE (no loop over volumes)
if is3D
    noise_vec_base = noisy4D(:,:,:,1);
    noise_vec_base(noise_mask~=1) = [];
    base_noise_std = std(noise_vec_base);
else
    base_noise_std = zeros(1, nVol);
    for n = 1:nVol
        noise_vec_base = noisy4D(:,:,:,n);
        noise_vec_base(noise_mask~=1) = [];
        base_noise_std(n) = std(noise_vec_base); % standard deviation of the noise
    end
end

%% Run VST (forward)
sigma = mean(base_noise_std);                  % std of noise
VST_ABC = 'B';                                 % needed in VST code
Im = Global.RiceOptVST.riceVST(noisy4D, sigma, VST_ABC); % apply VST
sigmafz = 1; %#ok<NASGU>                        % standard deviation of noise in f(z) (kept)

% Compute VST noise SD:
% If 3D, do this ONCE; if 4D, keep per-volume estimate then average
if is3D
    VST_noise_base = Im(:,:,:,1);
    VST_noise_base(noise_mask==0) = [];
    VST_std_noise_base = std(VST_noise_base);
    vstSigma = mean(VST_std_noise_base);
else
    VST_std_noise_base = zeros(1, nVol);
    for n = 1:nVol
        VST_noise_base = Im(:,:,:,n);
        VST_noise_base(noise_mask==0) = [];
        VST_std_noise_base(n) = std(VST_noise_base);
    end
    vstSigma = mean(VST_std_noise_base);
end

%% Run GLHOSVD Denoising
VST_denoised_images = Global.GL_HOSVD.glhosvd_flexible(Im, vstSigma, kglobal, klocal, patchsize, step, sw);
disp('denoised completed')

%% Inverse VST
nu_hat = Global.RiceOptVST.riceVST_EUI(VST_denoised_images, sigma, VST_ABC); % exact unbiased inverse
nu_hat = squeeze(nu_hat);

% Return same dimensionality as input
denoised_images = nu_hat;           % X×Y×Z×N


end
