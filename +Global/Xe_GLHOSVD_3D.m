function denoised_images = Xe_GLHOSVD_3D(noisy_images, lung_mask, noise_mask, kglobal, klocal, patchsize, sw, step)
%% GLHOSVD denoising function for 3D Xe MRI images
%
%
% Inputs:
%
% noisy_images: original image data (3 dimensional)
% lung_mask: mask containing lung parenchyma, bright airways omitted
% noise_mask: mask of background noise with any artifacts excluded
%
%
% Outputs:
% denoised_images: 3D denoised image set
%
% Questions? email stephanie.soderlund@cchmc.org
%
%

%% If denoising parameters not specified, use these defaults

if ~exist('kglobal','var')
    kglobal = 0.4;
end

if ~exist('klocal','var')
    klocal = 0.8;
end

if ~exist('patch','var')
    patchsize = 6;
end
if ~exist('sw','var')
    sw = 15;
end

if ~exist('step','var')
    step = 2;
end

%% Calculate SNR for denoising input


    %Calculating the siganl vector
    signal_vec_base = noisy_images(:,:,:);
    signal_vec_base(lung_mask==0)=[];
    base_signal = median(signal_vec_base);
    %Calculating the noise vector
    noise_vec_base = noisy_images(:,:,:);
    noise_vec_base(noise_mask~=1)=[]; 
    base_noise_std= std(noise_vec_base);      %standard deviation of the noise
    base_SNR = base_signal / base_noise_std;  %signal to noise ratio

%% Run VST (forward)
psnr = base_SNR; 
sigma = mean(base_noise_std);                  %std of noise
VST_ABC = 'B';                                 %needed in VST code
Im = riceVST(noisy_images,sigma,VST_ABC);           %apply variance-stabilizing transformation (fz is result)
sigmafz = 1;                                   %standard deviation of noise in f(z)

    %Calculating the siganl vector
    VST_signal_base = squeeze(Im(:,:,:));
    VST_signal_base(lung_mask==0)=[];
    VST_signal = median(VST_signal_base);
    %Calculating the noise vector
    VST_noise_base = squeeze(Im(:,:,:));
    VST_noise_base(noise_mask==0)=[];               
    VST_std_noise_base = std(VST_noise_base);            %the standard deviation of the noise
    VST_SNR_base = VST_signal / VST_std_noise_base; %signal to noise ratio

%% Run GLHOSVD Denoising

tmp = (Im(:,:,:));
VST_denoised_images = glhosvd_flexible(tmp, VST_std_noise_base, kglobal, klocal, patchsize, step, sw);

disp('denoised')

%% Inverse VST
nu_hat = riceVST_EUI(VST_denoised_images,sigma,VST_ABC);   %% apply exact unbiased inverse for estimating nu
denoised_images = squeeze(nu_hat); 
end