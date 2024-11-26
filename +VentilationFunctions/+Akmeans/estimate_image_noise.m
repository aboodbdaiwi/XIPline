function [sampled_noise,noise_mean,noise_std] = estimate_image_noise(gas,lungsmask)
%
% W. Zha @ August 2016
lungbox  =  zeros(size(gas));
plotting_range = get_plotting_range(lungsmask,10);
        lungbox(plotting_range.y,plotting_range.x,plotting_range.z)=1;
        bkgbox  =   zeros(size(gas));
        noise_range  =  get_plotting_range(lungsmask,20);
        
        bkgbox(noise_range.y,noise_range.x,noise_range.z)=1;
        bkgmask  =  bkgbox-lungbox;
        sampled_noise  =  gas(bkgmask(:)>0);
        noise_mean = mean(sampled_noise);
        noise_std = std(sampled_noise);