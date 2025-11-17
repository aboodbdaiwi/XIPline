function SNR = calcSNR(vent_img)
    % read image
    img = squeeze(single(dicomread(vent_img)));
 
    % initial segmentation to get noise
    [~, threshs] = imsegkmeans3(img(img~=0),2);
    [noise_thresh, idx] = min(threshs);
    init_noise_mask = img<noise_thresh;
    init_noise_mask(img==0) = 0;
 
    % determine noise distribution
    noise_sig = mean(img(init_noise_mask));
    noise_std = std(img(init_noise_mask));
 
    % determine threshold (rose criterion = SNR = 5)
    % SNR = mean_sig/std_noise; 5*std_noise = min_sig
    mask = img>((noise_std*5)+noise_sig);
    noise_mask = imerode(~mask, strel('cube',3));
    noise_mask(img == 0) = 0; % account for zeroed regions from grad warp
    sig_mask = imerode(mask, strel('cube',1));
 
    % calculate SNR (0.66 factor for non-gaussian noise)
    SNR = mean(img(sig_mask))/std(img(noise_mask))*0.66;
end