function lung_mask = CCHMC_diff_thresh(NoiK,diffimg,SE,nzcof)
% Threshold ADC
    nzcof=nzcof/10000;
    noise_thresh=NoiK/max(NoiK);
    noise_thresh = max(abs(fft(noise_thresh)));
    Bzimg = diffimg(:,:,:,1);
    Bzimg=Bzimg-min(Bzimg(:));
    Bzimg=Bzimg/max(Bzimg(:));
    mask =zeros(size(Bzimg)); %make binary mask
    mask(Bzimg>nzcof*noise_thresh)=1;
    
    %errosion/dilation
    [x2,y2]=meshgrid(-SE:SE,-SE:SE); 
    nhood2=x2.^2+y2.^2<=SE^2;
    se1=strel('arbitrary',nhood2);
    mask_erode=imerode(mask,se1);
    lung_mask=imdilate(mask_erode,se1); % apply mask

end
