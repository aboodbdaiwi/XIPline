function lung_mask = ADC_thresh2(diffimg,SE,nzcof)
% Threshold ADC

noise_thresh=max(diffimg(:));
% noise_thresh = max(abs(fft(NoiK)));
Bzimg = diffimg(:,:,:,1);
mask =zeros(size(Bzimg)); %make binary mask
mask(Bzimg>nzcof*noise_thresh)=1;

%errosion/dilation
[x2,y2]=meshgrid(-SE:SE,-SE:SE); 
nhood2=x2.^2+y2.^2<=SE^2;
se1=strel('arbitrary',nhood2);
mask_erode=imerode(mask,se1);
lung_mask=imdilate(mask_erode,se1); % apply mask

end
