function lung_mask = SegmentLungthresh(Image,SE,nzcof)

% SE: structuring element size for erode/dilate
Image = Image(:,:,:,1);
Image = Image./max(Image(:));
[counts,~] = imhist(Image,16);
T = otsuthresh(counts);

% noise_thresh=max(Image(:));
% noise_thresh = max(abs(fft(NoiK)));
% Bzimg = Image(:,:,:,1);
% mask =zeros(size(Bzimg)); %make binary mask
% mask(Bzimg>nzcof*noise_thresh)=1;

mask = imbinarize(Image,T.*nzcof);

%errosion/dilation
[x2,y2]=meshgrid(-SE:SE,-SE:SE); 
nhood2=x2.^2+y2.^2<=SE^2;
se1=strel('arbitrary',nhood2);
mask_erode=imerode(mask,se1);
lung_mask=imdilate(mask_erode,se1); % apply mask

end