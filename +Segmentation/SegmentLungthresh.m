function lung_mask = SegmentLungthresh(Image,SE,nzcof)
%   Inputs:
%      
%   Outputs:
%                   
%   Package: https://github.com/aboodbdaiwi/HP129Xe_Analysis_App
%
%   Author: Abdullah S. Bdaiwi
%   Work email: abdullah.bdaiwi@cchmc.org
%   Personal email: abdaiwi89@gmail.com
%   Website: https://www.cincinnatichildrens.org/research/divisions/c/cpir
%
%   Please add updates at the end. Ex: 3/10/24 - ASB: update .... 

% SE: structuring element size for erode/dilate
Image = Image(:,:,:,1);
Image = Image./max(Image(:));
[counts,~] = imhist(Image,256);
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