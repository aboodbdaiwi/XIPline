function MASK = SegmentXeLung(V,MASK)
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

% if ~any(MASK(:))
    % No initial mask, generate one
    nzcof = 1; 
    SE = 1;
    Image = V(:,:,:,1);
    Image = Image./max(Image(:));
    [counts,~] = imhist(Image,16);
    T = otsuthresh(counts);
    mask = imbinarize(Image,T.*nzcof);

    %errosion/dilation
    [x2,y2]=meshgrid(-SE:SE,-SE:SE); 
    nhood2=x2.^2+y2.^2<=SE^2;
    se1=strel('arbitrary',nhood2);
    mask_erode=imerode(mask,se1);
    
% end

MASK = imdilate(mask_erode,se1); % apply mask    

%--------------------------------------------------------------------------
