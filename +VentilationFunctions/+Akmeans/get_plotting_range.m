function [plotting_range,crop_mask]= get_plotting_range(lungmask,offset)
%
% Copyright Wei Zha @ 2014


if nargin<2
    offset=10;
end
[nRows,nCols,nSlices]=size(lungmask);

lungmask_sum=sum(lungmask,3);
lungmask_sum(lungmask_sum>0)=1;

nonzero_pixels = find(lungmask_sum==1);

[lung_x,lung_y]=meshgrid(1:nCols,1:nRows);
row_indices = lung_x(nonzero_pixels);
col_indices = lung_y(nonzero_pixels);

plotting_range.x = max(min(row_indices)-offset,1):min(max(row_indices)+offset,nCols);
plotting_range.y = max(min(col_indices)-offset,1):min(max(col_indices)+offset,nRows);
plotting_range.z = find(squeeze(sum(sum(lungmask,1),2)));
% plotting_range.z = max(plotting_range.z(1)-offset,1): min(plotting_range.z(end),nSlices);
crop_mask = zeros(nRows,nCols,nSlices);
crop_mask(plotting_range.y,plotting_range.x,plotting_range.z) = 1;
