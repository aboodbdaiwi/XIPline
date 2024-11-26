function [separated_lung,lung] = separate_connected_lung(lung_binary2D)
%% Separate connected lung by # of erosion and then restore the
% original size by the same # of dilation.
% W. Zha @ 7/25/2016
[nRows, nCols] = size(lung_binary2D);
se1 = strel('disk',1);
eroded_lung = imerode(lung_binary2D,se1);
nErosion = 1;
[labeled_lung, nMid] = bwlabel(eroded_lung);
while nMid == 1 % erode by 1 pixel until separated
    eroded_lung = imerode(eroded_lung,se1);
    nErosion = nErosion + 1;
    [labeled_lung, nMid] = bwlabel(eroded_lung);
    if sum(labeled_lung(:))==0
        
        separated_lung = lung_binary2D;
        lung.right = lung_binary2D;
        lung.left = lung_binary2D;
        return;
    end
end
pixels = zeros(nMid,1);

for np = 1:nMid
    pixels(np) = sum(labeled_lung(:) == np);
end
[~, sorted_ind] = sort(pixels, 'descend');

lung.right = zeros(nRows,nCols);
lung.right(labeled_lung==sorted_ind(1)) = 1;
seN = strel('disk',nErosion);
lung.right = imdilate(lung.right,seN);
lung.left = zeros(nRows,nCols);
lung.left(labeled_lung==sorted_ind(2)) = 1;
lung.left = imdilate(lung.left,seN);
se2 = strel('disk',2);
lung.left(imdilate(lung.right,se2)==1) = 0;
separated_lung = lung.right + lung.left;
