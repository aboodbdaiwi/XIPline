function oimg = miprmdebrisn(image,excSize)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function oimg = miprmdebrisn(image,excSize)
% This function removes structures (areas) smaller than "EXCSIZE" in the 
% binary image "IMAGE". The output OIMG is also a binary image.
% 
% Author:
%   Omer Demirkaya, Musa Asyali, Prasana Shaoo, ... 9/1/06
%   Medical Image Processing Toolbox
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L = bwlabeln(image,18);
S = regionprops(L,'Area');
oimg = ismember(L,find([S.Area] >= excSize));
