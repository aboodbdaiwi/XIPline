%REDGREENBLUE
%
%REDGREENBLUE returns and rgb image with the selected colors in the input
%parameters
%
%INPUT PARAMETERS
%   blackandwhite: base image. It it the original image that will be grey
%       because it will be added to the three components (red, green and blue)
%       with the same intensity
%   red: part of the image that will be added only to the red component.
%       Previously, it will have been subtracted from the blackandwhite image
%   blue: part of the image that will be added only to the blue component.
%       Previously, it will have been subtracted from the blackandwhite image
%
%OUTPUT PARAMETERS
%   rgb: rgb image containing the input parameters to conform the total
%   image
%
%The "out-of-Matlab-Toolboxes" function(s) used by REDGREENBLUE is/are:
%   
%
%REDGREENBLUE is used by the following function(s):
%   AUTOMATICDEFECTS1slice.m

function rgb = plentyofcolors(blackandwhite,red,blue,green,yellow,pink,purple,orange)

bwnoborders=blackandwhite.*(1-and(blackandwhite,red)).*(1-and(blackandwhite,blue)).*(1-and(blackandwhite,green)).*(1-and(blackandwhite,yellow)).*(1-and(blackandwhite,pink)).*(1-and(blackandwhite,purple)).*(1-and(blackandwhite,orange));
rgb(:,:,1)=bwnoborders+red+yellow+purple+0.7*orange;
rgb(:,:,2)=bwnoborders+green+yellow+pink;
rgb(:,:,3)=bwnoborders+blue+pink+purple+0.3*orange;