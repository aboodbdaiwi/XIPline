% SHOWIMAGESLICE  - Displays a 2D slice of a 3D volume. 
%
% Displays a 2D image from a 3D volume at a given slice index in a given
% dimmension.
%
% hImage = showImageSlice(vol, slice, dim)
%
% hImage = showImageSlice(vol, slice, dim, hImage)
%
% Copyright: 2012 Scott Haile Robertson.
% Website: www.ScottHaileRobertson.com
%
% $Revision: 1.0 $  
% $Date: Nov 28, 2012 $
function hImage = showImageSlice(vol, sliceIdx, dim, varargin)

%Calculate the slice in the given dimmension
slice_im = Globale.calcImageSlice(vol, sliceIdx, dim);

%Show the slice in the given image if its given
if(nargin > 3)
    hImage = varargin{1};
    set(hImage,'CData',slice_im);
    axis image;
else
    %No image given, so just make a new one
    hImage = imagesc(slice_im);
end