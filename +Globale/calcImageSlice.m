% CALCIMAGESLICE  - Calculates a 2D slice of a 3D volume. 
%
% Calculates a 2D image from a ND volume at a given slice index in a given
% dimmension.
%
% slice_im = calcImageSlice(vol, sliceIdx, dim)
%
% Copyright: 2012 Scott Haile Robertson.
% Website: www.ScottHaileRobertson.com
%
% $Revision: 1.0 $  
% $Date: Nov 30, 2012 $
function slice_im = calcImageSlice(vol, sliceIdx, dim, varargin)

%Calculate the slice in the given dimmension
slice_cmd = '';
switch dim(1)
    %use permute(image_ud.srs(:,:,z,timePts(n)),[2 1 3 4])
    case 1
        slice_cmd = 'slice_im = permute(vol.Data(sliceIdx(1),:,:';
        perm_cmd = '),[2 3 1';
    case 2
        slice_cmd = 'slice_im = permute(vol.Data(:,sliceIdx(1),:';
        perm_cmd = '),[3 1 2';
    case 3
        slice_cmd = 'slice_im = permute(vol.Data(:,:,sliceIdx(1)';
        perm_cmd = '),[1 2 3';
    otherwise
        error('Invalid dimension.');
end

num_extra_dims = length(dim);
for i=2:num_extra_dims
    slice_cmd = [slice_cmd ',' num2str(sliceIdx(i))];
    perm_cmd = [perm_cmd ' ' num2str(dim(i))];
end
eval_cmd = [slice_cmd perm_cmd ']);'];

try
    eval(eval_cmd);
catch
    error(['Could not run command (calcImageSlice): ' eval_cmd]);
end
end

