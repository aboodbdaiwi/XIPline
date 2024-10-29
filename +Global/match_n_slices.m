
function [matched_im] = match_n_slices(Im1,Im2)
%   Inputs:
%     Im1 is the fixed one
%     Im2 is the changeable one
%   Outputs:
% 
%   Example: 
% 
%   Package: 
%
%   Author: Abdullah Bdaiwi 
%   Work email: abdullah.bdaiwi@cchmc.org
%   Personal email: abdaiwi89@gmail.com
%   Website: https://cpir.cchmc.org/

    Im1_slice_indx = 1:size(Im1,3);   
    Ratio = size(Im2, 3) ./ size(Im1, 3);

    % -- Have # of slices of Test2 equal Image3D
    scale_factor = 0.001; 
    for i = 1:500
        Img = squeeze(Im2(:,:,:,1));
        Image_size = size(Img); %size of the original image
        largest_dim = max(Image_size); %find the max dim
        [rows,cols,slices] = size(Img);
        % find steps for the meshgrid
        stepX=rows/largest_dim; 
        stepY=cols/largest_dim;
        stepZ = Ratio; % ratio between              
        [X,Y,Z] = meshgrid(1:Image_size(2), 1:Image_size(1), 1:Image_size(3));
        [X2,Y2,Z2] = meshgrid(1:stepX:rows, 1:stepY:cols, 1:stepZ:slices);
        Image_3D = interp3(X, Y, Z, Img, X2, Y2, Z2, 'linear', 0);
        if size(Image_3D,3) == size(Im1,3)
            break;
        elseif size(Image_3D,3) > size(Im1,3)
            Ratio = Ratio + scale_factor;  
        elseif size(Image_3D,3) < size(Im1,3)
            Ratio = Ratio - scale_factor;
        end
    end

    matched_im = zeros(size(Im2,1),size(Im2,2),size(Im1,3));
    matched_im(1:min(size(Im2,1),size(Image_3D,1)),1:1:min(size(Im2,2),size(Image_3D,2)),Im1_slice_indx) = Image_3D(1:min(size(Im2,1),size(Image_3D,1)),1:1:min(size(Im2,2),size(Image_3D,2)),Im1_slice_indx);         
%     figure; Global.imslice(matched_im);
end