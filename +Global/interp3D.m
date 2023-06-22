
function [Image_3D] = interp3D(viewing_img)
        Img = squeeze(viewing_img(:,:,:,1));
        Image_size = size(Img); %size of the original image
        largest_dim=max(Image_size); %find the max dim
        [rows,cols,slices] = size(Img);
        % find steps for the meshgrid
        stepX=rows/largest_dim; 
        stepY=cols/largest_dim;
        stepZ=slices/(largest_dim+slices/2);             
        [X,Y,Z] = meshgrid(1:Image_size(2), 1:Image_size(1), 1:Image_size(3));
        [X2,Y2,Z2] = meshgrid(1:stepX:rows, 1:stepY:cols, 1:stepZ:slices);
        Image_3D = interp3(X, Y, Z, Img, X2, Y2, Z2, 'linear', 0);
end