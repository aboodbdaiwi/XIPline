function [CT_resized,ct_mask_resized] = resample_CT_to_proton(CTdata, CT_header,proton_header,ct_mask)
% resize CT to proton resolution
% W. Zha@4/1/2016
%%
if nargin<4
    ct_mask = zeros(size(CTdata));
end




proton_x  = double(proton_header.Width);
proton_y = double(proton_header.Height);
ct_x_mm     = double(round(CT_header.Width*CT_header.PixelSpacing(1)));
ct_y_mm    = double(round(CT_header.Width*CT_header.PixelSpacing(2)));


proton_resX  = double(proton_header.PixelSpacing(1));
proton_resY = double(proton_header.PixelSpacing(2));

% Calculate the desired number of voxels for each dim
ct_x_pixel_new = round(ct_x_mm / proton_resX);
ct_y_pixel_new = round(ct_y_mm / proton_resY);

if mod(ct_x_pixel_new,2)==1
    ct_x_pixel_new = ct_x_pixel_new - 1;
end

if mod(ct_y_pixel_new,2)==1
    ct_y_pixel_new = ct_y_pixel_new - 1;
end

nSlices=size(CTdata,3);
CT_resized = zeros(proton_x,proton_y,nSlices);
ct_mask_resized = zeros(proton_x,proton_y,nSlices);
for nsl = 1:nSlices
    
    ct_slice = CTdata(:, :, nsl);
    sl_resized = imresize(ct_slice,[ct_x_pixel_new,ct_y_pixel_new]);
    rowpad = round((proton_x-size(sl_resized,1))/2);
    colpad = round((proton_y-size(sl_resized,2))/2);
    CT_resized(:,:,nsl) = padarray(sl_resized,[rowpad, colpad],-2000,'both');
    
    mask_resized = imresize(ct_mask(:,:,nsl),[ct_x_pixel_new,ct_y_pixel_new],'nearest');
    ct_mask_resized(:,:,nsl) = padarray(mask_resized,[rowpad, colpad],0,'both');

end
