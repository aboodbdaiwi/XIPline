function BinaryMask = FreehandSegment(Image, BinaryMask)
% % SegmentLungParenchyma - Segments lung parenchyma using user-defined ROI
%   Inputs:
%       Image - 3D image data (HxWxS)
%       BinaryMask - (Optional) Initial binary mask to overlay on the image (HxWxS)
%   Outputs:
%       ImageMask - Binary mask indicating the segmented lung parenchyma
%   Package: https://github.com/aboodbdaiwi/HP129Xe_Analysis_App
%   Author: Abdullah S. Bdaiwi
%   Work email: abdullah.bdaiwi@cchmc.org
%   Personal email: abdaiwi89@gmail.com
%   Website: https://www.cincinnatichildrens.org/research/divisions/c/cpir
%   Please add updates at the end. Ex: 3/10/24 - ASB: update ....

% Check if only one argument is passed
if nargin < 2
    BinaryMask = zeros(size(Image));
end

% Normalize the image to range [0, 1]
Image = double(Image);
Image = (Image - min(Image(:))) / (max(Image(:)) - min(Image(:)));

b1img = Image(:,:,:,1);
szimg = size(b1img);
% ImageMask = zeros(szimg);

slice_num = 1;

while slice_num <= szimg(3)
    airway = 1;
    % if isempty(BinaryMask)
    %     mask = zeros(szimg(1), szimg(2));
    % else
    %     mask = BinaryMask(:,:,slice_num);
    % end
    mask = BinaryMask(:,:,slice_num);
    while airway == 1
        % Close any previous figures
        close all;
        
        % Create overlay image
        % mask = mask > 0;
        overlayedImage = imoverlay(b1img(:,:,slice_num), BinaryMask(:,:,slice_num), [1 0 0], 0.5);
        
        figure;
        imshow(overlayedImage, 'InitialMagnification', 'fit', 'DisplayRange', [0 1]);
        title(gca, ['Image Slice ', int2str(slice_num)]);
        
        % Construct a quest with options
        msg = "Select an action";
        opts = ["Add" "Remove" "Go Back Slice" "Next Slice"];
        choice = menu(msg,opts);
        % choice = questdlg('Select an action', ...
        %     'Action Choice', 'Add', 'Remove', 'Go Back Slice', 'Next Slice');

        % Handle response
        switch choice
            case 1 %'Add'
                disp('Adding region containing lung parenchyma.');
                h = drawfreehand('Color', 'r');
                BW = createMask(h);
                roi_mask = ones(szimg(1), szimg(2)) - BW;
                mask = mask + BW;
                % b1img(:,:,slice_num) = roi_mask .* double(b1img(:,:,slice_num));
                b1img(:,:,slice_num) = double(b1img(:,:,slice_num));
                delete(h);
                mask(mask > 0) = 1;
                BinaryMask(:,:,slice_num) = mask;
            case 2 %'Remove'
                disp('Removing region from lung parenchyma.');
                h = drawfreehand('Color', 'r');
                BW = createMask(h);
                mask = mask - BW;
                mask(mask < 0) = 0; % Ensure no negative values
                delete(h);
                mask(mask > 0) = 1;
                BinaryMask(:,:,slice_num) = mask;
            case 3 %'Go Back Slice' 
                disp('Going back to previous slice.');
                airway = 0;
                slice_num = max(slice_num - 1, 1);
            case 4 %'Next Slice'
                disp('Moving to next slice.');
                airway = 0;
                slice_num = slice_num + 1;                
        end
    end

    close;  

    if slice_num > szimg(3)
        disp('Segmenting lung parenchyma complete.');
        break;
    end
end
BinaryMask = BinaryMask > 0;
end

function overlayedImage = imoverlay(image, mask, color, alpha)
% Create an overlayed image with a specified color and alpha transparency.
% image - the grayscale image to overlay on
% mask - binary mask where the overlay will be applied
% color - RGB triplet specifying the overlay color
% alpha - alpha transparency value for the overlay

% Ensure the input image is in the range [0, 1]
image = double(image);
image = (image - min(image(:))) / (max(image(:)) - min(image(:)));

% Convert the grayscale image to RGB
overlayedImage = repmat(image, [1 1 3]);

% Create the color overlay
colorOverlay = bsxfun(@times, mask, reshape(color, 1, 1, 3));

% Apply the overlay with the specified alpha
overlayedImage = bsxfun(@times, overlayedImage, 1 - alpha * mask) + alpha * colorOverlay;

end
