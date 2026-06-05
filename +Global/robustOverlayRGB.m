function RGB = robustOverlayRGB(B,F,climB,overlayColor,alpha)
%ROBUSTOVERLAYRGB Create RGB overlay image without figure/axes/getframe.
%
% B = grayscale background image
% F = binary or numeric foreground mask
% climB = [low high] display range for background
% overlayColor = [R G B], values 0-1
% alpha = overlay opacity, 0-1

if nargin < 3 || isempty(climB)
    climB = [min(B(:)) max(B(:))];
end
if nargin < 4 || isempty(overlayColor)
    overlayColor = [1 0 1]; % magenta
end
if nargin < 5 || isempty(alpha)
    alpha = 1;
end

B = double(squeeze(B));
F = double(squeeze(F));

if ~isequal(size(B),size(F))
    F = imresize(F,size(B),'nearest');
end

Bgray = mat2gray(B,climB);
RGB = repmat(Bgray,[1 1 3]);

mask = F > 0;

for c = 1:3
    tmp = RGB(:,:,c);
    tmp(mask) = (1-alpha).*tmp(mask) + alpha.*overlayColor(c);
    RGB(:,:,c) = tmp;
end

RGB = im2uint8(RGB);
end