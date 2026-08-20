% RH: new file — same purpose as Global.robustOverlayRGB (build the composited
% RGB image via pure array math, no figure/axes/getframe involved) but for a
% multi-color colormap overlay (e.g. VentilationFunctions.calculate_LB_VDP's
% 6-bin defect map), mirroring what Global.imoverlay does with a cmap+climF
% argument. Needed because getframe() on an invisible/minimized axes returned
% a cropped/zoomed-in capture on macOS (see calculate_LB_VDP.m's BinnedVentmap
% and calculate_VDP_CCHMC.m's VentDefectmap for the reports this broke).
function RGB = robustColormapOverlayRGB(B,F,climB,climF,cmap,alpha)
%ROBUSTCOLORMAPOVERLAYRGB Overlay a colormap-indexed value array onto a
%grayscale image using pure array math (no figure/axes/getframe).
%
% B = grayscale background image
% F = foreground value array (colormap-indexed, e.g. discrete bin labels)
% climB = [low high] display range for background
% climF = [low high] value range mapped across cmap
% cmap = Nx3 colormap (rows = colors, columns = R,G,B in 0-1)
% alpha = overlay opacity where F >= climF(1), 0-1 (default 1)

if nargin < 6 || isempty(alpha)
    alpha = 1;
end

B = double(squeeze(B));
F = double(squeeze(F));

if ~isequal(size(B),size(F))
    F = imresize(F,size(B),'nearest');
end

Bgray = mat2gray(B,climB);
RGB = repmat(Bgray,[1 1 3]);

nColors = size(cmap,1);
% Map F values onto colormap rows the same way imagesc/caxis does.
idx = round((F - climF(1)) / (climF(2) - climF(1)) * (nColors-1)) + 1;
idx = min(max(idx,1),nColors);

mask = F >= climF(1);

for c = 1:3
    channel = cmap(:,c);
    overlayChannel = reshape(channel(idx(:)), size(F));
    tmp = RGB(:,:,c);
    tmp(mask) = (1-alpha).*tmp(mask) + alpha.*overlayChannel(mask);
    RGB(:,:,c) = tmp;
end

RGB = im2uint8(RGB);
end
