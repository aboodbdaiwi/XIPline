function image_overlay(bg,fg,alpha_val,fg_thresh,mabg,mafg,havoshi)
%IMAGE_OVERLAY Overlay fg (colormap jet) to bg (colormap gray)
% image_overlay(bg,fg,alpha_val,fg_thresh,mabg,mafg,havoshi)
%        bg  Background image (colormap gray)
%        fg  Foreground image (colormap jet)
% alpha_val  Transparency: 0=bg; 1=fg                                 (0.5)
% fg_thresh  Threshold for displaying fg image (1-64)                  ([])
%      mabg  Scaling of background                             (max(bg(:)))
%      mafg  Scaling of foreground                             (max(fg(:)))
%   havoshi  Half voxel shift                                       (false)
%
%  5/2024 Rolf Schulte
if (nargin<1), help(mfilename); return; end


%% parameters
if ~exist('alpha_val','var'), alpha_val = []; end
if isempty(alpha_val), alpha_val = 0.5; end
if ~exist('fg_thresh','var'), fg_thresh = []; end
if ~exist('mabg','var'),      mabg = []; end
if ~exist('mafg','var'),      mafg = []; end
if ~exist('havoshi','var'),   havoshi = []; end
if isempty(havoshi),          havoshi = false; end


%% checks
if length(size(bg))~=2, error('bg not 2D'); end
if length(size(fg))~=2, error('fg not 2D'); end
if ~isreal(bg), warning('bg complex; taking abs'); bg = abs(bg); end
if ~isreal(fg), warning('fg complex; taking abs'); fg = abs(fg); end
if isinteger(bg), bg = double(bg); end
if isinteger(fg), fg = double(fg); end


%% half voxel shift
if havoshi
    fprintf('Shifting images by half a pixel\n');
    for ll=1:2
        kk = ((-size(fg,ll)/2):(size(fg,ll)/2-1))/size(fg,ll);
        if ll==1, kk = kk(:); end
        phafu = exp(1i*pi*kk);
        fg = ifftshift(ifft(fg,[],ll),ll);
        fg = bsxfun(@times,fg,phafu);
        fg = real(fft(fftshift(fg,ll),[],ll));
    end
    for ll=1:2
        kk = ((-size(bg,ll)/2):(size(bg,ll)/2-1))/size(bg,ll);
        if ll==1, kk = kk(:); end
        phafu = exp(1i*pi*kk);
        bg = ifftshift(ifft(bg,[],ll),ll);
        bg = bsxfun(@times,bg,phafu);
        bg = real(fft(fftshift(bg,ll),[],ll));
    end
end


%% normalise and convert to indexed
if isempty(mabg), mabg = max(abs(bg(:))); end
fprintf('mabg = %g \tmabg/max(abs(bg(:))) = %g\n',mabg,mabg/max(abs(bg(:))));
bg = 63*bg/mabg;bg = bg+1;bg = round(bg);
bg_rgb = ind2rgb(bg,gray);

if isempty(mafg), mafg = max(abs(fg(:))); end
fprintf('mafg = %g \tmafg/max(abs(fg(:))) = %g\n',mafg,mafg/max(abs(fg(:))));
fg = 63*fg/mafg;
fg = fg+1;fg = round(fg);
fg_rgb = ind2rgb(fg,jet(64));
if ~isempty(fg_thresh)
    ind = fg<fg_thresh;
    for l=1:3
        tmp1 = fg_rgb(:,:,l); 
        tmp2 = bg_rgb(:,:,l);
        tmp1(ind) = tmp2(ind);
        fg_rgb(:,:,l) = tmp1; 
    end
end


%% axes
sibg = size(bg);
sifg = size(fg);

xbg = linspace(-1,1,sibg(1)); 
ybg = linspace(-1,1,sibg(2)); 
xfg = linspace(-1,1,sifg(1)); 
yfg = linspace(-1,1,sifg(2)); 

% xbg = ((-sibg(1)/2):(sibg(1)/2-1))/sibg(1); 
% ybg = ((-sibg(2)/2):(sibg(2)/2-1))/sibg(2); 
% xfg = ((-sifg(1)/2):(sifg(1)/2-1))/sifg(1); 
% yfg = ((-sifg(2)/2):(sifg(2)/2-1))/sifg(2); 


%% actual plotting and overlaying
image(xbg,ybg,bg_rgb);
hold on
iim2 = image(xfg,yfg,fg_rgb);
set(iim2,'AlphaData',alpha_val);
hold off
axis image off; axis([-1 1 -1 1]);

end      % image_overlay.m
