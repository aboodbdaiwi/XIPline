function image_overlay(bg,fg,alpha_val,fg_thresh,mabg,mafg)
%image_overlay(bg,fg,alpha_val,fg_thresh,mabg,mafg)
%  Overlay fg (colormap jet) to bg (colormap gray)
%
%        bg  Background image (colormap gray)
%        fg  Foreground image (colormap jet)
% alpha_val  Transparency (0=bg; 1=fg; optional; default=0.5)
% fg_thresh  Threshold for displaying fg image (1-64)
%      mabg  Scaling of background (optional)
%      mafg  Scaling of foreground (optional)
%
% 10/2007 Rolf Schulte
if (nargin<1), help(mfilename); return; end;


%% parameters
if ~exist('alpha_val','var'), alpha_val = []; end
if isempty(alpha_val), alpha_val = 0.5; end
if ~exist('fg_thresh','var'), fg_thresh = []; end
if ~exist('mabg','var'), mabg = []; end
if ~exist('mafg','var'), mafg = []; end
if length(size(bg))~=2, error('bg not 2D'); end
if length(size(fg))~=2, error('fg not 2D'); end
%if any(size(bg)~=size(fg)), warning('sizes bg/fg mismatch'); end
if ~isreal(bg), error('bg complex'); end
if ~isreal(fg), error('fg complex'); end
if isinteger(bg), bg = double(bg); end
if isinteger(fg), fg = double(fg); end


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
if ~isempty(fg_thresh),
    ind = fg<fg_thresh;
    for l=1:3, 
        tmp1 = fg_rgb(:,:,l); 
        tmp2 = bg_rgb(:,:,l);
        tmp1(ind) = tmp2(ind);
        fg_rgb(:,:,l) = tmp1; 
    end
end

% axes
% sibg = size(bg);
% sifg = size(fg);
% sima = max([sibg ; sifg]);
% xbg = 1:(sima(1)/sibg(1)):sima(1);
% ybg = 1:(sima(2)/sibg(2)):sima(2);
% xfg = 1:(sima(1)/sifg(1)):sima(1);
% yfg = 1:(sima(2)/sifg(2)):sima(2);

xbg = linspace(-1,1,size(bg,1)); 
ybg = linspace(-1,1,size(bg,2)); 
xfg = linspace(-1,1,size(fg,1)); 
yfg = linspace(-1,1,size(fg,2)); 


%% actual plotting and overlaying
image(xbg,ybg,bg_rgb);
hold on
iim2 = image(xfg,yfg,fg_rgb);
set(iim2,'AlphaData',alpha_val);
hold off
axis image off; axis([-1 1 -1 1]);
