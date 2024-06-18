function image_subplot(b,cm,title_str)
%IMAGE_SUBPLOT  Multiple 2D imagesc in subplots 
%               Dim1+2=imagesc; dim>2=subplot
% image_subplot(b,cm,title_str)
%         b  Real-valued images [x,y,...]
%        cm  colormap scale [min max]
% title_str  Title string for subplots (in cell-array)
%
% 10/2007 Rolf Schulte
if (nargin<1), help(mfilename); return; end;
if ~isreal(b(:)), warning('image_subplot:real','b must be real'); end
if ~exist('cm','var'), cm = []; end
if isempty(cm), 
    bmax = max(b(:));
    bmin = min(b(:));
    if bmin<0, 
        cm = [-1 1]*max([bmax -bmin]);
    else
        cm = [0 bmax]; 
    end
end

clf;

n = size(b);
switch length(n),
    case 0, error('bild empty');
    case 1, error('bild 1D');
    case 2, ns = 1;
    otherwise,
        ns = prod(n(3:end));
end
ncol = ceil(sqrt(ns));
b = reshape(b,n(1),n(2),ns);
for l=1:ns
    if ns>1, subplot(ceil(ns/ncol),ncol,l); end
    imagesc(b(:,:,l));
    if isnumeric(cm), 
        caxis(cm); 
    else
        tmp = b(:,:,l); tmp = max(tmp(:));
        text(2,39,sprintf('%2.2g',tmp),'color','w','FontSize',6)
    end
    axis image off; 
    pos = get(gca,'Position');
    set(gca,'Position',[pos(1:2)-0.025 1.2*pos(3:4)]);
    if exist('title_str','var'), 
        % title(title_str{l}); 
        text(2,4,title_str{l},'color','w','FontSize',8)
    end
end
drawnow;
