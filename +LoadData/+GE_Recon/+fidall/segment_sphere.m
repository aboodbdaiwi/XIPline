function [mask] = segment_sphere(b,centre,radius,fac_thresh)
%SEGMENT_SPHERE  Segment image data by thresholding and building a sphere
% [mask] = segment_sphere(b,centre,radius,fac_thresh)
%          b  Image (any dimension)
%     centre  Location of centre of sphere  (units=pixel)
%     radius  Radius of sphere              (units=pixel)
% fac_thresh  Factor on threshold           (b>(mean(b(:))*fac_thresh)
%       mask  Resulting mask (logical; same dimension as b)
%
% 2/2016 Rolf Schulte
if nargin<1, help(mfilename); return; end

%% misc parameter checks + defaults
if ~exist('fac_thresh','var'), fac_thresh = []; end
if isempty(fac_thresh),        fac_thresh = 1; end

if ~isreal(b), b = abs(b); end
n = size(b);
dim = length(n);

centre = centre(:).';
if length(centre)~=dim,
    error('length(centre)(=%g)~=length(size(b))',length(centre),dim); 
end

radius = radius(:).';
if length(radius)==1, radius = repmat(radius,[1 dim]); end
if length(radius)~=dim, 
    error('length(radius)(=%g)~=%g',length(radius),dim);
end
radius = radius./n;

if any(centre>n), 
    warning('centre>n'); 
    fprintf('centre = \n'); disp(centre);
    fprintf('n = \n'); disp(n);
end
if any(centre<1), warning('centre<1'); fprintf('centre = \n'); disp(centre); end
shft = (centre-n/2)./n;

%% generate spherical mask based on geometrical input
mask = zeros(n);
for l=1:dim,
    m1 = ones(1,dim);
    m1(l) = n(l);
    ax = linspace(-0.5,0.5,n(l)).'-shft(l);
    ax = reshape(ax,m1);
    
    m2 = n; m2(l) = 1;
    mask = mask + repmat(ax.^2/radius(l)^2,m2);
    % imagesc(mask); colorbar; drawnow; pause(1);
    % fprintf('%g\n',l); input('');
end
mask = mask<=1;

%% add image thresholding mask
thresh = mean(b(:))*fac_thresh;
mask = mask & (b>thresh);

