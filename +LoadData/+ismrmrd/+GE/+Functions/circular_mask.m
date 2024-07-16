function [mask] = circular_mask(mtx,centre,radius)
%CIRCULAR_MASK  Generate elliptical (spherical) mask
% [mask] = circular_mask(mtx,centre,radius)
% All units in matrix pixel
%        mtx  Matrix size (any dimension)
%     centre  Location of centre of sphere  (optional; default=mtx/2)
%     radius  Radius of sphere              (optional; default=mtx/2)
%       mask  Resulting mask (logical; size(mtx))
%
% 9/2016 Rolf Schulte
if nargin<1, help(mfilename); return; end

%% misc parameter checks + defaults
if ~exist('centre','var'), centre = []; end
if isempty(centre), centre = mtx/2; end
if ~exist('radius','var'), radius = []; end
if isempty(radius), radius = mtx/2; end

dim = length(mtx);
if dim<2, error('dim(=%g)<2',dim); end

centre = centre(:).';
if length(centre)~=dim,
    error('length(centre)(=%g)~=length(size(b))',length(centre),dim); 
end

radius = radius(:).';
if length(radius)==1, radius = repmat(radius,[1 3]); end
if length(radius)~=dim, 
    error('length(radius)(=%g)~=%g',length(radius),dim);
end
radius = radius./mtx;

if any(centre>mtx), 
    warning('centre>mtx');
    fprintf('centre = \n'); disp(centre);
    fprintf('mtx = \n'); disp(mtx);
end
if any(centre<1), warning('centre<1'); fprintf('centre = \n'); disp(centre); end
shft = (centre-mtx/2)./mtx;


%% generate spherical mask based on geometrical input
mask = zeros(mtx);
for l=1:dim,
    m1 = ones(1,dim);
    m1(l) = mtx(l);
    ax = linspace(-0.5,0.5,mtx(l)).'-shft(l);
    ax = reshape(ax,m1);
    
    m2 = mtx; m2(l) = 1;
    mask = mask + repmat(ax.^2/radius(l)^2,m2);
    % imagesc(mask); colorbar; drawnow; pause(1);
    % fprintf('%g\n',l); input('');
end
mask = mask<=1;


