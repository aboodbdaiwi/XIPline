function se = strel_sphere(r,dim)
%STREL_SPHERE  Create spherical structuring element
%  se = strel_sphere(r,dim)
%    r  Radius 
%  dim  Spatial dimensions (2D or 3D)
%   se  Spehrical structuring element
%
%See also IMDILATE.
% 7/2016 Rolf Schulte
if nargin<1, help(mfilename); return; end

%% 3D structuring element: sphere with radius r
switch dim
    case 1, [x1] = ndgrid(-r:r); x3 = zeros(size(x1)); x2 = x3;
    case 2, [x1,x2] = ndgrid(-r:r); x3 = zeros(size(x1));
    case 3, [x1,x2,x3] = ndgrid(-r:r);
    otherwise, error('dim(=%g)~=1 or 2 or 3',dim);
end

se = strel(sqrt(x1.^2 + x2.^2 + x3.^2) <=r);
