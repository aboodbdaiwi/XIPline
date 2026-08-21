function dd=zte_image_fid(d,what,ind1,ind4,ind6,cm)
%ZTE_IMAGE_FID  Plot FID of ZTE data via imagesc_row
%zte_image_fid(d,what,ind1,ind4,ind6,cm)
%                                                          (default)
%    d  raw data
% what  what to plot                                       (1)
%       1=log(abs(d)); 2=abs(d); 3=real(d); 4=imag(d); 5=angle(d);
% ind1  index list for first dimension  (spokes)
% ind4  index list for fourth dimension (echoes)
% ind6  index list for sixth dimension  (coils)            (1)
%   cm  range of caxis
%
% 11/2018 Rolf Schulte
if (nargin<1), help(mfilename); return; end

if ~exist('what','var'), what = []; end
if isempty(what),        what = 1; end
if ~exist('ind1','var'), ind1 = []; end
if ~exist('ind4','var'), ind4 = []; end
if ~exist('ind6','var'), ind6 = []; end
if isempty(ind6),        ind6 = 1; end
if length(ind6)~=1, warning('length(ind6)(=%g)~=1',length(ind6)); end
if ~exist('cm','var'),   cm = []; end


%% reshape + parameters
n = size(d);
if ((length(n)~=5) && (length(n)~=6))
    warning('length(size(d))(=%g)~=5 or 6',length(n));
end

dd = zeros(n(1)*n(5),n(2),1,n(4));
for l4=1:n(4)
    for l5=1:n(5)
        dd((1:n(1))+(l5-1)*n(1),:,1,l4)=d(:,:,1,l4,l5,ind6);
    end
end
if isempty(ind1), ind1 = 1:(n(1)*n(5)); end
if isempty(ind4), ind4 = 1:n(4); end
dd = dd(ind1,:,:,ind4);

switch what
    case 1, dd = log(abs(dd)); t_str = 'log(abs(dd))';
    case 2, dd = abs(dd);      t_str = 'abs(dd)';
    case 3, dd = real(dd);     t_str = 'real(dd)';
    case 4, dd = imag(dd);     t_str = 'imag(dd)';
    case 5, dd = angle(dd);    t_str = 'angle(dd)';
    otherwise, error('what (=%g) unknown',what);
end

if any(isnan(dd(:))), warning('dd isnan; zeroing'); dd(isnan(dd)) = 0; end
if any(isinf(dd(:))), fprintf('dd isinf; zeroing\n'); dd(isinf(dd)) = 0; end
if isempty(cm)
    max_val = max(abs(dd(:)));
    if any(dd(:)<0)
        cm = max_val*[-1 1];
    else
        cm = [0 max_val];
    end
    if what==5, cm = [-pi pi]; end
end

%% actual plotting
clf
imagesc_row(dd,cm);
title(t_str); colorbar
axis normal

end   % main function zte_image_fid.m
