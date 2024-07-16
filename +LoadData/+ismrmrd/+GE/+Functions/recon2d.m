function b=recon2d(d,res,fftshiftdir)
%b=recon2d(d,res,fftshiftdir)
% Reconstruct raw MRI data with 2d IFFT
%
%           d  data (as returned from read_MR_rawdata) 
%              [ky, kx, echo, slice, rcvr]
%         res  anticipated resolution (optinal)
% fftshiftdir  1=1;2=2;3=all;rest=none
%           b  reconstructed image
%
% 10/2008 Rolf Schulte
if (nargin<1), help(mfilename); return; end;
si = size(d);
if ~exist('res','var'), res=[]; end
if ~exist('fftshiftdir','var'), fftshiftdir = 0; end
if isempty(res),
    zf = false;
    res = [si(1) si(2)]; 
else
    zf = true;
    switch length(res)
        case 1, res = [res res];
        case 2,
        otherwise, error('length(res) ~=1 or 2');
    end
    if any(res<[si(1) si(2)]), warning('interpolating down!???'); end
end

switch length(si),
    case 2, n_coils = 1;
    case 5, n_coils = si(5);
        d = squeeze(d);
        if si(3)~=1, si,error('only 1 echo implemented'); end
        if si(4)~=1, si,error('only 1 slice implemented'); end
    case 6, n_coils = si(6);
        d = squeeze(d);
        if any(si(3:5)~=1), si,error('any(si(3:5))~=1'); end
    otherwise, error('length(si)~=2 or 5 (6)');
end

% reconstruction of images
%b = zeros(si(2),si(1),n_coils);
b = zeros(res(2),res(1),n_coils);
for l=1:n_coils,
    % transpose -> 1.dim: x-axis; 2.dim: y-axis
    if zf,
        tmp = zero_fill(d(:,:,l).',[res(2) res(1)]);
    else
        tmp = d(:,:,l).';
    end
    tmp = ifft2(ifftshift(tmp));
    switch fftshiftdir
        case 1, b(:,:,l) = fftshift(tmp,1);
        case 2, b(:,:,l) = fftshift(tmp,2);
        case 3, b(:,:,l) = fftshift(tmp);
        otherwise, b(:,:,l) = tmp;
    end
end

