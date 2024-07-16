function mm = gridding(d,k,dcf,res)
%GRIDDING Interface to grid_interp_mex.c
% m = gridding(d,k,dcf,res)
%
% INPUT:
%     d  Original k-space data sampled along an abitrary trajectory
%        [#reps,#k-space points]
%     k  Arbitrary k-space trajectory specified as k=kx+i*ky;
%   dcf  Local area surrounding the samples along the abitrary 
%        k-space trajectory -> e.g. k_area=voronoi_area(k_positions)                   
%   res  desired Cartesian resoltion
% OUTPUT:
%     m  Reconstructed image including Deapodization
%
% core gridding stuff by Philip Beatty
%
% 11/2020 Rolf Schulte
% See also GRID_INTERP, DEAPODIZE
if (nargin<1), help(mfilename); return; end

% Initialize Parameters
oversamp = 1.5;
% nreps = size(d,1);
switch length(res)
    case 1, res = [res res];
    case 2
    otherwise, error('length(res)~=1 or 2');
end
if length(size(d))~=2, error('length(size(d)~=2'); end
if size(d,2)~=length(k), error('size(d,2)~=length(k)'); end
if length(k)~=length(dcf)
    error('length(k)~=length(dcf)'); 
end

N01 = round(res(1));
N02 = round(res(2));
N11 = 2*ceil(res(1)*oversamp/2); % round to even value (fftshift=iffthift)
N12 = 2*ceil(res(2)*oversamp/2);

if isreal(d), d = complex(d); end
if ~isa(d,'double'), d = double(d); end
if isreal(k)
    if length(size(k))~=3, error('length(size(k))~=3'); end
    if size(k,1)~=1, error('size(k,1)~=1'); end
    k = k(1,:,1) +1i*k(1,:,2);
end
traj = double([real(k(:)).' ; imag(k(:)).' ; dcf(:).']); % adjust to grid_interp input
deapod_im = deapodize([N11 N12]);      % deapodisation kernel
% mm = complex(zeros(N01,N02,nreps));    % pre-allocate
ctr1_index = N11/2+1;                  % centre indices
ctr2_index = N12/2+1;
ind1 = int64(ctr1_index-N01/2:ctr1_index+N01/2-1);
ind2 = int64(ctr2_index-N02/2:ctr2_index+N02/2-1);

% % loop over coils: major speedup possible by moving loop into grid_interp
% for lc=1:nreps
%     % M = grid_interp(double(d(lc,:)),traj,[N11 N12],[N11 N12]);
%     M = grid_interp_mex(d(lc,:).',traj,[N11 N12],[N11 N12]);
%     m = deapod_im.*fftshift(ifft2(ifftshift(M)));
%     m = m(ind1,ind2);
%     mm(:,:,lc) = m*oversamp^2;  % inverse FFT scaled by 1/N => correct oversampling
% end
M = grid_interp_mex(d.',traj,[N11 N12],[N11 N12]);
m = bsxfun(@times,deapod_im,fftshift(fftshift(ifft2(ifftshift(ifftshift(M,1),2)),2),1));
mm = m(ind1,ind2,:)*oversamp^2;  % inverse FFT scaled by 1/N => correct oversampling


end      % main function gridding.m
