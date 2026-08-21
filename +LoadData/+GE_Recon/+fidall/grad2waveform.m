function wf = grad2waveform(fname,h,verb,real_k)
%GRAD2WAVEFORM  Load fidspiral.grad & convert to recon_spiral.m waveform file
%      wf = grad2waveform(fname,h,verb,real_k)
%   fname   Input k-space filename                       ('fidspiral.grad')
%       h   Data header structure: check against fidspiral.grad pars   ([])
%    verb   0=off, 1=verbose, 2=verbose+plotting                        (0)
%  real_k   Real wf.k values; 2nd dimension in size(k,3)             (true)
%      wf   Waveform structure
%
% 12/2024 Rolf Schulte
if (nargin<1) && (nargout<1), help(mfilename); return; end


%% iput args
if ~exist('fname','var'),  fname = ''; end
if isempty(fname),         fname = 'fidspiral.grad'; end
if ~exist(fname,'file'),   error('fname(=''%s'') not found',fname); end
if ~exist('h','var'),      h = []; end
if ~exist('verb','var'),   verb = []; end
if isempty(verb),          verb = 1; end
if ~exist('real_k','var'), real_k = []; end
if isempty(real_k),        real_k = true; end

gdt = 4d-6;                       % gradient update time [s]


%% reading fidspiral.grad data
fid    = fopen(fname,'r','b');    % big endian by default (written on ice)
ng     = fread(fid,1,'int');      % #grad points per interleave
nk     = fread(fid,1,'int');      % #acquired samples per interleave
nviews = fread(fid,1,'int');      % #spiral arms (aka interleaves)
mtx    = fread(fid,1,'int');      % reconstruction matrix size

if (ng<1) || (nk<1) || (nviews<1) || (mtx<1) || (mtx>2048) || (nviews>64)
    warning('Strange numbers: ng=%d nk=%d nviews=%d mtx=%d: trying little endian',...
        ng,nk,nviews,mtx);
    fclose(fid);
    fid = fopen(fname,'r','l');   % little endian (written on host)
    ng     = fread(fid,1,'int');
    nk     = fread(fid,1,'int');
    nviews = fread(fid,1,'int');
    mtx    = fread(fid,1,'int');
end
bw     = fread(fid,1,'float');    % half sampling BW       [kHz]
fov    = fread(fid,1,'float');    % field of view          [mm]
GAM    = fread(fid,1,'float');    % gyromagnetic ratio     [Hz/G]
gmax   = fread(fid,1,'float');    % max gradient strength  [G/cm]
gx = fread(fid,ng,'short');
gy = fread(fid,inf,'short');
if ~feof(fid), warning('~feof(fid)'); end
if length(gy)~=ng
    warning('length(gy)(=%g)~=ng(=%g)',length(gy),ng);
    gy = gy(1:ng);
end
fclose(fid);


%% info
if verb>0
    fprintf('grad2waveform.m\n');
    fprintf('ng=%g nk=%g nviews=%g mtx=%g bw=%g fov=%g GAM=%g gmax=%g\n',...
        ng,nk,nviews,mtx,bw,fov,GAM,gmax);
end


%% convert to SI units
bw  = bw*2000;                    % sampling BW            [Hz]
% fov = fov*1d-3;                   % field of view          [m]
gamma = 2*pi*GAM*1d4;             % gyromagnetic ratio [rad/(s T)]
gmax   = gmax*1d-2;               % max gradient strength  [T/m]

if mtx<1, warning('mtx(=%g)<1; setting to 100',mtx); mtx = 100; end
kdt = round(1d6/bw)*1d-6;         % sampling dwell time    [s]
res = fov*1d-3/mtx;               % pixel size=resolution  [m]


%% assembling & scaling grad
grad = gx(:).' +1i*gy(:).';       % (1,pts); gx=real; gy=imag
grad = grad*gmax/2^15;


%% interpolating to kdt
tgrad = (0:(ng-1))*gdt;           % gradient time [s]
tksp  = (0:(nk-1))*kdt;           % ksp/acq time [s]
grad = interp1(tgrad,grad,tksp);
grad = grad/2/pi*res*gamma*kdt;   % for k-space normalised to +-0.5


%% convert to k-space
ksp = cumsum(grad,2);             % integrate for k-space


%% density compensation function (similar results as voronoi_area.m)
dcf = abs(real(ksp).*real(grad) + imag(ksp).*imag(grad));
dcf = dcf/max(dcf);               % scaling unknown; normalising


%% adding interleaves
if nviews>1
    ksp = bsxfun(@times,exp(1i*2*pi*(0:(nviews-1))/nviews).',ksp);
    if verb>1, plot(ksp.'); grid on; end
    ksp = reshape(ksp.',[1 nk*nviews]);
else
    if verb>1, plot(ksp.'); grid on; end
end


%% wf fields
if real_k
    wf.k = cat(3,real(ksp),imag(ksp));
else
    wf.k = ksp;
end
wf.dcf  = repmat(dcf,[1 nviews]);
wf.ind  = true(nviews,nk);
wf.t    = repmat(tksp,[nviews 1]);
wf.npix = mtx;
wf.mtx  = mtx;


%% re-check with header
if ~isempty(h)
    % ng,nk,nviews,mtx,bw,fov,GAM,gmax
    if (nk~=h.image.dim_X)
        warning('nk(=%g) ~= h.image.dim_X(=%g)',nk,h.image.dim_X);
    end
    if (nviews~=h.image.dim_Y)
        warning('nviews(=%g) ~= h.image.dim_Y(=%g)',nviews,h.image.dim_Y);
    end
    % if (mtx ~= h.rdb_hdr.im_size)
    %     warning('mtx(=%g) ~= h.rdb_hdr.im_size(=%g)',mtx,h.rdb_hdr.im_size);
    % end
    if abs(bw-h.rdb_hdr.bw*2d3)>1
        warning('bw(=%g) ~= h.rdb_hdr.bw*2d3(=%g)',bw,h.rdb_hdr.bw*2d3);
    end
    if abs(fov-h.rdb_hdr.fov)>0.1
        warning('fov(=%g)~=h.rdb_hdr.fov(=%g)',fov,h.rdb_hdr.fov);
    end
    gg = gyrogamma(h.image.specnuc);
    if abs(gamma-gg)>1d4
        warning('gamma(=%g) ~= gyrogamma(%d)(=%g)',gamma,h.image.specnuc,gg);
    end
    if (abs(res*1d3-h.image.user48)>1d-3)
        warning('res(=%g)~=h.image.user48(=%g)',res*1d3,h.image.user48);
    end
end


%% info
kmax = max(abs(ksp(:)));          % max k-space extend
if kmax>0.51, warning('kmax(=%g)>0.51',kmax); end
if verb>0, fprintf('kmax = %g\n',kmax); end



end      % grad2waveform.m
