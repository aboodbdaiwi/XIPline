function wf = calc_spiral_waveform(h,verb)
%CALC_SPIRAL_WAVEFORM Calculate waveform structure for spiral.e
% k-space via kspace_calcMEX
%   wf = calc_spiral_kspace(h,verb)
%    h   Header structure
% verb   Verbose                                                     (true)
%   wf   Waveform structure
%
% 12/2024 Rolf Schulte
if (nargin<1), help(mfilename); return; end

if ~exist('verb','var'),   verb = []; end
if isempty(verb),          verb = true; end


%% extract info from header
fprintf('calc_spiral_waveform.m\n');
fov   = h.rdb_hdr.grid_scan_fov;
t_delta = h.rdb_hdr.grid_a2d_time;
gmax    = h.rdb_hdr.grid_max_gradient;
slewmax = h.rdb_hdr.grid_max_slew;
% mtx    = ceil(h.rdb_hdr.grid_scan_fov/h.image.pixelSizeX/2)*2;
mtx    = floor(h.rdb_hdr.grid_scan_fov/h.image.pixelSizeX);
mtx1    = (h.rdb_hdr.grid_scan_fov/h.image.pixelSizeX);
nk     = h.image.dim_X;
nviews = h.image.dim_Y;
bw = h.image.vbw*2d3;
mtx_reco = h.rdb_hdr.im_size;
if (nk ~= h.rdb_hdr.da_xres), warning('xres'); end
if (nviews ~= (h.rdb_hdr.da_yres - 1)), warning('yres'); end


%% calculate spiral.e k-space
ksp = kspace_calcMEX(nviews,nk,fov,t_delta,gmax,slewmax);
ksp = ksp(1:nviews*nk*3,1);            % crop waveform
ksp = reshape(ksp,[3,nviews*nk]);      % reshape

%% k-space scaling
% kscale = 256/mtx1;
kscale = 256/mtx;


%% density compensation function
% analytical calculation (specific to spiral design)
% same formula, but more accurate than dcf from kspace_calcMEX
grad = diff([zeros(3,1) ksp],[],2);
dcf = abs(ksp(1,:).*grad(1,:) + ksp(2,:).*grad(2,:));
dcf = dcf/max(dcf);                    % scaling unknown; normalising


%% populate wf structure
wf.k = cat(3,ksp(1,:),ksp(2,:))*kscale;
% wf.dcf  = ksp(3,:);
wf.dcf  = dcf;
wf.ind  = true(nviews,nk);             % index list for reconstruction
wf.t = repmat((0:(nk-1))/bw,[nviews 1]);    % time [s]
wf.npix = mtx;                         % acq. matrix size
wf.mtx  = mtx;


%% print info
if verb
    kmax = max(sqrt(sum(wf.k(1,:,1).^2 + wf.k(1,:,2).^2,3)));
    fprintf('mtx = %g; kmax = %g\n',mtx,kmax);
    if kmax>0.51, warning('kmax(=%g)>0.51',kmax); end
end

end      % calc_spiral_waveform.m
