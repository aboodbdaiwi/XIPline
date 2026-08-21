function [k,dcf,npts_low,npts_high,t,xtra]=zte_header2kspace(h,verb,...
    mode,burst_half,max_waspi)
%ZTE_HEADER2KSPACE Create ZTE 3dradial k-space trajectory for recon_zte.m
%[k,dcf,npts_low,npts_high,t,xtra] = ...
%    zte_header2kspace(h,verb,mode,burst_half,max_waspi)
%                                                                 (default)
%         h  p-file header
%      verb  verbose mode (0=off; 1=verbose; 2=verbose+plotting)        (1)
%      mode  k-space trajectory mode                                    (0)
%            0=truncate centre                              (ZTE/WASPI)
%            1=full trajectory (echo-top; centre of k)      (BURST)
%            2=full traj. (max-echo; full acq.; xtra=0)     (BURST)
%            3=full traj. (echo-top); discard lowres data   (BURST)
%            4=full traj. (max-echo); discard lowres data   (BURST)
%burst_half  half k-space data (0=full; 1=upper half (hires); 2=symmetric
% max_waspi  Maximise WASPI/crop highres data for slow TR switches   (false)
%
%         k  k-space trajectory (range +-0.5)
%       dcf  density compensation function (analytically calculated)
%  npts_low  #points used for recon - inner part of k-space
% npts_high  #points used for recon - outer part of k-space
%         t  time [s]
%      xtra  Data samples in front for discarding
%
% 8/2020  Rolf Schulte
if (nargin<1), help(mfilename); return; end

if ~exist('verb','var'), verb = []; end
if isempty(verb),        verb = 1; end
if verb>1, plt = true; 
else,      plt = false;
end
if ~exist('mode','var'), mode = []; end
if isempty(mode),        mode = 0; end
if ~exist('burst_half','var'), burst_half = []; end
if isempty(burst_half),       burst_half = false; end
if ~exist('max_waspi','var'), max_waspi = []; end
if isempty(max_waspi),        max_waspi = false; end


%% checks
if h.rdb_hdr.noncart_dual_traj~=1
    warning('waspi_flag=%g (no separate low+high res)',...
        h.rdb_hdr.noncart_dual_traj);
end


%% HEADER info
over_samp = h.rdb_hdr.oversamplingfactor; % oversampling along readout 
xres = h.rdb_hdr.da_xres;         % #acq points
yres = h.rdb_hdr.da_yres;         % #exc stored in yres and nslc
nslc = h.rdb_hdr.nslices;
xtra = h.rdb_hdr.feextra;         % #acq points at FID start to discard
npts = xres-xtra;                 % #acq points used for recon 
fov = h.rdb_hdr.fov;              % field in view [mm]
slthick = h.image.slthick;        % slice thickness = voxel size [mm]
mtx = floor(fov/slthick+0.5);     % nominal matrix size of acquisition
nrcv = (h.rdb_hdr.dab(2)-h.rdb_hdr.dab(1))+1;   % #coils
nspokes_low  = h.rdb_hdr.nspokes_lowres;        % #spokes in k-space centre
nspokes_high = h.rdb_hdr.nspokes_highres;       % #spokes outside
waspi_factor = h.rdb_hdr.noncart_traj_kmax_ratio; % kmax_high/kmax_low
if isfield(h.image,'spokesPerSeg')
    spokesPerSeg = h.image.spokesPerSeg;
end
dt = 1/(2d3*h.rdb_hdr.bw*h.rdb_hdr.oversamplingfactor); % sampling dwell time
if isfield(h.rdb_hdr,'noncart_traj_mode')
    trajectory_mode = h.rdb_hdr.noncart_traj_mode;
else
    trajectory_mode = 0;
end
if trajectory_mode>0
    error('trajectory_mode(=%g)>0: not yet implemented',trajectory_mode);
end

if (xres>(2*xtra+2)) && (mode>0)
    warning('#acq pts (=xres=%g)>(2*xtra(=%g)+2): increase burst_prerf_factor cv17(=%g)',...
        xres,xtra,h.rdb_hdr.user17);
end
if (xres<(2*xtra-4))
    warning('#acq pts (xres=%g)<(2*xtra(=%g)-4): reduce burst_prerf_factor cv17(=%g)',...
        xres,xtra,h.rdb_hdr.user17);
end



%% merge bounds for ZTE lowres/highres data
[~,mi] = min(abs(h.rdb_hdr.bw-[15.625 31.25 62.5 125]));
switch mi
    case 1, merge_start=2;  merge_end=3;    % 15.625
    case 2, merge_start=2;  merge_end=3;    % 31.25
    case 3, merge_start=6;  merge_end=8;    % 62.5
    case 4, merge_start=12; merge_end=14;   % 125
    otherwise, error('no minimum found (bw=%g)',h.rdb_hdr.bw);
end
% merge_start = h.rdb_hdr.noncart_traj_merge_start;
% merge_end = h.rdb_hdr.noncart_traj_merge_end;
if over_samp~=2
    merge_start = floor(merge_start*over_samp/2+0.5);
    merge_end = floor(merge_end*over_samp/2+0.5);
end
if ((burst_half==1) && (mode==0))
    merge_end = min([merge_end*2 floor(npts/waspi_factor)]);
end
if (max_waspi && (mode==0))
    if verb
        fprintf('maximum lowres & minimum highres data\n'); 
    end
    merge_end = floor(npts/waspi_factor);
    merge_start = merge_end-1;
end


%% k-space trajectory
% direction of spokes
switch burst_half
    case 0
        KTrajHigh = sub_zte_trajectory(nspokes_high,false);
    case 1
        KTrajHigh = sub_zte_trajectory_half(nspokes_high,false);
    case 2
        KTrajHigh = sub_zte_trajectory_symmetric(nspokes_high,false);
    otherwise
        error('burst_half~=0,1 or 2');
end

KTrajLow  = sub_zte_trajectory(nspokes_low,false);
KTrajLow  = KTrajLow/2;
KTrajHigh = KTrajHigh/2;


switch mode
    case 0    % truncate centre                            (ZTE/WASPI)
        npts_low  = merge_end*waspi_factor;      % #points used for recon
        if (npts_low+xtra)>xres
            warning('npts_low(=%g)+xtra(=%g))>npts(=%g); reducing npts_low to %g',...
                npts_low,xtra,xres,xres-xtra);
            npts_low = xres-xtra;
        end
        npts_high = npts-merge_start;
    case 1    % full trajectory (echo-top; centre of k)    (BURST)
        npts_low = npts;                         % #points used for recon
        npts_high = npts;
        merge_start = 0;
    case 2    % full traj. (max-echo; full acq.; xtra=0)   (BURST)
        npts_low = xres;
        npts_high = xres;
        merge_start = 0;
    case 3    % full traj. (echo-top); discard lowres data (BURST)
        npts_low = 0;
        npts_high = npts;
        merge_start = 0;
        nspokes_low = 0;                         % skip low-res data
    case 4    % full traj. (max-echo); discard lowres      (BURST)
        npts_low = 0;
        npts_high = xres;
        merge_start = 0;
        nspokes_low = 0;                         % skip low-res data
    otherwise
        error('mode (=%g) unknown',mode);
end

% if npts_low>npts, 
%     warning('npts_low(=%g)>npts(=%g); setting equal',npts_low,npts);
%     npts_low = npts;
% end
if ((mode~=2) && (mode~=4))
    ks_low  = (0:npts_low-1).'/npts/waspi_factor;    % k-space for single spoke
    ks_high = (merge_start:npts-1).'/npts;
    t_low  = (0:npts_low-1)*dt;
    t_high = (merge_start:npts-1)*dt;
else
    ks_low = ((0:xres-1)-xtra).'/npts/waspi_factor;
    ks_high = ((0:xres-1)-xtra).'/npts;
    t_low  = ((0:xres-1)-xtra)*dt;
    t_high = t_low;
    % xtra = 0;
end
k_low = zeros(nspokes_low*npts_low,3);           % init for full k-space
k_high = zeros(nspokes_high*npts_high,3);

for l=1:nspokes_low
    ii = (l-1)*npts_low+(1:npts_low);
    k_low(ii,:) = ks_low*KTrajLow(:,l).';
end
for l=1:nspokes_high
    ii = (l-1)*npts_high+(1:npts_high);
    if ii(end)>(npts_high*nspokes_high)
        warning('ii(end)(=%g)>%g',ii(end),(npts_high*nspokes_high)); 
    end
    k_high(ii,:) = ks_high*KTrajHigh(:,l).';
end
k = [k_low ; k_high];   % size(k) = (nk,3)
k = shiftdim(k,-1);     % size(k) = (1,nk,3)


%% density compensation function
% 4/3*pi*(r1^3-r0^3)/n
dcfs_low = (ks_low(2:end).^3-ks_low(1:end-1).^3)/nspokes_low;
dcfs_low = [dcfs_low;dcfs_low(end)];
if mode==2, dcfs_low = dcfs_low/10; end
if (mode==0) && (burst_half==1)
    dcfs_low = dcfs_low/2;
end
dcf_low  = repmat(dcfs_low,[1 nspokes_low]);
dcf_low  = reshape(dcf_low,[1 npts_low*nspokes_low]);
dcfs_high = (ks_high(2:end).^3-ks_high(1:end-1).^3)/nspokes_high;
dcfs_high=[dcfs_high;dcfs_high(end)];
if ((mode==2)||(mode==4))
    if xres>2*xtra
        dcfs_high(1:2*xtra) = dcfs_high(1:2*xtra)/2;
    else
        dcfs_high = dcfs_high/2;
    end
end
dcf_high = repmat(dcfs_high,[1 nspokes_high]);
dcf_high = reshape(dcf_high,[1 npts_high*nspokes_high]);
dcf = 1/6*pi*mtx^3*[dcf_low dcf_high];
% scaling: max(ks_high)=1 -> divide by 8, so 4/3*pi/8


%% time
t = [repmat(t_low,[1 nspokes_low]) , repmat(t_high,[1 nspokes_high])];


%% print info
if verb
    fprintf('BW = +- %g [kHz]; dt = %g [us]\n',h.rdb_hdr.bw,dt*1d6);
    if exist('spokesPerSeg','var')
        fprintf('spokesPerSeg = %g; nseg_low = %g; nseg_low = %g\n',...
            spokesPerSeg,nspokes_low/spokesPerSeg,nspokes_high/spokesPerSeg);
    end
    fprintf('merge_start = %g; merge_end = %g\n',merge_start,merge_end);
    fprintf('xres = %g; yres = %g; nslc = %g; nrcv = %g\n',xres,yres,nslc,nrcv);
    fprintf('xtra = %g; npts = xres-xtra = %g\n',xtra,npts);
    fprintf('waspi_factor  = %g; over_samp = %g\n',waspi_factor,over_samp);
    fprintf('npts_low  = %g; npts_high = %g; npts = %g\n',npts_low,npts_high,npts);
    fprintf('nspokes_low  = %g; nspokes_high = %g; nspokes = %g\n',...
        nspokes_low,nspokes_high,nspokes_low+nspokes_high);
    fprintf('(yres-1)*nslc = %g\n',(yres-1)*nslc);
    fprintf('kmax = %g\n',max(sqrt(sum(k.^2,3))));
end


%% plotting
if plt
    figure(11); clf; plot(KTrajLow.');  title('KTrajLow');
    figure(12); clf; plot(KTrajHigh.'); title('KTrajHigh');
    figure(13); clf; plot(squeeze(k)); title('k');
end

if ((mode==2)||(mode==4)), xtra = 0; end

end      % zte_header2kspace.m


%% sub-function to create actual radial trajectory
function vec = sub_zte_trajectory(n, plotOn)
%   n, the number of spokes, needs to be even
vec = zeros(3,n);
nOver2 = n/2;
dt = 1/nOver2;

% calculate trajectory "velocity" squared
c2 = n*4*pi;

% calculate a hemisphere of points
phiOld = 0.0;
vec(1,nOver2+1) = 1;
vec(2,nOver2+1) = 0;
vec(3,nOver2+1) = 0;
for m=1:nOver2-1
    cT = m*dt;
    cT2 = cT*cT;
    sT2 = 1-cT2;
    sT = sqrt(sT2);
    dP = 0.5*dt*sqrt((1/sT2)*(c2-(1/sT2)));
    phiNew = phiOld+dP;
    vec(1,nOver2+1-m) = cos(phiNew)*sT;
    vec(2,nOver2+1-m) = sin(phiNew)*sT;
    vec(3,nOver2+1-m) = cT;
    phiOld = phiNew;
end
vec(1,1) = 0;
vec(2,1) = 0;
vec(3,1) = 1;
% reflect the second hemisphere
for m=1:nOver2-1
    vec(1,m+nOver2+1) = vec(1,nOver2+1-m);
    vec(2,m+nOver2+1) = -vec(2,nOver2+1-m);
    vec(3,m+nOver2+1) = -vec(3,nOver2+1-m);
end
if plotOn
    sel = 1:n;
    figure, plot3(vec(1,sel),vec(2,sel),vec(3,sel)); axis equal;
end

end  % sub_zte_trajectory


%% sub-function to create half radial trajectory
function vec = sub_zte_trajectory_half(n, plotOn)
%   n, the number of spokes, needs to be even
% n = n*2;
vec = zeros(3,n);
% nOver2 = n/2;
nOver2 = n;
dt = 1/nOver2;

% calculate trajectory "velocity" squared
c2 = n*4*pi;
% c2 = c2/2;

% calculate a hemisphere of points
phiOld = 0.0;
% vec(1,nOver2+1) = 1;
% vec(2,nOver2+1) = 0;
% vec(3,nOver2+1) = 0;
for m=1:nOver2-1
    cT = m*dt;
    cT2 = cT*cT;
    sT2 = 1-cT2;
    sT = sqrt(sT2);
    dP = 0.5*dt*sqrt((1/sT2)*(c2-(1/sT2)));
    phiNew = phiOld+dP;
    vec(1,nOver2+1-m) = cos(phiNew)*sT;
    vec(2,nOver2+1-m) = sin(phiNew)*sT;
    vec(3,nOver2+1-m) = cT;
    phiOld = phiNew;
end
vec(1,1) = 0;
vec(2,1) = 0;
vec(3,1) = 1;
% reflect the second hemisphere
% for m=1:nOver2
%     vec(1,m+nOver2) = -vec(1,m);
%     vec(2,m+nOver2) = -vec(2,m);
%     vec(3,m+nOver2) = -vec(3,m);
% end
if plotOn
    sel = 1:n;
    figure, plot3(vec(1,sel),vec(2,sel),vec(3,sel)); axis equal;
end

end      % sub_zte_trajectory_half


%% sub-function to create symmetric radial trajectory
function vec = sub_zte_trajectory_symmetric(n, plotOn)
%   n, the number of spokes, needs to be even

vec = zeros(3,n);
nOver2 = n/2;
dt = 1/nOver2;

% calculate trajectory "velocity" squared
c2 = n*4*pi;

% calculate a hemisphere of points
phiOld = 0.0;
vec(1,nOver2+1) = -1;
vec(2,nOver2+1) = 0;
vec(3,nOver2+1) = 0;
for m=1:nOver2-1
    cT = m*dt;
    cT2 = cT*cT;
    sT2 = 1-cT2;
    sT = sqrt(sT2);
    dP = 0.5*dt*sqrt((1/sT2)*(c2-(1/sT2)));
    phiNew = phiOld+dP;
    vec(1,nOver2+1-m) = cos(phiNew)*sT;
    vec(2,nOver2+1-m) = sin(phiNew)*sT;
    vec(3,nOver2+1-m) = cT;

    
    phi2 = phiNew+dP/2;
    vec(1,nOver2+1+m) = -cos(phi2)*sT;
    vec(2,nOver2+1+m) = -sin(phi2)*sT;
    vec(3,nOver2+1+m) = -cT;
    
    phiOld = phiNew;
end
vec(1,1) = 0;
vec(2,1) = 0;
vec(3,1) = 1;

if plotOn
    sel = 1:n;
    figure, plot3(vec(1,sel),vec(2,sel),vec(3,sel)); axis equal;
end

end      % sub_zte_trajectory_symmetric

