function epsi = header2epsi(h,verb,calc_traj)
%HEADER2EPSI Generate EPSI structure from pressepsi/fidepsi cv parameters 
%            saved in p-file header
% epsi = header2epsi(h,verb,calc_traj)
%        h  GE p-file header with epsi parameters (in h.rdb_hdr.maxcoefXX)
%     verb  Verbose mode                                   (default=true)
%calc_traj  Calculate gradient and k-space trajectories    (default=false)
%     epsi  Structure for EPSI reconstruction
%           ind  Index list for data selection
% 1/2018 Rolf Schulte
if (nargin<1), help(mfilename); return; end;

if ~exist('verb','var'),   verb = []; end
if isempty(verb),          verb = true; end
if ~exist('calc_traj','var'),   calc_traj = []; end
if isempty(calc_traj),          calc_traj = false; end


%% misc parameters
gdt            = 4;                    % gradient update time          [us]
epsi.mode      = h.rdb_hdr.user22;     % epsi mode: 0=off; 1=flyback; 2=bidirectional
epsi.intlv     = h.rdb_hdr.maxcoef1a;  % #EPSI interleaves (time shifted excitations)
% kdt          = h.rdb_hdr.maxcoef1b;  % sampling dwell time           [us]
kdt            = 1d3/(2*h.rdb_hdr.bw); % sampling dwell time           [us]
epsi.lobepts   = h.rdb_hdr.maxcoef1c;  % #acq samples of single EPSI lobe (readout+flyback+ramps)
epsi.res       = h.rdb_hdr.maxcoef1d;  % #spatial samples along EPSI direction
epsi.dir       = h.rdb_hdr.maxcoef2a;  % Logical axis for EPSI readout gradient (0=x,1=y,2=z)
epsi.loc       = h.rdb_hdr.maxcoef2b;  % off-iso location in epsi_dir  [mm]
epsi.pw_gepsi  = h.rdb_hdr.maxcoef2c;  % dur flat EPSI part (readout)  [us]
epsi.pw_gepsia = h.rdb_hdr.maxcoef2d;  % dur ramps to flat part        [us]
epsi.pw_gflyb  = h.rdb_hdr.maxcoef3a;  % dur flat part of flyback      [us]
epsi.pw_gflyba = h.rdb_hdr.maxcoef3b;  % dur of flyback ramps          [us]
epsi.a_gepsi   = h.rdb_hdr.maxcoef3c;  % amplitude of EPSI readout   [G/cm]
epsi.a_gflyb   = h.rdb_hdr.maxcoef3d;  % amplitude of flyback        [G/cm]
epsi.samplepts = h.rdb_hdr.frame_size; % #acquired samples per excitation
% specpts        = h.image.user1;      % #prescribed spectral samples
specpts        = h.rdb_hdr.user1;      % #acquired spectral samples
% specwidth      = h.rdb_hdr.spectral_width; % spectral BW [Hz]
specwidth      = h.rdb_hdr.user0;      % spectral BW [Hz]
ispecpts       = specpts/epsi.intlv;   % #spectral points for one EPSI interleave
if (ispecpts-floor(ispecpts))>0, 
    warning('ispecpts-floor(ispecpts)(=%g)>0; rounding up to %g',...
        (ispecpts-floor(ispecpts)),ceil(ispecpts));
    ispecpts = ceil(ispecpts);
end


%% misc checks
if ~any(epsi.mode==[1 2])
    warning('h.image.user22(=%g)~=1 or 2; really EPSI?',epsi.mode);
end
if epsi.intlv<1, warning('epsi.intlv(=%g) < 1',epsi.intlv); end
if epsi.intlv>4, warning('epsi.intlv(=%g) > 4',epsi.intlv); end

% sampling dwell time (1/BW)
if abs(kdt-floor(kdt))>1d-5, 
    warning('epsi.kdt (=%g) not a natural number; rounding');
    % kdt = round(kdt);
end
if kdt~=round(kdt), kdt = round(kdt); end
    
if abs(1d6*epsi.intlv/(epsi.lobepts*kdt)-specwidth)>1d-1,
    warning('abs(1d6*epsi.intlv(=%g)/(epsi.lobepts(=%g)*kdt(=%g))(=%g) - specwidth(=%g))>1d-1',...
        epsi.intlv,epsi.lobepts,kdt,1d6*epsi.intlv/(epsi.lobepts*kdt),specwidth);
end

if ((abs(epsi.loc)>1d-7)&&(epsi.mode==2))
    warning('bidirectional epsi should use epsi_iso=2');
end


%% number of spectral samples
nk_gepsi  = epsi.pw_gepsi/kdt;          % on flat part for sampling 
nk_gepsia = epsi.pw_gepsia/kdt;         % ramps of flat epsi part


%% index list for data for reconstruction
switch epsi.mode
    case 1              % flyback epsi (no ramp sampling)
        ind = false(1,epsi.samplepts);
        for l=1:ispecpts
            li = (l-1)*epsi.lobepts + round(nk_gepsia) + (1:nk_gepsi);
            ind(1,li) = true;
        end
    case 2              % bidirectional epsi with ramp sampling
        ind = true(1,epsi.samplepts);
    case 3              % bidirectional epsi no ramp sampling
        warning('not yet implemented');
end

if size(ind,2)>epsi.samplepts, 
    warning('size(ind,2)(=%g) > epsi.samplepts(=%g): truncating ind',...
        size(ind,2),epsi.samplepts);
    ind = ind(1,1:epsi.samplepts);
end
if size(ind,2)<epsi.samplepts, 
    warning('size(ind,2)(=%g) < epsi.samplepts(=%g): filling up',...
        size(ind,2),epsi.samplepts);
    ind = [ind false(1,(epsi.samplepts-size(ind,2)))];
end


%% calculate phase function to compensate for phase error
% calculate phase difference between grad lobes due to off-isocentre 
% EPSI grid shifted using setfrequency (instead of feeding gx into OMEGA)
% -> shift in spectral f0-frequency
if abs(epsi.loc)>1d-7,
    epsi.f0 = (epsi.loc*1d-3)*(epsi.a_gepsi*10d-3)*abs(gyrogamma(h.image.specnuc))/2/pi;
    epsi.dphi = 2*pi*epsi.f0/specwidth;
    pha = repmat((0:specpts-1)*epsi.dphi,[epsi.res 1]);    
    % tilt k-space lines back; compensate for time evolution
    % See mrm54_1286_2005_cunningham.pdf equation 3
    % How to verfiy??? Effect seems small...
    % pha = pha + 2*pi*(0:epsi.res-1).'*(0:specpts-1)*kdt*1d-6*specwidth/specpts;
    % pha = pha.';
    epsi.phafu = exp(-1i*pha(:));
else
    epsi.phafu = [];
end

%% put into epsi structure
epsi.npts = [epsi.res specpts];             % #pts (spatial,spectral)
epsi.bw  = [h.rdb_hdr.bw*2d3 specwidth];    % BW (spat,spect) [Hz]
epsi.ind = ind;
epsi.gdt = gdt;
epsi.kdt = kdt;


%%
if sum(ind)~=epsi.res*ispecpts,
    warning('sum(ind)(=%g) ~= epsi.res*ispecpts (=%g*%g=%g)',...
        sum(ind),epsi.res,ispecpts,epsi.res*ispecpts);
end

%% print some info
if verb,
    fprintf('epsi.mode      = %g      \t(0=off; 1=flyback; 2=bidirectional)\n',epsi.mode);
    fprintf('epsi.intlv     = %g      \t(#EPSI interleaves)\n',epsi.intlv);
    fprintf('epsi.kdt       = %g [us] \t(sampling dwell time)\n',epsi.kdt);
    fprintf('epsi.lobepts   = %g      \t(#samples of single EPSI lobe (readout+flyback+ramps)\n',epsi.lobepts);
    fprintf('epsi.pw_gepsi  = %g [us] \t(time on plateau)\n',epsi.pw_gepsi);
    fprintf('epsi.pw_gepsia = %g [us] \t(ramp time)\n',epsi.pw_gepsia);
    fprintf('epsi.samplepts = %g    \t(#acq samples per excitation)\n',epsi.samplepts);
    fprintf('epsi.ispecpts  = %g      \t(#spectral points for one EPSI interleave)\n',ispecpts);
    fprintf('epsi.res       = %g      \t(Spatial samples along EPSI readout direction)\n',epsi.res);
    fprintf('epsi.dir       = %g      \t(Logical axis for EPSI readout gradient (0=x,1=y,2=z)\n\n',epsi.dir);
    
    fprintf('user0          = %g [Hz] \t(BW)\n',h.image.user0);
    fprintf('user1          = %g      \t(samples)\n',h.image.user1);
    if specpts>h.image.user1,
        fprintf('   specpts     = %g (due to epsi_intlv rounding)\n',specpts);
    end
    fprintf('csi_dims       = %g\n',h.rdb_hdr.csi_dims);
    fprintf('xcsi           = %g\n',h.rdb_hdr.xcsi);
    fprintf('ycsi           = %g\n',h.rdb_hdr.ycsi);
    fprintf('zcsi           = %g\n',h.rdb_hdr.zcsi);
    fprintf('sum(ind)       = %g\n\n',sum(ind));
end


%% calculate trajectory
if calc_traj
    if epsi.mode~=1
        warning('calc_traj only implemented for flyback');
        return
    end
    
    %% calculate gradients
    % number of gradient points
    ng_gepsi  = epsi.pw_gepsi/gdt;          % on flat part for sampling
    ng_gepsia = epsi.pw_gepsia/gdt;         % ramps of flat epsi part
    ng_gflyb  = epsi.pw_gflyb/gdt;          % flat part of flyback
    ng_gflyba = epsi.pw_gflyba/gdt;         % flyback ramps
    ng_pts = ng_gepsi + 2*ng_gepsia + ng_gflyb + 2*ng_gflyba;
    if rem(ng_gepsi,1),  warning('ng_gepsi not a natural number'); end
    if rem(ng_gepsia,1), warning('ng_gepsia not a natural number'); end
    if rem(ng_gflyb,1),  warning('ng_gflyb not a natural number'); end
    if rem(ng_gflyba,1), warning('ng_gflyba not a natural number'); end
    
    if verb,
        fprintf('ng_gepsi=%g; ng_gepsia=%g; ng_gflyb=%g; ng_gflyba=%g; ng_pts=%g\n',...
            ng_gepsi,ng_gepsia,ng_gflyb,ng_gflyba,ng_pts);
    end
    
    if (ng_pts ~= (epsi.lobepts*kdt/gdt))
        warning('ng_pts(=%g) ~= (epsi.lobepts(=%g)*kdt(=%g)/gdt(=%g)(=%g)',...
            ng_pts,epsi.lobepts,kdt,gdt,epsi.lobepts*kdt/gdt);
    end
    
    % adjust + check amplitudes
    a_gepsi = epsi.a_gepsi*10d-3;         % grad amplitude epsi    [G/cm -> T/m]
    a_gflyb = epsi.a_gflyb*10d-3;         % grad amplitude flyback [G/cm -> T/m]
    if a_gflyb>0,
        warning('a_gflyb(=%g) > 0; inverting sign',a_gflyb);
        a_gflyb = -a_gflyb;
    end
    
    % calculate gradient sublobes
    gepsia = linspace(0,1,ng_gepsia+1);
    gepsid = linspace(1,0,ng_gepsia+1);
    gepsi  = a_gepsi*[gepsia ones(1,ng_gepsi-1) gepsid];
    gflyba = linspace(0,1,ng_gflyba+1);
    gflybd = linspace(1,0,ng_gflyba+1);
    gflyb  = a_gflyb*[gflyba ones(1,ng_gflyb-1) gflybd];
    gintlv = [gepsi(1:end-1) gflyb(1:end-1)];
    
    if length(gintlv)~=ng_pts,
        warning('length(gintlv)(=%g) ~= ng_pts(=%g)',...
            length(gintlv),ng_pts);
    end
    
    area_gepsi = sum(gepsi)*gdt*1d-6;
    area_gepsipre = area_gepsi/2;   % area of prewinder [T/m*s]

    % print info
    if verb,
        fprintf('area_gepsi    = %g [T/m*s]\n',area_gepsi);
        fprintf('area_gflyb    = %g [T/m*s]\n',sum(gflyb)*gdt*1d-6);
        fprintf('area_gintlv   = %g [T/m*s]\n',sum(gintlv)*gdt*1d-6);
        fprintf('area_gepsipre = %g [T/m*s]\n\n',area_gepsipre);
    end
    
    % checks
    if abs(sum(gintlv))>1d-5, warning('gradient moments not nulled'); end
    if abs(area_gepsi-area_gepsipre/2)>1d-5, warning('prewinder'); end
    
    % assemble sublobes into full gradient waveform
    grad = zeros(1,epsi.samplepts*kdt/gdt);
    for l=1:ispecpts
        li = (l-1)*ng_pts + (1:ng_pts);
        grad(1,li) = gintlv;
    end
    
    %% convert grad to k-space    
    ktmp = gyrogamma('1h')*(gdt*1d-6*cumsum(grad)-area_gepsipre)/2/pi*epsi.res;
    dki = kdt/gdt;
    if dki>1, ktmp = ktmp(1:dki:length(ktmp)); end
    if abs(dki-0.5)<1d-10,
        ktmp = interp1((1:length(ktmp)),ktmp,(1:0.5:length(ktmp)+0.5),'spline');
    end
    k = zeros(1,size(ktmp,2),2);
    k(1,:,1) = ktmp;
    k(1,:,2) = (0:size(ktmp,2)-1)*kdt*1d-6;
    
    
    %% assign to epsi structure
    epsi.grad = grad;
    epsi.k = k;
    epsi.t = (0:size(grad,2)-1)*gdt;
end
