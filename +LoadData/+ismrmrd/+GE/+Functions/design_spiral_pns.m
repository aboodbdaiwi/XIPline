function [k,dcf,t,grad,ind,out]=design_spiral_pns(fov,npix,arms,kdt,...
    fname,gmax,smax,nucleus,pthresh,chronaxie,rheobase,alpha,...
    acq_round2n,do_rot_file,balanced,gdt)
%DESIGN_SPIRAL_PNS Design a spiral with delayed acq for fast CSI
%[k,dcf,t,grad,ind,out]=design_spiral_pns(fov,npix,arms,kdt,...
%    fname,gmax,smax,nucleus,pthresh,chronaxie,rheobase,alpha,...
%    acq_round2n,do_rot_file,balanced)
%
%        fov  Field of view                                       [mm]
%       npix  #pixels (Cartesian resolution after gridding)
%       arms  #spatial interleaves                         (1)
%        kdt  k-space sampling time                        (16)   [us]
%      fname  Filename: if given->write_ak_wav 
%             if logical true: generate name               (false)
%       gmax  Max gradient amp                             (50)   [mT/m]
%       smax  Max slew rate                                (200)  [T/m/s]
%    nucleus  Nucleus                                      ('1H')
%    pthresh  PNS Threshold (80=normal mode; 100=1st+2nd controlled mode
%  chronaxie  Chronaxie time constant (cfrfact [us])              [s]
%             alternatively coil XRM,XRMW,WHOLE,ZOOM,HRMW
%   rheobase  Rheobase scale factor   (cfrinf [T/s])              [T/s]
%      alpha  Effective coil length   (cfdbdtdx/y/z [cm])         [m]
%             SRmin=rheobase/alpha
%acq_round2n  Round up #acq pts to to 2^n                  (false)
%do_rot_file  Write single static spiral into gradient 
%             waveform file & rotate via vap_phiXX.fdl     (false)
%   balanced  Balancing gradient area                      (false)
%
%          k  k-space trajectory  [-0.5..0.5]
%        dcf  Density compensation function (calculated with vornoi_area)
%          t  Time  (nexc x recon-pts)
%       grad  Gradient waveform [T/m] with dt=4us
%        ind  Index (2nd dim) for k-space points on spiral (excl ramps, etc)
%        out  Output structure of wrt_wavs
%
% 10/2020 Rolf Schulte
if (nargin<1), help(mfilename); return; end


%% fixed parameters
fufa_gmax = 0.99;
fufa_smax = 0.99;
fufa_gmax_rew = 0.99;
fufa_smax_rew = 0.4;
fufa_pthresh = 0.995;
nk_off = 2;                            % skip first acq points (artefact)


%% input parameters
if ~exist('arms','var'),        arms = []; end   % #spatial interleaves
if isempty(arms),               arms = 1; end
if ~exist('kdt','var'),         kdt = []; end
if isempty(kdt),                kdt = 16; end    % k-space dwell time [us]
if kdt<2, warning('kdt(=%g)<2',kdt); input(''); end
if ~exist('fname','var'),       fname = []; end   % for waveform file
if ~exist('gmax','var'),        gmax = []; end
if isempty(gmax),               gmax = 50; end    % [mT/m]
if ~exist('smax','var'),        smax = []; end
if isempty(smax),               smax = 200; end   % [T/m/s]
if ~exist('nucleus','var'),     nucleus = []; end
if isempty(nucleus),            nucleus = '1H'; end

if ~exist('pthresh','var'),     pthresh = []; end
if ~exist('chronaxie','var'),   chronaxie = []; end
if ~exist('rheobase','var'),    rheobase = []; end
if ~exist('alpha','var'),       alpha = []; end

if ~exist('acq_round2n','var'), acq_round2n = []; end
if isempty(acq_round2n),        acq_round2n = false; end
if ~exist('do_rot_file','var'), do_rot_file = []; end
if isempty(do_rot_file),        do_rot_file = false; end
if ~exist('balanced','var'),    balanced = []; end
if isempty(balanced),           balanced = false; end
if ~exist('gdt','var'),         gdt = []; end
if isempty(gdt),                gdt = 4; end


%% parameters for various GE gradient coils
if ischar(chronaxie)
    switch lower(chronaxie)
        case 'xrmw',  chronaxie=360d-6; rheobase=20.0; alpha=0.324;
        case 'xrm',   chronaxie=334d-6; rheobase=23.4; alpha=0.333;
        case 'whole', chronaxie=370d-6; rheobase=23.7; alpha=0.344;
        case 'zoom',  chronaxie=354d-6; rheobase=29.1; alpha=0.309;
        case 'hrmw',chronaxie=642.4d-6; rheobase=17.9; alpha=0.310;
        otherwise, error('gradient coil (%s) unkown',chronaxie);
    end
end


%% checks
if gmax<1, warning('gmax(=%g)<1',gmax); input(''); end
if fov<1,  warning('fov(=%g)<1',fov);   input(''); end
if gdt<4,  warning('gdt(=%g)<4',gdt); end
if ~isempty(pthresh)
    if pthresh<0, error('pthresh(=%g)<0',pthresh); end
end
dki = kdt/gdt;
if dki<0.5
    error('kdt: dki<0.5');
end
if (abs(dki-round(dki))>1d-10)
    if abs(dki-0.5)>1d-10
        error('kdt(=%g[us]) not multiple of gdt(=%g[us])',kdt,gdt);
    end
else
    dki = round(dki);
end
ng_off = nk_off*dki;


%% convert to standard SI units
fov = fov*1d-3;           % [mm] -> [m]
kdt = kdt*1d-6;           % [us] -> [s]
gdt = gdt*1d-6;           % [us] -> [s]
gmax = gmax*1d-3;         % [mT/m] -> [T/m]
res = fov/npix;           % [m] resolution


%% design spiral
rgamma = abs(gyrogamma('1h')/gyrogamma(nucleus));
fov_mns = fov/rgamma;
res_mns = res/rgamma;
fprintf('gamma ratio=%g\n',rgamma);
fprintf('scaled FOV (fov_mns) = %g [m]\n',fov_mns); 
fprintf('scaled res (res_mns) = %g [mm]\n',res_mns*1d3);
gmax_nyquist = 2*pi/(gyrogamma('1h')*kdt*fov_mns);
fprintf('gmax_nyquist = %g [mT/m]\n',gmax_nyquist*1d3);
if (gmax_nyquist<gmax)
    fprintf('Attention: approaching sampling BW limited regime\n');
    if (kdt>gdt)
        fprintf('!!! Undersampling will occur: reduce kdt !!!\n'); 
        input('press key to continue');
    end
end
pause(1);

figure(1);
SRmin = rheobase/alpha;
[~,g1,~,t1] = vds_pns(fufa_smax*smax*1e2,fufa_gmax*gmax*1e2,gdt, ...
    arms,[fov_mns*1e2,0],1/(2*res_mns*1e2),fufa_pthresh*pthresh,chronaxie,SRmin);


%% calculate single interleave
gspir = 1e-2*g1;                % convert to SI unit [T/m]
kspir = gdt*cumsum(gspir);      % [s*T/m]


%% rewinders
if balanced
    [grewx,ngx] = gradient_lobe(-real(kspir(1,end)),real(gspir(1,end)),...
        fufa_gmax_rew*gmax/sqrt(2),fufa_smax_rew*smax/sqrt(2),gdt);
    [grewy,ngy] = gradient_lobe(-imag(kspir(1,end)),imag(gspir(1,end)),...
        fufa_gmax_rew*gmax/sqrt(2),fufa_smax_rew*smax/sqrt(2),gdt);
    
    grew = zeros(1,max([ngx,ngy])+1);
    grew(1:ngx) = 1*grewx;
    grew(1:ngy) = grew(1:ngy)+1i*grewy;
else
    % nrew = ceil(gmax/smax/gdt/fufa_gmax_rew);
    nrew = ceil(abs(gspir(1,end))/smax/gdt/fufa_gmax)+1;
    grew = linspace(real(gspir(1,end)),0,nrew) + ...
        1i*linspace(imag(gspir(1,end)),0,nrew);
end


%% construct full (not yet rotated) gradient waveform
if isodd(ng_off+size(gspir,2)+size(grew,2))
    neven = 1;
else
    neven = 0;
end
gg = [zeros(1,ng_off),gspir,grew,zeros(1,neven)].';
ng = size(gg,1);
nk = ceil(length(kspir)/dki);


%% k-space trajectory (complex -> 2D)
if dki>1
    kc = kspir(1:dki:length(kspir))/2/pi*res_mns*gyrogamma('1H'); 
else
    kc = kspir/2/pi*res_mns*gyrogamma('1H');
end
if abs(dki-0.5)<1d-10
    kc = interp1((1:length(kc)),kc,(1:0.5:length(kc)+0.5),'spline');
end


%% rotate trajectories
phi = 0;                               % rotations [deg]
if arms>1, phi = mod(-(0:arms-1)/arms*360,360); end

nintlv = size(phi,2);                  % #interleaves
ng2 = size(phi,2);                     % #interleaves stored in grad
if do_rot_file, ng2 = 1; end
grad = zeros(ng,ng2,2);                % init grad

k = bsxfun(@times,kc,exp(-1i*phi*pi/180).').';
k = k(:).';
if ~do_rot_file
    gg = bsxfun(@times,gg,exp(-1i*phi*pi/180));
end
grad(:,:,1) = real(gg);
grad(:,:,2) = imag(gg);
if arms==1, phi = []; end

% calculate density compensation function (dcf)
% appr dcf via analytical radial expression to determine threshold
n2 = arms; mtx = npix;
tmp = pi/n2*(k(1,2:end).*conj(k(1,2:end))-k(1,1:end-1).*conj(k(1,1:end-1)))*mtx^2;
dcf_thresh = max(tmp);
dcf = voronoi_area(k*npix,'','',dcf_thresh);


%% calculate time and index list
if acq_round2n
    % round to 2^n
    acq_pts = 2^ceil(log2(ng/dki));
else
    % round to multiple of 32
    acq_pts = ceil((length(gspir)+ng_off)/dki/32)*32;
end
if acq_pts>16384
    warning('#sampling pts/exc (=%g)>16384: exceeds fidall limit',acq_pts);
end
inds = false(1,acq_pts);               % index for acq data (single arm)
inds(1,(nk_off+(1:nk))) = true;
ind = repmat(inds,[nintlv 1]);         % index acq all interleaves
ts = (0:nk-1)*kdt;                     % time of single trajectory acq [s]
t = repmat(ts,[nintlv 1]);             % time of all trajectories [s]


%% plot pns
if ~isempty(chronaxie)
    % figure(2);
    [~,~,ptmax,gmax_act,smax_act] = pns(grad,chronaxie,rheobase,alpha,gdt);
    ptmax = max(ptmax);
else
    gmax_act = max(max(abs(gg)));
    smax_act = max(max(abs(diff(gg,[],1))))/gdt;
end


%% print info about waveform
fprintf('actual gmax  = %g [mT/m]\n',gmax_act*1d3);
if (gmax_act>gmax), warning('gmax exceeded'); end
fprintf('actual smax  = %g [T/m/s]\n',smax_act);
if (smax_act>smax), warning('smax exceeded'); end

fprintf('Sequence details\n');
fprintf('Acq BW = %g [kHz] (full)\n',1d-3/kdt);
fprintf('gdt = %g [us]; kdt = %g [us]\n',gdt*1d6,kdt*1d6);
fprintf('g_pts = %gx%g; k_pts = %d; acq_pts/exc = %g\n',...
    size(gg,1),size(gg,2),size(k,2),acq_pts);
t_arm = t1(1,end); t_rew = (length(grew)+1)*gdt;
fprintf('t_arm = %g [ms]; t_rew = %g [ms]\n',t_arm*1d3,t_rew*1d3);
fprintf('t_seq = %g [ms]; nintlv = %d\n',ng*gdt*1d3,nintlv);
fprintf('actual PNS = %g[%%]\n',ptmax);
if ptmax>pthresh, warning('PTmax(=%g)>pthresh(=%g)',ptmax,pthresh); end


%% checks
if size(k,1)~=1,   error('size(k,1)(=%g)~=1',  size(k,1)); end
if size(dcf,1)~=1, error('size(dcf,1)(=%g)~=1',size(dcf,1)); end
if size(dcf,2)~=size(k,2)
    warning('size(dcf,2)~=size(k,2): interpolating dcf'); 
    dcf = interp1(linspace(0,1,size(dcf,2)),dcf,linspace(0,1,size(k,2)));
end
if size(k)~=sum(ind(:)),     error('size(k)~=sum(ind(:))'); end
if (prod(size(t))~=size(k,2))
    warning('prod(size(t))~=size(k,2)'); 
end

%% generate filename
if (islogical(fname) || isnumeric(fname))
    if fname
        if isnumeric(nucleus)
            warning('Please enter string for nucleus');
            nuc = num2str(nucleus);
        else
            nuc = nucleus;
        end
        balstr = ''; if balanced, balstr = '_blncd'; end
        fname = sprintf('spiral_pns%d_%s_fov%g_mtx%g_arms%g_kdt%g_gmax%g_smax%g_dur%g%s', ...
            round(ptmax),nuc,round(fov*1d3),npix,arms,kdt*1d6,...
            round(gmax_act*1d3),round(smax_act),...
            round(ng*gdt*1d4)/10,balstr);
        fname = strrep(fname,'.','p');
        if abs(gdt*1d6-4)>1d-10, fname = [fname '_gdt' num2str(gdt*1d6)]; end
    else
        fname = [];
    end
end


%% export waveforms
if ~isempty(fname)
    desc = sprintf('design_spiral_pns(%g,%g,%g,%g,fname,%g,%g,%s,%g,%g,%g,%g,%g)\n',...
        fov*1d3,npix,arms,kdt*1d6, gmax*1d3,smax,nuc,pthresh,...
        chronaxie,rheobase,alpha);
    
    fprintf('Saving gradient + mat file to\n\t%s\n',fname);
    fprintf('%s\n',desc);
    
    mtx = npix;
    out = write_ak_wav([fname '.wav'],grad,1/kdt,fov_mns,desc,acq_pts,...
        [],[],gdt*1d6);
    save(fname,'out','k','dcf','t','ind','fov','mtx',...
        'npix','gdt','kdt','gmax','smax','nucleus',...
        't_rew','t_arm','ng_off','gmax_nyquist','rgamma','fov_mns','phi');
else
    out.gmax = gmax;
    out.smax = smax;
    out.gdt = gdt;
    out.kdt = kdt;
    out.grad = gg;
    out.bw = 1/kdt;
    out.fov = fov;
    out.ptmax = ptmax;
end

end      % main function design_spiral_pns.m
