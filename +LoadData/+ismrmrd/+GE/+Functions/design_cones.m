function [k,dcf,t,grad,ind,bw] = design_cones(fov,mtx,tacq,fname,gmax,smax,...
    nucleus,sampfac,gdt,order,dcf_method,verb)
%DESIGN_CONES  Design isotropic 3D Cones sampling trajectory
% [k,dcf,t,grad,ind] = design_cones(fov,mtx,tacq,fname,gmax,smax,...
%                  nucleus,sampfac,order,balanced,dcf_method,verb)
%                                                        [unit]    (default)
%        fov  Field-of-view                              [mm]
%        mtx  Matrix size                                [1]
%       tacq  Length of acquisition window               [ms]      (15)
%      fname  Filename (string, or logical)                        (false)
%       gmax  Maximum gradient strength                  [mT/m]    (33)
%       smax  Maximum slewrate                           [T/m/s]   (120)
%    nucleus  Nucleus                                              ('1H')
%    sampfac  Sampling factor (>1 over-, <1 under-sampled)         (2 2 2)
%             1=frequency encoding direction (speed in k-space)
%             2=spacing between cones along theta
%             3=factor for #intlv/theta
%             Oversampling required for decent k-space distribution
%        gdt  Gradient update time                       [us]      (4)
%      order  Cones ordering: 1=subsequent,2=random                (1)
% dcf_method  0=ones(size(k)); 1=voronoi_area; 2+3=calc_sdc3       (4)
%             1+2->cap_dcf(expensive); 3+4->empirically reduced dcf
%       verb  Verbose mode                                         (1)
%
%          k  k-Space trajectory (normalised to +-0.5)
%             #acq-pts #dim
%        dcf  Density compensation function (Zwart/Pipe)
%          t  Time vector for f0 shift in recon     [s]
%       grad  Grad waveform (#pts/intlv,#intlv,#dim)     [T/m/s]
%        ind  Index for acquired data
%
%Literature:  Gurney PT et al. MRM 2006, 55:575-582 (cones design)
%             Riemer F et al.  MAGMA 2014, 27:35-46 (sodium)
%             Zwart NR et al. MRM 2012, 67:701-710  (density compensation)
%  Download:  http://mrsrl.stanford.edu/~ptgurney/cones.php
% 12/2018  Frank Riemer (f.riemer@web.de) (initial version)
%  8/2019  Rolf Schulte
if (nargin<1), help(mfilename); return; end


%% default parameters
timerVal = tic;              % record execution time
fufa_gmax = 0.99;            % factor to derate gmax
fufa_smax = [0.89 0.4];      % factor to derate smax_xy and smax_z
% obey smax limits for oblique: sqrt(sum(diff(grad,[],1).^2,3))/gdt
fufa_dcf = [1 0.1; 0.5 0.6]; % scale down dcf in steps (dcf_method==3)
fufa_dcf2 = 1;               % scale down dcf linearly (dcf_method==4)
if ~exist('mtx','var'),      mtx = []; end
if isempty(mtx),             mtx = 32; end
if ~exist('tacq','var'),     tacq = []; end
if isempty(tacq),            tacq = 15; end
if ~exist('fname','var'),    fname = []; end
if ~exist('gmax','var'),     gmax = []; end
if isempty(gmax),            gmax = 33; end
if ~exist('smax','var'),     smax = []; end
if isempty(smax),            smax = 120; end
if ~exist('nucleus','var'),  nucleus = ''; end
if isempty(nucleus),         nucleus = '1H'; end
if ~exist('sampfac','var'),  sampfac = []; end
if isempty(sampfac),         sampfac = [2 2 2]; end
if ~exist('gdt','var'),      gdt = []; end
if isempty(gdt),             gdt = 4; end
if ~exist('order','var'),    order = []; end
if isempty(order),           order = 1; end
if ~exist('dcf_method','var'), dcf_method = []; end
if isempty(dcf_method),      dcf_method = 4; end
if ~exist('verb','var'),     verb = []; end
if isempty(verb),            verb = 1; end


%% checks
if length(fov)~=1,     error('length(fov)~=1'); end
if fov<1,              warning('fov(=%g)<1',fov); end
if length(mtx)~=1,     error('length(mtx)~=1'); end
if gmax<1,             warning('gmax(=%g)<1',gmax); end
if gmax<1d-3,          error('gmax(=%g)<1d-3',gmax); end
sampfac = sampfac(:).';
if length(sampfac)~=3, error('length(sampfac)(=%g)~=3',length(sampfac)); end
if gdt<4,              error('gdt(=%g)<4[us]',gdt); end
if abs(rem(gdt,4))>1d-10, error('gdt(=%g) not multiple of 4[us]',gdt); end
gdt = gdt*1d-6;              % gradient update time [s]


%% misc parameters
nk_off = 2;                            % skip first acq points (artefact)
RFS_MAX_NSCANS = 16382;                % maximum number of scans
fov = fov*1d-3;                        % [mm] -> [m]
tacq = tacq*1d-3;                      % [ms] -> [s}
gmax = gmax*1d-3;                      % [mT/m] -> [T/m]
GAM = abs(gyrogamma(nucleus))/2/pi;    % gyromagnetic ratio [Hz/T]
res = fov/mtx;                         % resolution [m]
res_mns = res*GAM/(gyrogamma('1H')/2/pi); % resolution scaled to nucleus [m]
ng_init = ceil((tacq+1d-3)/gdt);       % #waveform points to initialise


%% calculate 3D Cones waveforms
if verb
    fprintf('Calculating 3D Cones waveforms via function findcone\n');
end
ncones = ceil(pi/4*mtx*sampfac(1,2));  % #individual cones for half k-space
if verb>0, fprintf('ncones=%g (half)\n',ncones); end
ng_cone = zeros(1,ncones);             % #waveform points per cone and ramp
% ng_max = 0;                          % longest waveform (cone only)
nint_arr = zeros(1,ncones);            % #intlvs per cone (not rounded up)
grad_arr = zeros(ng_init,3,ncones);    % cones waveforms not yet rotated
theta = ((ncones:-1:1)-0.5)/ncones*pi/2;    % polar angles [rad]

for lcon = 1:ncones
	% Calculation of specific waveform for this cone surface	
	[grd_cone,~,nint] = findcone(res_mns*1d3,fov*100,...
        round(tacq/gdt),theta(1,lcon)*[1 1],0.1,gdt,1.0,...
        smax*100*fufa_smax,gmax*100*fufa_gmax,1);
	
    nint_arr(1,lcon) = nint;           % #interleaves per cone
    % if size(grd_cone,1)>ng_max, ng_max = size(grd_cone,1); end
    
    % straight ramp-down of gradients
    ng_ramp = ceil(max(abs(grd_cone(end,:)))/smax/gdt*sqrt(3)/100);
    grd_rew = linspace(1,0,ng_ramp).'*grd_cone(end,:);
    
    grd = [grd_cone ; grd_rew];        % assemble gradient
	ng_cone(1,lcon) = size(grd,1);
    if ng_cone(1,lcon)>ng_init
        warning('ng_lcone(1,%g)(=%g)>ng_init(=%g)',...
            lcon,ng_cone(1,lcon),ng_init); 
    end
    grad_arr((1:ng_cone(1,lcon)),:,lcon) = grd;
end
ng1 = ceil((max(ng_cone)+1)/2)*2;      % actually used grad waveform points
grad_arr = grad_arr((1:ng1),:,:);      % truncate to actually used size
if verb
    fprintf('nint=%.4g+-%.4g\n',mean(nint_arr),std(nint_arr)); 
end


%% rotate waveform interleaves and assemble final grad waveform
if verb, fprintf('Rotating and assembling grad waveforms\n'); end
nexc = sum(ceil(sampfac(1,3)*nint_arr),2);
if verb>0, fprintf('nexc=%g (half)\n',nexc); end
if nexc*2>RFS_MAX_NSCANS
    error('nexc(=%g)*2>RFS_MAX_NSCANS(=%g)',nexc,RFS_MAX_NSCANS); 
end
grad = zeros(ng1,nexc,3);              % init final grad waveform
nint_arr_exc = zeros(1,nexc);
% golden_angle = 180/((sqrt(5)+1)/2);
tau = (1+sqrt(5))/2;
golden_angle = 180/(tau+2);            % tiny golden angle

lexc = 0;
for lcon2 = 1:2*ncones
    if lcon2<=ncones
        lcon = lcon2;
        sgn = 1;
    else
        lcon = 2*ncones-lcon2+1;
        sgn = -1;
    end
    % fprintf('lcon2=%d; lcon=%d\n',lcon2,lcon);
    nint = ceil(sampfac(1,3)*nint_arr(1,lcon));    % #rotations per cone
    dphi = 360/nint;                   % phase steps [deg]
    phi_off = mod(lcon2*golden_angle,360);
    % phi_off = mod(lcon*golden_angle,360);
    % phi_off = mod(randn(1,1)*1000,360);
    phi = phi_off+(0:(nint-1))*dphi;
    if verb>1
        fprintf('lcon=%g; nint=%g; dphi=%g; phi_off=%g\n',...
            lcon,nint,dphi,phi_off);
    end
    ii = lexc+(1:nint);
    grad(:,ii,1) =  grad_arr(:,1,lcon)*cosd(phi) + grad_arr(:,2,lcon)*sind(phi);
    grad(:,ii,2) = -grad_arr(:,1,lcon)*sind(phi) + grad_arr(:,2,lcon)*cosd(phi);
    grad(:,ii,3) = sgn*repmat(grad_arr(:,3,lcon),[1 nint]);
    
    nint_arr_exc(1,ii) = nint_arr(1,lcon);
    
    lexc = lexc+nint;
end
grad = grad/100;                       % [G/cm] -> [T/m]
% grad = [grad , grad(:,(nexc:-1:1),:)]; % mirror grad along xy-plane
% ii = nexc+(1:nexc);
% grad(:,ii,3) = -grad(:,ii,3);

nexc = 2*nexc;
% nint_arr_exc = [nint_arr_exc , nint_arr_exc(1,end:-1:1)];

if order==2
    if verb>0, fprintf('Randomising views\n'); end
    ii = randperm(nexc);
    grad = grad(:,ii,:);
end


%% sampling parameters
gmax_act = max(max(sqrt(sum(grad.^2,3))));
bw = gmax_act*fov*GAM*sampfac(1,1);
kdt = floor(1d6/bw/2)*2d-6;            % round-down to multiple of 2[us]
if kdt<2d-6
    warning('kdt(=%g)<2[us]; setting kdt=2[us]',kdt*1d6);
    kdt = 2d-6;
end
bw = 1/kdt;                            % full (+-) acq bandwidth [Hz]
if bw<5d3, warning('bw(=%g)<5000[Hz]',bw); end
ng_off = nk_off*kdt/gdt;               % #grad pts offset for start of traj
if abs(ceil(ng_off)-ng_off)>1d-10
    warning('ng_off(=%g) not even; rounding',ng_off);
end
ng_off = round(ng_off);
if verb>0
    fprintf('gdt=%g[us]; kdt=%g[us]; ',gdt*1d6,kdt*1d6);
    fprintf('nk_off=%g; ng_off=%g\n',nk_off,ng_off);
    fprintf('size(grad)=(%g,%g,%g)\n', size(grad));
end
grad = [zeros(ng_off,nexc,3) ; grad];  % adding offset
if isodd(size(grad,1))
    grad = [grad ; zeros(1,nexc,3)];   % even grad length
end


%% caculate k-space trajectory
if verb, fprintf('Calculating k-space trajectory\n'); end
kk = cumsum(grad(((ng_off+1):ng1),:,:),1)*gdt*GAM*res; % normalised to +/-0.5
nkg = size(kk,1);                      % #k-pts relative to gdt
nkk = floor(gdt/kdt*nkg);              % #k-pts relative to kdt
% resampling kk from gdt to kdt
kk = interp1((0:(nkg-1)).'*gdt,kk,(0:(nkk-1)).'*kdt,'spline');

k = zeros(1,nkk*nexc,3); 
t = zeros(nkk*nexc,1);
n_ind = ceil((tacq/kdt+nk_off+2)/32)*32;
ind = false(nexc,n_ind);

lwarn1 = 0; lwarn2 = 0;                % collect warnings
lk = 0;                                % loop count interleaves
for lexc=1:nexc                        % loop through excitations
    abs_kc = sqrt(sum(kk(:,lexc,:).^2,3));  % k-space radius
    % truncate k-space at radius=0.5
    if max(abs_kc)<0.5
        lwarn1 = lwarn1+1;
        ii = nkk;
    else
        if max(abs_kc)>0.7, lwarn2 = lwarn2+1; end
        ii = find(abs_kc>0.5);
        ii = ii(1);
    end
    % concatenate excitations
    k(1,(lk+(1:ii)),:) = reshape(kk(1:ii,lexc,:),[1 ii 3]);
    t((lk+(1:ii)),1) = (0:(ii-1)).'*kdt;    % time [s]
    ind(lexc,(nk_off+(1:ii))) = true;       % index for data
    
    lk = lk+ii;
end
if lwarn1>0, warning('max(abs_kc)<0.5: %g times',lwarn1); end
if lwarn2>0, warning('max(abs_kc)>0.7: %g times',lwarn2); end
if any(any(abs(k(1,(lk+1):end,:))))>0
    warning('k(4truncating) not zero');
end
if sum(ind(:))~=lk, warning('sum(ind(:))(=%g)~=lk(=%g)',sum(ind(:)),lk); end
k = k(1,(1:lk),:);                     % truncating to actually used size
t = t((1:lk),1);                       % truncating to actually used size
if verb, fprintf('size(k)=(%g,%g,%g)\n', size(k)); end


%% calculating density-compensation function (dcf)
switch dcf_method
    case 1
        if verb, fprintf('Calculating dcf via voronoi_area\n'); end
        dcf = voronoi_area(k*mtx);
    case 2
        if verb, fprintf('Calculating dcf via calc_sdc3\n'); end
        dcf = calc_sdc3(k,mtx,15);
    case 3
        if verb, fprintf('Calculating dcf via calc_sdc3\n'); end
        dcf = calc_sdc3(k,mtx,15);
        krad = sqrt(sum(k.^2,3));
        dk = gmax*kdt*GAM*res;
        for l=1:size(fufa_dcf,2)
            k_thresh = 0.5-fufa_dcf(1,l)*dk;
            if verb
                fprintf('\tscaling dcf: k krad>%g->*%g\n',...
                    k_thresh,fufa_dcf(2,l));
            end
            ii = krad>k_thresh;
            dcf(1,ii) = dcf(1,ii)*fufa_dcf(2,l);
        end
    case 4
        if verb, fprintf('Calculating dcf via calc_sdc3\n'); end
        dcf = calc_sdc3(k,mtx,8);
        krad = sqrt(sum(k.^2,3));
        dk = gmax*kdt*GAM*res;
        k_thresh = 0.5-fufa_dcf2*dk;
        ii = krad>k_thresh;
        att = (0.5-krad(1,ii))/(fufa_dcf2*dk);
        att(att<0) = 0;
        dcf(1,ii) = dcf(1,ii).*att;
        
    otherwise
        if verb, fprintf('Skipping dcf -> set to one\n'); end
        dcf = ones(size(k));
end
if ((dcf_method>0)&&(dcf_method<3))
    if verb, fprintf('Capping dcf\n'); end
    dcf = cap_dcf(dcf,k,mtx);
end


%% print info
ss = sqrt(sum(diff(grad,[],1).^2,3))/gdt;   % slewrate worst case (rotated)
smax_act = max(ss(:));                      % max actual slew rate [T/m/s]
if verb
    fprintf('gmax: actual(=%.3g,%.3g,%.3g->%.4g); nominal(=%g) [mT/m]\n',...
        squeeze(max(max(abs(grad),[],1),[],2))*1d3,gmax_act*1d3,gmax*1d3);
    fprintf('smax: actual(=%.4g,%.4g,%.4g->%.5g); nominal(=%g) [T/m/s]\n',...
        squeeze(max(max(abs(diff(grad,[],1)/gdt),[],1),[],2)),...
        smax_act,smax);
    fprintf('acq: BW=%g[kHz]; #pts=%g; duration=%g[ms]\n',...
        1d-3/kdt,n_ind,n_ind*kdt*1d3);
end
if gmax_act>(gmax+1d-10)
    warning('gmax: actual(=%g)>nominal(=%g) [mT/m]',gmax_act,gmax); 
end
if smax_act>smax
    warning('smax: actual(=%g)>nominal(=%g) [T/m/s]',smax_act,smax); 
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
        gdt_str = '';
        if gdt>4d-6, gdt_str = sprintf('_gdt%d',round(gdt*1d6)); end
        fname = sprintf('cone_%s_fov%g_mtx%g_intlv%d_kdt%g_gmax%g_smax%g_dur%g%s', ...
            nuc,round(fov*1d3),mtx,nexc,kdt*1d6,round(gmax_act*1d3),...
            round(smax_act),round(tacq*1d4)/10,gdt_str);
        fname = strrep(fname,'.','p');
    else
        fname = [];
    end
end


%% Saving waveform
if ~isempty(fname)
    fov1H = fov/abs(gyrogamma('1h'))*abs(gyrogamma(nucleus));
    desc = sprintf('design_cones(%g,%g,%g,%g,fname,%g,%g,%s,[%g,%g,%g],%g)\n',...
        fov*1d3,mtx,kdt*1d6,tacq*1d3,gmax*1d3,smax,nuc,sampfac,verb);
    npix = mtx;
    fprintf('Saving gradient + mat file to\n\t%s\n',fname);
    fprintf('%s\n',desc);
    out = write_ak_wav([fname,'.wav'],grad,bw,fov1H,desc,n_ind,'','',gdt*1d6);
    dcf = single(dcf); k = single(k); t = single(t);
    if isfield(out,'wave'), out = rmfield(out,'wave'); end
    if isfield(out,'grad'), out = rmfield(out,'grad'); end
    save(fname,'out','k','t','ind','fov','fov1H','mtx','gmax','smax',...
        'nucleus','dcf','n_ind','npix','nexc','gmax_act','smax_act');
end

fprintf('design_cones.m finished: ');
fprintf('total execution  time = %.2f [s]\n',toc(timerVal));


end      % main function design_cones.m
