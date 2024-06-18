function [k,dcf,t,grad,ind,out,fname]=design_spiral3d(fov,npix,arms,kdt,...
    fname,gmax,smax,nucleus,mode3d,acq_round2n,do_rot_file,rewinder,...
    sampfac,do_inout,gdt,save_reduced)
%DESIGN_SPIRAL3D Design a spiral with delayed acq for fast CSI
%[k,dcf,t,grad,ind,out,fname]=design_spiral3d(fov,npix,arms,kdt,...
%    fname,gmax,smax,nucleus,mode3d,acq_round2n,do_rot_file,rewinder,...
%    sampfac,do_inout,gdt,save_reduced)
%
%        fov  Field of view                                [mm]
%       npix  #pixels (Cartesian resolution after gridding)
%       arms  #spatial interleaves                                      (1)
%             mode3d>0 -> #spokes=arms(2)
%        kdt  k-space sampling time ([]->automatic)        [us]        ([])
%      fname  Filename: if given->write_ak_wav 
%             if logical true: generate name                        (false)
%       gmax  Max gradient amp                             [mT/m]      (33)   
%       smax  Max slew rate                                [T/m/s]    (120)  
%    nucleus  Nucleus                                                ('1H')
%     mode3d  3D encoding (0=2D; >1->3D rotated spirals)                (0)
%             1=consecutive, 2=golden angle along radial
%acq_round2n  Round up #acq pts to to 2^n                           (false)
%do_rot_file  Write single static spiral into gradient               (true)
%             waveform file & rotate via vap_phiXX.fdl
%   rewinder  Gradient ramp down method:                                (0)
%             0=fast, 1=return to centre of k-space, 2=crushing
%             2->[2 area-factor(default=1)]
%    sampfac  Sampling factor                                   ([1,1,0,1])
%             Directions: readout,arm,vds,radial: >1->over,<1->under sample
%             vds=variable density spiral:0=normal,-1..0=undersample high k
%   do_inout  In-out spiral (2D only)                               (false)
%        gdt  Gradient update time (must be multiple of 4) [us]         (4)
%save_reduced Reduce mat file size                                      (2)
%
%          k  k-space trajectory  [-0.5..0.5]
%        dcf  Density compensation function (calculated with vornoi_area)
%          t  Time  (nexc x recon-pts)
%       grad  Gradient waveform [T/m] with dt=4us
%             (#pts/intlv, #intlv, #dim)
%        ind  Index (2nd dim) for k-space points on spiral (excl ramps,etc)
%        out  Output structure of wrt_wavs
%      fname  Filename; if nargout==1 -> return fname
%
% 2/2021 Rolf Schulte
if (nargin<1), help(mfilename); return; end


%% fixed parameters
fufa_gmax = 0.99;
fufa_smax = 0.99;
nk_off = 2;                  % skip first acq points (artefact)


%% input parameters
if abs(rem(npix,1))>1d-10
    warning('npix(=%g) not natural number; rounding to %d',npix,round(npix));
    npix = round(npix);
end
if ~exist('arms','var'),        arms = []; end
if isempty(arms),               arms = 1; end
if ~exist('kdt','var'),         kdt = []; end
if ~exist('fname','var'),       fname = []; end
if ~exist('gmax','var'),        gmax = []; end
if isempty(gmax),               gmax = 33; end
if ~exist('smax','var'),        smax = []; end
if isempty(smax),               smax = 120; end
if ~exist('nucleus','var'),     nucleus = []; end
if isempty(nucleus),            nucleus = '1H'; end
if ~exist('mode3d','var'),      mode3d = []; end
if isempty(mode3d),             mode3d = 0; end
if length(mode3d)<2,            mode3d = [mode3d 1]; end
if ~exist('acq_round2n','var'), acq_round2n = []; end
if isempty(acq_round2n),        acq_round2n = false; end
if ~exist('do_rot_file','var'), do_rot_file = []; end
if isempty(do_rot_file),        do_rot_file = true; end
if ~exist('rewinder','var'),    rewinder = []; end
if isempty(rewinder),           rewinder = 0; end
if length(rewinder)<2,          rewinder = [rewinder 1]; end
if ~exist('sampfac','var'),     sampfac = []; end
if isempty(sampfac),            sampfac = 1; end
if length(sampfac)<2,           sampfac = [sampfac 1]; end
if length(sampfac)<3,           sampfac = [sampfac 0]; end
if length(sampfac)<4,           sampfac = [sampfac 1]; end
if ~exist('do_inout','var'),    do_inout = []; end
if isempty(do_inout),           do_inout = false; end
if ~exist('gdt','var'),         gdt = []; end
if isempty(gdt),                gdt = 4; end
if gdt<4, error('gdt(=%g)<4',gdt); end
if ~exist('save_reduced','var'),save_reduced = []; end
if isempty(save_reduced),       save_reduced = 2; end


%% checks
if gmax<1, warning('gmax(=%g)<1',gmax); input(''); end
if fov<1,  warning('fov(=%g)<1',fov);   input(''); end
if isempty(kdt)
    kdt = 2d12*pi/(gmax*abs(gyrogamma(nucleus))*fov);
    kdt = floor(kdt/4/sampfac(1))*4;
    if kdt<4, kdt = 2; end
end
if kdt<2, warning('kdt(=%g)<2',kdt); input(''); end
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
fov = fov*1d-3;              % [mm] -> [m]
kdt = kdt*1d-6;              % [us] -> [s]
gdt = gdt*1d-6;              % [us] -> [s]
gmax = gmax*1d-3;            % [mT/m] -> [T/m]
res = fov/npix;              % [m] resolution


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
    pause(1);
    if (kdt>gdt)
        fprintf('!!! Undersampling will occur: reduce kdt !!!\n'); 
        input('press key to continue');
    end
end
if do_inout, sampfac(2) = sampfac(2)/2; end
figure(100); clf;
[~,g1,~,t1] = vds(fufa_smax*smax*1e2,fufa_gmax*gmax*1e2,gdt,...
    arms(1)/sampfac(2),fov_mns*1e2*[1,sampfac(3)],1/(2*res_mns*1e2));
gspir = 1e-2*g1;             % convert to SI unit [T/m]


%% in-out spiral
if do_inout
    ga_pre = -gdt*sum(gspir);
    [gprex,ngx] = gradient_lobe(real(ga_pre),[0 real(gspir(1,end))],...
        fufa_gmax*gmax/sqrt(2),fufa_smax*smax/sqrt(2),gdt);
    [gprey,ngy] = gradient_lobe(imag(ga_pre),[0 imag(gspir(1,end))],...
        fufa_gmax*gmax/sqrt(2),fufa_smax*smax/sqrt(2),gdt);
        
    ng_pre = max([ngx,ngy])+1;
    if dki>1, ng_pre = ceil(ng_pre/dki)*dki; end
    gpre = [zeros(1,ng_pre-ngx) gprex] + 1i*[zeros(1,ng_pre-ngy) gprey];
    gspir = [gpre gspir(1,end:-1:1),gspir];
    nk_off = round(nk_off+ng_pre/dki);
end


%% calculate k-space
kspir = gdt*cumsum(gspir);   % [s*T/m]
if do_inout, kspir = kspir(1,(ng_pre+1):end); end


%% rewinders
switch rewinder(1)
    case 0         % fast grad ramp down
        nrew = ceil(abs(gspir(1,end))/smax/gdt/fufa_gmax)+1;
        grew = linspace(real(gspir(1,end)),0,nrew) + ...
            1i*linspace(imag(gspir(1,end)),0,nrew);
        balstr = '';
    case 1         % return to centre of k-space
        [grewx,ngx] = gradient_lobe(-real(kspir(1,end)),real(gspir(1,end)),...
            fufa_gmax*gmax/sqrt(2),fufa_smax*smax/sqrt(2),gdt);
        [grewy,ngy] = gradient_lobe(-imag(kspir(1,end)),imag(gspir(1,end)),...
            fufa_gmax*gmax/sqrt(2),fufa_smax*smax/sqrt(2),gdt);
        
        ng = max([ngx,ngy])+1;
        grew = [grewx zeros(1,ng-ngx)] + 1i*[grewy zeros(1,ng-ngy)];
        balstr = '_blncd';
    case 2         % add crusher
        ga_rew = rewinder(2)*abs(kspir(1,end))*gspir(1,end)/abs(gspir(1,end));
        [grewx,ngx] = gradient_lobe(real(ga_rew),real(gspir(1,end)),...
            abs(real(gspir(1,end))),fufa_smax*smax/sqrt(2),gdt);
        [grewy,ngy] = gradient_lobe(imag(ga_rew),imag(gspir(1,end)),...
            abs(imag(gspir(1,end))),fufa_smax*smax/sqrt(2),gdt);
        
        ng = max([ngx,ngy])+1;
        grew = [grewx zeros(1,ng-ngx)] + 1i*[grewy zeros(1,ng-ngy)];
        balstr = '_crush';
    case 3         % return partially to centre of k-space
        ga_rew = -rewinder(2)*kspir(1,end);
        [grewx,ngx] = gradient_lobe(real(ga_rew),real(gspir(1,end)),...
            fufa_gmax*gmax/sqrt(2),fufa_smax*smax/sqrt(2),gdt);
        [grewy,ngy] = gradient_lobe(imag(ga_rew),imag(gspir(1,end)),...
            fufa_gmax*gmax/sqrt(2),fufa_smax*smax/sqrt(2),gdt);
        
        ng = max([ngx,ngy])+1;
        grew = [grewx zeros(1,ng-ngx)] + 1i*[grewy zeros(1,ng-ngy)];
        balstr = sprintf('_ret%g',rewinder(2));
    otherwise
        error('rewinder(1)(=%g)~=0,1 or 2',rewinder(1));
end


%% construct full (not yet rotated) gradient waveform
gg = [zeros(1,ng_off),gspir,grew].';
if mod(size(gg,1),2)         % adding 0 if odd
    gg = [gg ; 0];
end
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
phi = 0; theta = [];         % rotations [deg]
if arms(1)>1
    phi = mod(-(0:arms(1)-1)/arms(1)*360/(do_inout+1),360);
end
if mode3d(1)>0
    % 3D rotation in fidall:
    % rotation of theta around y then of phi around z
    % -> put spirals into xz-plane -> swap phi & theta
    kmax = 1/(2*res_mns)/2;
    if length(arms)>1
        Ns = arms(2);                            % #radial spokes
    else
        Ns = ceil(pi*kmax*fov*sampfac(4))*2;     % #radial spokes
    end
    fprintf('#spokes in 3D = %g\n',Ns);
    switch mode3d(1)
        case 1          % consecutive
            theta = repmat(phi,[1 Ns]);
            phi = (0:Ns-1)/Ns*180;
            phi = repmat(phi,[arms(1) 1]);
            phi = phi(:).';
        case {2,3}      % golden angle along radial
            theta = zeros(1,arms(1)*Ns);
            for ll=1:arms(1)
                ii1 = (1:Ns)+(ll-1)*Ns;
                theta(1,ii1) = mod((-(0:Ns-1)/Ns-(ll-1)/arms(1))*360,360);
            end
            ga = golden_angle(mode3d(2),true);
            phi = mod((0:Ns-1)*ga,360);     % MRP uses 180!?
            phi = repmat(phi,[1 arms(1)]);
            if mode3d(1)==3
                ii = randperm(Ns*arms(1));
                phi = phi(1,ii);
                theta = theta(1,ii);
            end
        case 4          % 2D golden angle
            % Fyhrdahl, MAGMA 2020
            ga1 = sub_fibonacci_series(1,100)/...
                sub_fibonacci_series(mode3d(2),102);
            ga2 = sub_fibonacci_series(1,100)/...
                sub_fibonacci_series(mode3d(2),101);
            fprintf('2D Golden angle-%d: ga1=%g; ga2=%g\n',mode3d(2),ga1,ga2);
            theta = acosd(mod(ga1*(0:(Ns*arms(1)-1)),2)-1);
            phi = 360*mod(ga2*(0:(Ns*arms(1)-1)),1);
            
            
        otherwise
            error('mode3d(=%g)~=1 or 2',mode3d(1));
    end
end
nintlv = size(phi,2);                  % #interleaves
ng2 = size(phi,2);                     % #interleaves stored in grad
if do_rot_file, ng2 = 1; end
grad = zeros(ng,ng2,2+(mode3d(1)>0));  % init grad

if (mode3d(1)==0)  % regular 2D spirals
    k = bsxfun(@times,kc,exp(-1i*phi*pi/180).').';
    k = k(:).';
    if ~do_rot_file
        gg = bsxfun(@times,gg,exp(-1i*phi*pi/180));
    end
    grad(:,:,1) = real(gg);
    grad(:,:,2) = imag(gg);
    if arms(1)==1, phi = []; end
    
    % calculate density compensation function (dcf)
    % appr dcf via analytical radial expression to determine threshold
    tmp = pi/arms(1)*(k(1,2:end).*conj(k(1,2:end))-k(1,1:end-1).*conj(k(1,1:end-1)))*npix^2;
    dcf_thresh = max(tmp);
    dcf = voronoi_area(k*npix,'','',dcf_thresh);
    % dcf = voronoi_area(k*npix);    % density compensation function
    k_pts = size(k,2);
else               % rotated 3D spirals
    if do_rot_file
        grad(:,:,1) = real(gg);
        grad(:,:,3) = imag(gg);
    else
        grad(:,:,1) =  real(gg)*(cosd(theta).*cosd(phi)) - ...
            imag(gg)*(sind(theta).*cosd(phi));
        grad(:,:,2) =  -real(gg)*(cosd(theta).*sind(phi)) + ...
            imag(gg)*(sind(theta).*sind(phi));
        grad(:,:,3) = real(gg)*sind(theta) + ...
            imag(gg)*cosd(theta);
    end
    
    kc = kc.';
    
    if ((save_reduced==0) || (nargout>1))
        k = zeros(length(kc),nintlv,3);
        k(:,:,1) =  real(kc)*(cosd(theta).*cosd(phi)) - ...
            imag(kc)*(sind(theta).*cosd(phi));
        k(:,:,2) =  -real(kc)*(cosd(theta).*sind(phi)) + ...
            imag(kc)*(sind(theta).*sind(phi));
        k(:,:,3) = real(kc)*sind(theta) + ...
            imag(kc)*cosd(theta);
        k = reshape(k,[1 length(kc)*size(phi,2) 3]);
    end
    
    % calculate density compensation function (dcf)
    % appr dcf via analytical radial expression to determine threshold
    if ((save_reduced<2) || (nargout>1))
        if mode3d(1)==1
            fprintf('Calculating dcf analytically\n');
            k2d = bsxfun(@times,kc.',exp(-1i*(0:arms(1)-1)/arms(1)*2*pi).').';
            k2d = k2d(:).';
            tmp = pi/arms(1)*(k2d(1,2:end).*conj(k2d(1,2:end))-k2d(1,1:end-1).*conj(k2d(1,1:end-1)))*npix^2;
            dcf_thresh = max(tmp);
            dcf2d = single(voronoi_area(k2d*npix,'','',dcf_thresh));
            dcf = 2*repmat(dcf2d.*abs(real(k2d)),[1 Ns]);
        else
            fprintf('Calculating dcf numerically via calc_sdc3\n');
            dcf = calc_sdc3(k,npix,5);
        end
    end
    
    k_pts = length(kc)*size(phi,2);
end


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
ts = (0:nk-1)*kdt;                     % time of single trajectory acq [s]
if ~((mode3d(1)>0) && (save_reduced>0))
    ind = repmat(inds,[nintlv 1]);         % index acq all interleaves
    t = repmat(ts,[nintlv 1]);             % time of all trajectories [s]
end


%% print info about waveform
gmax_act = max(max(abs(gg)));
smax_act = max(max(abs(diff(gg,[],1))))/gdt;
fprintf('actual gmax  = %g [mT/m]\n',gmax_act*1d3);
if (gmax_act>gmax), warning('gmax exceeded'); end
fprintf('actual smax  = %g [T/m/s]\n',smax_act);
if (smax_act>smax), warning('smax exceeded'); end

fprintf('Sequence details\n');
fprintf('Acq BW = %g [kHz] (full)\n',1d-3/kdt);
fprintf('gdt = %g [us]; kdt = %g [us]\n',gdt*1d6,kdt*1d6);
fprintf('g_pts = %gx%g; k_pts = %d; acq_pts/exc = %g\n',...
    size(gg,1),size(gg,2),k_pts,acq_pts);
t_arm = t1(1,end); t_rew = (length(grew)+1)*gdt;
fprintf('t_arm = %g [ms]; t_rew = %g [ms]\n',t_arm*1d3,t_rew*1d3);
fprintf('t_seq = %g [ms]; nintlv = %d\n',ng*gdt*1d3,nintlv);
if do_inout, fprintf('TE += %g [ms]\n',(t_arm+nk_off*kdt)*1d3); end


%% checks
if ~((mode3d(1)>0) && (save_reduced>0))
    if size(k,1)~=1, error('size(k,1)(=%g)~=1',size(k,1)); end
    if k_pts~=sum(ind(:))
        warning('k_pts(=%d)~=sum(ind(:))(=%d)',k_pts,sum(ind(:)));
    end
    if prod(size(t))~=size(k,2)
        warning('prod(size(t))~=size(k,2)');
    end
end
if exist('dcf','var')
    if size(dcf,1)~=1, error('size(dcf,1)(=%g)~=1',size(dcf,1)); end
    if size(dcf,2)~=k_pts
        warning('size(dcf,2)~=size(k,2): interpolating dcf');
        dcf = interp1(linspace(0,1,size(dcf,2)),dcf,linspace(0,1,size(k,2)));
    end
end


%% generate filename
if isnumeric(nucleus)
    warning('Please enter string for nucleus');
    nuc = num2str(nucleus);
else
    nuc = nucleus;
end
if (islogical(fname) || isnumeric(fname))
    if fname
        switch mode3d(1)
            case 0, str3d = '';
            case 1, str3d = sprintf('3D_intlv%g',nintlv); 
            case 2, str3d = sprintf('3Dga%d_intlv%g',mode3d(2),nintlv); 
            case 3, str3d = sprintf('3Dga%dR_intlv%g',mode3d(2),nintlv);
            case 4, str3d = sprintf('3Dga2_intlv%g',nintlv);
        end
        if do_inout, str3d = [str3d '_inout']; end
        
        fname = sprintf('spiral%s_%s_fov%g_mtx%g_arms%g_kdt%g_gmax%g_smax%g_dur%g%s', ...
            str3d,nuc,round(fov*1d3),npix,arms(1),kdt*1d6,...
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
    desc = sprintf('design_spiral3D(%g,%g,%g,%g,fname,%g,%g,%s,%g,%g',...
        fov*1d3,npix,arms(1),kdt*1d6, gmax*1d3,smax,nuc,mode3d(1),acq_round2n);
    desc = sprintf('%s,%g,[%g,%g],[%g,%g,%g,%g],%g,%g,%g)\n',desc,...
        do_rot_file,rewinder(1:2),sampfac(1:4),do_inout,gdt*1d6,save_reduced);
    fprintf('Saving gradient + mat file to\n\t%s\n',fname);
    fprintf('%s\n',desc);
    
    if exist([fname '.wav'],'file')
        warning('file (=''%s'') exists; overwriting',[fname '.wav']);
    end
    
    mtx = npix;
    out = write_ak_wav([fname '.wav'],grad,1/kdt,fov_mns,desc,acq_pts,...
        [],[],gdt*1d6);
    if isfield(out,'wave'), out = rmfield(out,'wave'); end
    if isfield(out,'grad'), out = rmfield(out,'grad'); end
    if ((mode3d(1)>0) && (save_reduced>0))
        nviews = nintlv;
        ks = zeros(1,nk,3);
        ks(1,:,1) = real(kc(1:nk,1)).';
        ks(1,:,3) = imag(kc(1:nk,1)).';
        if save_reduced>1
            save([fname '.mat'],'out','ks','ts','inds','nviews','fov','mtx',...
                'npix','gdt','kdt','gmax','smax','nucleus','t_rew','t_arm',...
                'ng_off','gmax_nyquist','rgamma','fov_mns','phi','theta');
        else
            dcf_range = max(dcf);
            dcf = single(dcf);
            save([fname '.mat'],'out','ks','ts','inds','nviews','fov','mtx',...
                'npix','gdt','kdt','gmax','smax','nucleus','t_rew','t_arm',...
                'ng_off','gmax_nyquist','rgamma','fov_mns','phi','theta',...
                'dcf','dcf_range');
        end
    else
        if (mode3d(1)>0) || (save_reduced>0)
            k = single(k);
            dcf = single(dcf);
            t = single(t);
        end
        save([fname '.mat'],'out','k','dcf','t','ind','fov','mtx',...
            'npix','gdt','kdt','gmax','smax','nucleus',...
            't_rew','t_arm','ng_off','gmax_nyquist','rgamma','fov_mns',...
            'phi','theta');
    end
    if do_rot_file
        if ~isempty(phi)
            write_fdl([fname '_phi.fdl'],phi,'phi'); 
        end
        if ~isempty(theta)
            write_fdl([fname '_theta.fdl'],theta,'theta');
        end
    end
else
    out.gmax = gmax;
    out.smax = smax;
    out.gdt = gdt;
    out.kdt = kdt;
    out.grad = gg;
    out.bw = 1/kdt;
    out.fov = fov;
end


%% plotting
dim = (2+(mode3d(1)>0));
figure(101); clf
yl = {'Gx [mT/m]','Gy [mT/m]','Gz [mT/m]'};
for l=1:dim
    subplot(dim,1,l);
    plot((1:size(grad,1))*gdt*1d3,grad(:,:,l)*1d3);
    xlabel('time [ms]');
    ylabel(yl{l})
    grid on
end


%% output arguments
if nargout==1
    k = fname;
end
if nargout>1
    if ~exist('dcf','var'), dcf = []; end
    if ~exist('t','var'),   t = ts; end
    if ~exist('ind','var'), ind = inds; end
end


end      % design_spiral3d.m


%% 3D golden means - tiny golden angle
function h3 = sub_fibonacci_series(N,n)
% Literature: Fyhrdahl, MAGMA 2020
% https://doi.org/10.1007/s10334-020-00859-z
% Eq. 5
h = [0,1,N];

for l=1:n
    h = [h(2),h(3),h(3)+h(1)];
end

h3 = h(3);

end

