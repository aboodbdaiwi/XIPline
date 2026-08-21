function [k,dcf,t,grad,ind,out,fname]=design_spiral_2d4d(fov,mtx,arms,kdt,...
    fname,gmax,smax,nucleus,dim,typ,sampfac,rewinder,order,do_rot_file,...
    redundant, man_corr_dcf, dda, num_echoes, rew_scale, grad_coil)

%DESIGN_SPIRAL2D4D Design a spiral with delayed acq for fast CSI
%[k,dcf,t,grad,ind,out,fname]=design_spiral_2d4d(fov,mtx,arms,kdt,...
%    fname,gmax,smax,nucleus,dim,typ,sampfac,rewinder,order,do_rot_file)
%        fov  Field-of-view                           [mm]
%        mtx  Matrix size                             [1]
%        kdt  Sampling dwell time (1/bw)              [us]      (16)
%      fname  Filename (string, or logical)                     (false)
%       gmax  Maximum gradient strength               [mT/m]    (33)
%       smax  Maximum slewrate                        [T/m/s]   (120)
%    nucleus  Nucleus                                           ('1H')
%        dim  Spatial dimensions (2=2D,4=2D+z-PE)          (3)
%        typ  Spiral variant:
%             0=out spiral, 1=in-out spiral (2D only)%          (0)
%             !!! TO BE IMPLEMENTED !!!
%    sampfac  Sampling factor (>1 over-, <1 under-sampled)      ([1,1,0,1])
%             Directions: readout,arm,vds,radial
%             vds=variable density spiral:0=normal,<0 undersample high k
%   rewinder  Gradient ramp down method:                        (0)
%             0=fast, 1=return to centre of k-space, 2=crushing
%             2->[2 area-factor(default=1)]
%      order  Spiral ordering:                                  (1)
%             1=subsequent, 2=golden anngle, 3=random golden angle
%             if 3->order(2)=type (input to golden_angle.m)     ([... 1])
%             if dim==4->order(3)=ordering in z-PE
%               order(3)=4 -> center out z phase encoding
%do_rot_file  Write single static spiral into gradient waveform file and 
%             rotate via vap_phiXX.fdl and vap_thetaXX.fdl      (true)
%        dda  Dummy scans to discard                            ([0 1])
%  num_echos  Number of echos                                   (1)
%  rew_scale  Scale (down) rewinder to reduce PNS               (1)
%             in a range [0, 1]
% grad_coil   Scanner  Gradient coil                            ([])
%                            chronaxie rheobase alpha  gmax   smax
%                            [s]       [T/s]    [m]    [mT/m] [T/m/s]
%   MR750    XRM             334d-6    23.4     0.333  50     200
%   MR750w   XRMW            360d-6    20.0     0.324  33     120
%   HDx      TRM WHOLE       370d-6    23.7     0.344  23     77
%   HDx      TRM ZOOM        354d-6    29.1     0.309  40     150
%   Rio      HRMW            642.4d-6  17.9     0.310  70     150
%   7T       HRMB            359d-6    26.5     0.370  100    200
%
%          k  k-Space trajectory (normalised to +-0.5)
%        dcf  Density compensation function (calculated with voronoi_area)
%          t  Time                                    [s]
%       grad  Gradient trajectory                     [T/m/s]
%        ind  Index for acquired data
%        out  Output structure of wrt_wavs
%      fname  Filename; if nargout==1 -> return fname
%
% 10/2021 Rolf Schulte
% 12/2021 Carolin Pirkl
% 07/2023 Carolin Pirkl
% 03/2024 Kylie Yeung (center out slice ordering)


if (nargin<1), help(mfilename); return; end

%% input + default parameters
fufa_gmax = 0.99;
fufa_smax = 0.99;
scale = 1;
%dda nk_off = 2;                  % skip first acq points (artefact)

gdt = 4e-6;                  % gradient update time [s]
if ~exist('mtx','var'),      mtx = []; end     % matrix size
if isempty(mtx),             mtx = 32; end
if ~exist('arms','var'),     arms = []; end    %spiral arms
if isempty(arms),            arms = 1; end
if length(arms) == 1,        arms = [arms, arms]; end
if ~exist('kdt','var'),      kdt = []; end     % sampling dwell time [us]
if isempty(kdt),             kdt = 16; end
if ~exist('fname','var'),    fname = []; end   % for waveform file
if ~exist('gmax','var'),     gmax = []; end
if isempty(gmax),            gmax = 33; end    % max grad amplitude [mT/m]
if ~exist('smax','var'),     smax = []; end
if isempty(smax),            smax = 120; end   % max slew rate [T/m/s]
if ~exist('nucleus','var'),  nucleus = ''; end
if isempty(nucleus),         nucleus = '1H'; end
if ~exist('dim','var'),      dim = []; end
if isempty(dim),             dim = 3; end
if ~exist('typ','var'),      typ = []; end
if isempty(typ),             typ = 0; end
if ~exist('sampfac','var'),     sampfac = []; end
if isempty(sampfac),            sampfac = 1; end
if length(sampfac)<2,           sampfac = [sampfac 1]; end
if length(sampfac)<3,           sampfac = [sampfac 0]; end
if length(sampfac)<4,           sampfac = [sampfac 1]; end
if ~exist('rewinder','var'), rewinder = []; end
if isempty(rewinder),        rewinder = 0; end
if length(rewinder)==1,      rewinder(1,2) = 1; end
rewinder = rewinder(:).';
if ~exist('order','var'),    order = []; end
if isempty(order),           order = 1; end
if length(order)<2,          order = [order 1]; end
if length(order)<3,          order = [order 1]; end
if ~exist('do_rot_file','var'), do_rot_file = []; end
if isempty(do_rot_file),     do_rot_file = true; end
%if dim==2,                   do_rot_file = false; end
if ~exist('dda','var'),      dda = []; end
if isempty(dda),             dda = 0; end
if length(dda)<2,            dda = [dda 1]; end
switch dim
    case {2},     save_reduced = false;
    case {4},     save_reduced = do_rot_file;  % save only single view
    otherwise,  error('dim(=%g)~=2 or 4',dim);
end
if ~exist('redundant', 'var'),  redundant = []; end
if isempty(redundant),      redundant = 0; end
if (redundant && typ==0),   redundant = 0;  fprintf('Setting redundant=0\n'); end

if ~exist('man_corr_dcf', 'var'),  man_corr_dcf = []; end
if isempty(man_corr_dcf),      man_corr_dcf = 0; end
if ~exist('num_echoes','var'),    num_echoes = []; end
if isempty(num_echoes),           num_echoes = 1; end
if ~exist('rew_scale','var'),    rew_scale = []; end
if isempty(rew_scale),           rew_scale = [1,1]; end
if length(rew_scale) == 1,       rew_scale = [rew_scale, rew_scale]; end
if ~exist('grad_coil','var'),    grad_coil = []; end

%% checks
if gmax<1, warning('gmax(=%g)<1',gmax); input(''); end
if fov<1,  warning('fov(=%g)<1',fov);   input(''); end
if length(mtx)<1,      error('length(mtx)<1'); end
if gmax<1,             warning('gmax(=%g)<1',gmax); input(''); end
if abs(gdt-4d-6)>1d-10, warning('gdt(=%g[us])~=4[us]',gdt*1d6); end
if kdt<1.9,             warning('kdt(=%g)<1.9[us]',kdt); end
if abs(round(kdt/2)*2-kdt)>1d-10
    warning('kdt(=%g) must be multiple of 2[us]',kdt); 
end
max_acq_pts = 16382;                   % maximum number of scans


%% convert to standard SI units
fov = fov*1d-3;                        % [mm] -> [m]
kdt = kdt*1d-6;                        % [us] -> [s]
bw = 1/kdt;                            % full (+-) acq bandwidth [Hz]
if bw<5d3, warning('bw(=%g)<5000[Hz]',bw); end
gmax = gmax*1d-3;                      % [mT/m] -> [T/m]

res = fov./mtx;                        % resolution[m]
gamma = abs(gyrogamma(nucleus))/2/pi;  % gyromagnetic ratio [Hz/T]
fprintf('gamma =%g\n',gamma);
fprintf('FOV (fov) = %g [m]\n',fov); 
fprintf('res (res) = %g [mm]\n',res*1d3);

rgamma = abs(gyrogamma('1h')/gyrogamma(nucleus));
fov_mns = fov./rgamma;
res_mns = res./rgamma;
fprintf('gamma ratio=%g\n',rgamma);
fprintf('scaled FOV (fov_mns) = %g [m]\n',fov_mns); 
fprintf('scaled res (res_mns) = %g [mm]\n',res_mns*1d3);
gmax_nyquist = 2*pi/(gyrogamma('1h')*kdt*fov_mns(1));
fprintf('gmax_nyquist = %g [mT/m]\n',gmax_nyquist*1d3);
if (gmax_nyquist<gmax)
    fprintf('Attention: approaching sampling BW limited regime\n');
    pause(1);
    if (kdt>gdt)
        fprintf('!!! Undersampling will occur: reduce kdt !!!\n'); 
        input('press key to continue');
    end
end

kmax = 1/(2*res_mns(1));
nk_off = 1;                            % skip first acq points (artefact)
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
ng_off = nk_off*dki;                   % #grad pts offset for start of traj
if abs(ceil(ng_off)-ng_off)>1d-10
    warning('ng_off(=%g) not even; rounding',ng_off);
end
ng_off = round(ng_off);
fprintf('gdt=%g[us]; kdt=%g[us]; ',gdt*1d6,kdt*1d6);
fprintf('nk_off=%g; ng_off=%g\n',nk_off,ng_off);


%% design single spiral
% spiral
figure(); clf;
[~,g1,~,t1] = vds(fufa_smax*smax*1e2,fufa_gmax*gmax*1e2,gdt,...
    arms(1)/sampfac(2),fov_mns(1)*1e2*[1,sampfac(3)],1/(2*res_mns(1)*1e2));
gspir = 1e-2*g1;                      % convert to SI unit [T/m]

if typ==1 %spiral in out
    fprintf('-----> spiral in out');
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
% calculate k-space for single spiral
kspir = gdt*cumsum(gspir);   % [s*T/m]
if typ==1
    kspir = kspir(1,(ng_pre+1):end); 
end


%rewindero
smax_rew = smax*rew_scale(1);
switch rewinder(1)
    case 0         % fast grad ramp down
        nrew = ceil(abs(gspir(1,end))/smax_rew/gdt/fufa_gmax)+1;
        grew = linspace(real(gspir(1,end)),0,nrew) + ...
            1i*linspace(imag(gspir(1,end)),0,nrew);
        balstr = '';
    case 1         % return to centre of k-space
        [grewx,ngx] = gradient_lobe(-real(kspir(1,end)),real(gspir(1,end)),...
            fufa_gmax*gmax/sqrt(2),fufa_smax*smax_rew/sqrt(2),gdt);
        [grewy,ngy] = gradient_lobe(-imag(kspir(1,end)),imag(gspir(1,end)),...
            fufa_gmax*gmax/sqrt(2),fufa_smax*smax_rew/sqrt(2),gdt);
        
        ng = max([ngx,ngy])+1;
        ng_rew = max([ngx,ngy]);
        grew = [grewx zeros(1,ng-ngx)] + 1i*[grewy zeros(1,ng-ngy)];
        balstr = '_blncd';
    case 2         % add crusher
        ga_rew = rewinder(2)*abs(kspir(1,end))*gspir(1,end)/abs(gspir(1,end));
        [grewx,ngx] = gradient_lobe(real(ga_rew),real(gspir(1,end)),...
            abs(real(gspir(1,end))),fufa_smax*smax_rew/sqrt(2),gdt);
        [grewy,ngy] = gradient_lobe(imag(ga_rew),imag(gspir(1,end)),...
            abs(imag(gspir(1,end))),fufa_smax*smax_rew/sqrt(2),gdt);
        
        ng = max([ngx,ngy])+1;
        ng_rew = max([ngx,ngy]);
        grew = [grewx zeros(1,ng-ngx)] + 1i*[grewy zeros(1,ng-ngy)];
        balstr = '_crush';
    case 3         % return partially to centre of k-space
        ga_rew = -rewinder(2)*kspir(1,end);
        [grewx,ngx] = gradient_lobe(real(ga_rew),real(gspir(1,end)),...
            fufa_gmax*gmax/sqrt(2),fufa_smax*smax_rew/sqrt(2),gdt);
        [grewy,ngy] = gradient_lobe(imag(ga_rew),imag(gspir(1,end)),...
            fufa_gmax*gmax/sqrt(2),fufa_smax*smax_rew/sqrt(2),gdt);
        
        ng = max([ngx,ngy])+1;
        ng_rew = max([ngx,ngy]);
        grew = [grewx zeros(1,ng-ngx)] + 1i*[grewy zeros(1,ng-ngy)];
        balstr = sprintf('_ret%g',rewinder(2));
    otherwise
        error('rewinder(1)(=%g)~=0,1,2 or 3',rewinder(1));
end


% construct full (not yet rotated) gradient waveform
gg = [zeros(1,ng_off),gspir,grew].';
if isodd(size(gg,1))                 % adding 0 if odd
    gg = [gg ; 0];
end
ng = size(gg,1);
nk = ceil(length(kspir)/dki);

%% construct 3D spiral through rotation matrices
switch dim
    % 3D rotation in fidall:
    % rotation of theta around y then of phi around z
    % -> put spirals into xz-plane -> swap phi & theta
   
    case {2,4} %2D & 2D+PEz -> phi
        % #spiral rotations (2D)
        if length(arms)>1
            Ns = arms(2); %*arms(2);                     
        else
            Ns = ceil(pi*kmax*fov(1)*sampfac(4))*2;     
        end
        
        
        fprintf('#spiral rotations in 2D = %g\n', Ns);
        %fprintf('#full spiral rotations = %g\n', arms(2));
        
        switch order(1)
            case 1 %consecutive
                phi = (0:Ns-1)/Ns*2*pi;
                phi = mod(phi,2*pi);
                phi = phi(:).';
                
            case {2,3} %golden angle + random
                ga = golden_angle(order(2));
                phi = (0:Ns-1)*ga;
                phi = mod(phi,2*pi);
                
            	if order(1)==3 %random
                    ii = randperm(Ns);
                    phi = phi(1,ii);
                end
        end
        
%         remain = mod(arms(2), arms(1));
%         repeats = fix(arms(2)/ arms(1));
%         tmp = zeros(1, arms(2));
%         tmp(1:repeats*arms(1))=repmat(phi,1,repeats);
%         tmp(repeats*arms(1)+1:end) = phi(1:remain);
%         phi = tmp;
end
    
if redundant
    fprintf('Designing redunant spiral shape\n');
    Ns = 2*Ns;
    phi = [phi, phi+pi];
end


%% add z-phase encoding (stack of spirals)
if dim==4
    if length(fov)<2, error('length(fov)(=%d)<2',length(fov)); end
    if length(mtx)<2, error('length(mtx)(=%d)<2',length(mtx)); end
    if res(2)<res(1)
        warning('res(2)(=%g)<res(1)(=%g)',res(2),res(1));
    end
    
    if mod(mtx(2),2) == 0
        agz = ndgrid((-mtx(2)/2:mtx(2)/2-1)/mtx(2)*2);
    else
        agz = -ndgrid((-(mtx(2)-1)/2:(mtx(2)-1)/2)/(mtx(2)-1)*2); %% added for symmetric sampling in z (KY 25/03/2024)
    end

    if order(3)==1      % subsequent
        agz = repelem(agz.',1,length(phi));
    elseif order(3)==4 %% center out slices sampling (KY 14/03/2024)
        if mod(length(agz),2)==1
            tmp=[0 agz(ceil(end/2)+1:end)';flip(agz(1:ceil(end/2)))']; % stick them together with 0 to match dim
            agz=tmp(:); %vectorize
            agz=repelem(agz(2:end).',1,length(phi)); % remove leading zero so right number of slices
        else
            tmp=[agz(ceil(end/2)+1:end)';flip(agz(1:ceil(end/2)))']; % stick them together
            agz=repelem(tmp(:).',1,length(phi)); %vectorize
        end
    else
        tmp2 = repelem(agz.',length(phi),1);
        for l1=1:length(phi)
            switch order(3)
                case 2  % random
                    ii = randperm(mtx(2));
                case 3  % progressing
                    ii = mod((0:mtx(2)-1)+l1-1,mtx(2))+1;
                otherwise, error('order(3)(=%g) unknown',order(3));
            end
            tmp2(l1,:) = tmp2(l1,ii);
        end
        agz = tmp2(:).';
    end
    ga_pha = mtx(2)/fov(2)/gamma/2;
    sder = res(1)/res(2)/sqrt(2);
    [G_pha,ng_pha] = gradient_lobe(ga_pha,0,gmax,sder*smax*rew_scale(2),gdt,true);
    %if ng_pha>ng_pre, warning('ng_pha(=%g)>ng_pre(=%g)',ng_pha,ng_pre); end
end



%% ordering + rotations
ng1 = length(gg);
nintlv = size(phi,2);                  % #interleaves
ng2 = size(phi,2);                     % #interleaves stored in grad

if ~exist('ng_rew', 'var')
    ng_rew = 0;
end
if ~exist('ng_pha', 'var')
    ng_pha = 0;
end

if dim==4, ng2 = ng2*mtx(2); end

if ng2>max_acq_pts*16
    warning('#interleaves(=%g)>%g*16 (fidall limit)',...
        ng2,max_acq_pts);
    input('press key to continue (might also crash matlab)');
end


%% assemble 2D/4D k-space
% k-space trajectory (complex -> 2D)
if dki>1
    kc = kspir(1:dki:length(kspir))/2/pi*res_mns(1)*gyrogamma('1H'); 
else
    kc = kspir/2/pi*res_mns(1)*gyrogamma('1H');
end
if abs(dki-0.5)<1d-10
    kc = interp1((1:length(kc)),kc,(1:0.5:length(kc)+0.5),'spline');
end
kc = kc.';
nks = size(kc,1); 
k = zeros(length(kc), ng2, dim-(dim==4));
switch dim
    case 2 
        k(:,:,1) = real(kc)*cos(phi) + imag(kc)*sin(phi);
        k(:,:,2) = - real(kc)*sin(phi) + imag(kc)*cos(phi);
 
    case 4
        k_tmp = zeros(nks, nintlv, 2);
        k_tmp(:,:,1) = real(kc)*cos(phi) + imag(kc)*sin(phi);
        k_tmp(:,:,2) = - real(kc)*sin(phi) + imag(kc)*cos(phi);
        k(:,:,1) = repmat(k_tmp(:,:,1),[1 mtx(2)]);
        k(:,:,2) = repmat(k_tmp(:,:,2),[1 mtx(2)]);
        kz = agz/2;
        k(:,:,3) = repmat(kz,[nks,1]);   
        theta = [];
end
 
k = reshape(k,[1 nks*ng2 size(k,3)]);
k_pts = size(k,2);


fprintf('Calculating dcf analytically\n');
switch dim
    case 2
        tmp = pi/arms(1)*(k(1,2:end).*conj(k(1,2:end))-k(1,1:end-1).*conj(k(1,1:end-1)))*mtx(1)^2;
        dcf_thresh = max(tmp);
        if redundant
            dcf = voronoi_area(k(:,1:end/2,:)*mtx(1),'','',dcf_thresh);
            dcf = [dcf, dcf];
        else
            dcf = voronoi_area(k*mtx(1),'','',dcf_thresh);
        end
        
        figure();
        plot(dcf(1,1:size(dcf,2)/Ns));
    case 4 
        remain = mod(arms(2), arms(1));
        repeats = fix(arms(2)/ arms(1));
        if remain>0
            repeats = repeats+1;
        end
        %only consider one fully sampled spiral (=arms(1)) after another
        for rep = 1:repeats
            if rep==repeats && remain>0
                arms_tmp = remain;
            else
                arms_tmp = arms(1);
            end

            if redundant
                k_tmp2 = k_tmp(:,(rep-1)*arms_tmp*2+1:rep*arms_tmp*2,:);
                k_tmp2 = reshape(k_tmp2, [1, nks*arms_tmp*2 size(k_tmp2,3)]);
            else
                k_tmp2 = k_tmp(:,(rep-1)*arms_tmp+1:rep*arms_tmp,:);
                k_tmp2 = reshape(k_tmp2, [1, nks*arms_tmp size(k_tmp2,3)]);
            end
            tmp = pi/arms_tmp*(k_tmp2(1,2:end).*conj(k_tmp2(1,2:end))-k_tmp2(1,1:end-1).*conj(k_tmp2(1,1:end-1)))*mtx(1)^2;
            dcf_thresh = max(tmp);
    
            if redundant
                dcf_tmp = voronoi_area(k_tmp2(:,1:end/2,:)*mtx(1),'','',dcf_thresh);
                dcf_tmp = [dcf_tmp, dcf_tmp];
            else
                dcf_tmp = voronoi_area(k_tmp2*mtx(1),'','',dcf_thresh);
            end
            if rep==1
                dcf = dcf_tmp;
            else
                dcf = [dcf, dcf_tmp];
            end
            clear dcf_tmp k_tmp2
        end

            
        %scale according to arms(2) <- voronoi_area calc breaks down for
        %dense spirals (e.g. with arms(2))
        %scale = arms(2)/arms(1);

        if man_corr_dcf
            fprintf('Manually cutting dcf!!!\n');
            dcf_tmp = reshape(dcf, size(k_tmp,1), size(dcf,2)/size(k_tmp,1));
            cut_idx = round(size(dcf_tmp,1)*0.85);
            dcf_tmp(cut_idx:end, :) = repmat(dcf_tmp(cut_idx,:), [size(dcf_tmp(cut_idx:end, :),1),1]);
            dcf_tmp = reshape(dcf_tmp, size(dcf));
            dcf  = dcf_tmp; 
            clear dcf_tmp
        end
        scale = arms(2)/arms(1);
        dcf = dcf/scale;
        
        dcf = repmat(dcf,[1 mtx(2)]);
        figure();
        plot(dcf(1,1:size(dcf,2)/size(agz,2)));
end


%% calculate time + index
if dim==4
    acq_pts = (length(gspir)+ng_off+ng_pha+ng_rew)/dki;
else
    acq_pts = (length(gspir)+ng_off)/dki;
end


inds = false(1,round(length(gg)/dki));                   % index for acq data (single arm)
inds(1,(nk_off+(1:nk))) = true;

%calculate echo time TE
te = nk_off*kdt*1d3; %[ms]
if typ>0
    te = te + t_arm*1d3;
end
if num_echoes>1
    del_te = kdt*size(gg,1)*1d3;
    te = (0:del_te:(num_echoes-1)*del_te) + te;
end
if dim==4
    te = te + kdt*ng_pha*1d3;
end

%repeat waveforms for multiecho acq
if num_echoes >1
    fprintf('Multiecho calculation\n');
  
    gtmp = gg;
    %itmp = inds;
        for i = 1:num_echoes-1
            gtmp = cat(1,gtmp, gg);
            %itmp = cat(2,itmp, inds);
        end
    gg = gtmp;
    %inds = itmp;  
    ng1 = size(gg,1);
end


if do_rot_file
    switch dim
        case 2
            grad = zeros(ng1, 1, 2);
            grad(:,:,1) = real(gg);
            grad(:,:,2) = imag(gg);
        case 4
            if typ>0
                grad = zeros(ng1+2*ng_pha,1,3);
                grad(ng_pha+1:ng_pha+ng1,:,1) = real(gg);          % waveform along x and y
                grad(ng_pha+1:ng_pha+ng1,:,2) = imag(gg);  
                grad(1:ng_pha,1,3) = G_pha; 
                if rewinder(1)==1
                    grad(ng_pha+ng1+1:end,:,3) = -G_pha;
                end
            else
                grad = zeros(ng1+ng_pha+max(0,ng_pha-ng_rew),1,3);
                grad(ng_pha+1:ng1+ng_pha,:,1) = real(gg);          % waveform along x and y
                grad(ng_pha+1:ng1+ng_pha,:,2) = imag(gg);                
                grad(1:ng_pha,1,3) = G_pha;      % phase encoding along z
                if rewinder(1)==1
                    grad(end-ng_pha+1:end,:,3) = -G_pha;
                end
            end

    end
else
   
    switch dim
        case 2
            grad = zeros(ng1,ng2,dim-(dim==4));
            grad(:,:,1) = real(gg)*cos(phi) + imag(gg)*sin(phi);
            grad(:,:,2) = - real(gg)*sin(phi) + imag(gg)*cos(phi);
        case 4
            grad = zeros(ng1+ng_pha-max(0,ng_pha-ng_rew),ng2,dim-(dim==4));
            grad(ng_pha+1:end-max(0,ng_pha-ng_rew),:,1) = repmat(real(gg)*cos(phi) + imag(gg)*sin(phi), [1 mtx(2)]);
            grad(ng_pha+1:end-max(0,ng_pha-ng_rew),:,2) = repmat(- real(gg)*sin(phi) + imag(gg)*cos(phi), [1 mtx(2)]);
            grad(1:ng_pha,:,3) = bsxfun(@times,G_pha.',agz);
            if rewinder(1)==1
                ii = (1:ng_pha)+ng1-abs(ng_pha-ng_rew);
                grad(ii,:,3) = -bsxfun(@times,G_pha.',agz);
            end
    end
    
end

nk_pha = round(ng_pha/dki);
inds = padarray(inds, [0,nk_pha], false,'pre');

if ~isempty(grad_coil)
    figure(); 
    pns(grad, grad_coil);
end


ts = (0:nk-1)*kdt;                         % time of single trajectory [s]


ind = repmat(inds,[nintlv 1]);             % index acq all interleaves
t = repmat(ts,[nintlv 1]);                 % time of all trajectories [s]

acq_pts = size(grad,1)/dki;
acq_pts = ceil(acq_pts/32)*32;             % round to mulmulttiple of 32

if acq_pts>max_acq_pts
    warning('#sampling pts/exc (=%g)>%g: exceeds fidall limit',acq_pts, max_acq_pts);
end

%% plot
switch dim
    case 2
        figure(); 
        plot(grad(:,1,1));
        hold on; 
        plot(grad(:,1,2));
        hold on; 
        plot(inds);
    otherwise
        figure(); 
        plot(grad(:,1));
        hold on; 
        plot(grad(:,2));
        hold on; 
        plot(grad(:,3));
        plot(inds);
end

%% print waveform info
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

if typ>0, fprintf('TE += %g [ms]\n',te); end


%% checks
if dim==2
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
if max(dcf)>1
    warning('max(dcf)(=%g)>1',max(dcf));
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
        switch dim
            case 4
                fov_str = sprintf('%g_%g',round(fov(1)*1d3),round(fov(2)*1d3));
                mtx_str = sprintf('%d_%d',mtx(1),mtx(2));
            otherwise
                fov_str = num2str(round(fov(1)*1d3));
                mtx_str = num2str(mtx(1));
        end
        
        if length(arms)>1
            arm_str = sprintf('%d_%d',arms(1),arms(2));
        else
            arm_str = num2str(arm(1));
        end
    
        str_typ = '';
        if typ >0
            str_typ = [str_typ 'inout'];
        else
            str_typ = [str_typ 'out'];
        end
        if redundant
            fname = sprintf('redundant_spiral%gD_%s_%s_fov%s_mtx%s_arms%s_kdt%g_gmax%g_smax%g_dur%g%s_numecho%s', ...
            dim,str_typ,nuc,fov_str,mtx_str,arm_str,kdt*1d6,...
            round(gmax_act*1d3),round(smax_act),...
            round(ng*gdt*1d4)/10,balstr, num2str(num_echoes));
            if ~(scale==1)
                fname = sprintf([fname '_scale']);
            end
        else
            fname = sprintf('spiral%gD_%s_%s_fov%s_mtx%s_arms%s_kdt%g_gmax%g_smax%g_dur%g%s_numecho%s', ...
            dim,str_typ,nuc,fov_str,mtx_str,arm_str,kdt*1d6,...
            round(gmax_act*1d3),round(smax_act),...
            round(ng*gdt*1d4)/10,balstr, num2str(num_echoes));
            if ~(scale==1)
                fname = sprintf([fname '_scale']);
            end
        end
        fname = strrep(fname,'.','p');
        if abs(gdt*1d6-4)>1d-10, fname = [fname '_gdt' num2str(gdt*1d6)]; end
    
        switch order(1)
            case 2, fname = [fname '_GA' num2str(order(2))];
            case 3, fname = [fname '_GA' num2str(order(2)) '_RND'];
        end
        if ((dim==4) && (order(3)>1))
            fname = [fname '_zPE' num2str(order(3))];
        end
        
    else
        fname = [];
    end
end

%% export waveform
if ~isempty(fname)   
    switch dim
        case 4
            fov_str = sprintf('%g_%g',round(fov(1)*1d3),round(fov(2)*1d3));
            mtx_str = sprintf('%d_%d',mtx(1),mtx(2));
        otherwise
            fov_str = num2str(round(fov(1)*1d3));
            mtx_str = num2str(mtx(1));
    end
    
    if length(arms)>1
        arm_str = sprintf('%d_%d',arms(1),arms(2));
    else
        arm_str = num2str(arm(1));
    end
    
    sampfac_str = [num2str(sampfac(1))];
    for i=2:length(sampfac)
        sampfac_str = [sampfac_str '_' num2str(sampfac(i))];
    end
    
    order_str = [num2str(order(1))];
    for i=2:length(order)
        order_str = [order_str '_' num2str(order(i))];
    end
    
    desc = sprintf('design_spiral_mrcp(%s,%s,%s,%g,%g,%g,%g,%s,%g,%g,%s,%g,%s,%g,%g)',...
        fov_str, mtx_str, arm_str, kdt*1d6, ~isempty(fname),gmax*1d3,...
        smax, nucleus, dim, typ, sampfac_str, rewinder(1), order_str, ...
        do_rot_file, dda(1));
        
    phi = phi*180/pi;                        % convert [rad] -> [deg]
    npix = mtx(1);
    if dim==2
        if length(mtx)==1, mtx = [mtx mtx]; end
    else
        switch length(mtx)
            case 1, mtx = repmat(mtx,[1 3]);
            case 2
                mtx = [mtx(1) mtx(1) mtx(2)];
                if length(fov)==2
                    fov = [fov(1) fov(1) fov(2)];
                else
                    warning('length(mtx)==2 && length(fov)~=2');
                end
            otherwise, warning('length(mtx)(=%d)~=1 or 2',length(mtx));
        end
    end
    
    fprintf('Saving gradient + mat file to\n\t%s\n',fname);
    fprintf('%s\n',desc);
   
    
    if exist([fname '.wav'],'file')
        warning('file (=''%s'') exists; overwriting',[fname '.wav']);
    end
    
    out = write_ak_wav([fname '.wav'],grad,1/kdt,fov_mns(1),desc,acq_pts,...
        [],[],gdt*1d6);
    if isfield(out,'wave'), out = rmfield(out,'wave'); end
    if isfield(out,'grad'), out = rmfield(out,'grad'); end
    
    if (save_reduced && dim==4)
        nviews = ng2;
        ks = zeros(1,nk,3);
        ks(1,:,1) = real(kc(1:nk,1)).';
        ks(1,:,2) = imag(kc(1:nk,1)).';
        
        
        dcf_range = max(dcf);
        dcf = single(dcf);
        
        save([fname '.mat'],'out','ks','ts','inds','nviews','fov','mtx',...
            'gdt','kdt','gmax','smax','nucleus','t_rew','t_arm',...
            'ng_off','gmax_nyquist','rgamma','fov_mns','phi',...
            'dcf','dcf_range', 'theta', 'kz', 'te', 'nk_pha');
    else
        dcf = single(dcf); k = single(k); t = single(t);
        switch dim
            case 2
                save([fname '.mat'],'out','k','dcf','t','ind','fov','mtx',...
                    'gdt','kdt','gmax','smax','nucleus','t_rew','t_arm',...
                    'ng_off','gmax_nyquist','rgamma','fov_mns','npix', 'phi', 'te');
            case 4
                save([fname '.mat'],'out','k','dcf','t','ind','fov','mtx',...
                    'gdt','kdt','gmax','smax','nucleus','t_rew','t_arm',...
                    'ng_off','gmax_nyquist','rgamma','fov_mns','npix', 'phi', ...
                    'theta', 'kz', 'te', 'nk_pha');
        end
    end
    

    if do_rot_file
        write_fdl([fname '_phi.fdl'],mod(phi,360),'phi');
        if dim==4
            write_fdl([fname '_agz.fdl'],agz,'agrad');
        end
    end

end    

%% output arguments
if nargout==1
    k = fname;
end


