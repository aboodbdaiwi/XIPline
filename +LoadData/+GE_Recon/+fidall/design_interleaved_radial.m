function [k,dcf,t,grad,ind,fname] = design_interleaved_radial(fov,mtx,kdt,fname,...
    gmax,smax,nucleus,dim,typ,sampfac,rewinder,order,do_rot_file,dda,flips,freqs,bonus_fids)
%DESIGN_INTERLEAVED_RADIAL  Design isotropic 2D or 3D radial sampling trajectory
% [k,dcf,t,grad,ind,fname] = design_interleaved_radial(fov,mtx,kdt,fname,gmax,smax,...
%            nucleus,dim,typ,sampfac,rewinder,order,do_rot_file,dda)
%                                                     [unit]    (default)
%        fov  Field-of-view                           [mm]
%        mtx  Matrix size                             [1]
%        kdt  Sampling dwell time (1/bw)              [us]      (16)
%      fname  Filename (string, or logical)                     (false)
%       gmax  Maximum gradient strength               [mT/m]    (33)
%       smax  Maximum slewrate                        [T/m/s]   (120)
%    nucleus  Nucleus                                           ('1H')
%        dim  Spatial dimensions (2=2D,3=3D,4=2D+z-PE)          (3)
%        typ  Radial variant:                                   (0)
%             0=centre-out constant; 1=centre-out density-adapted
%             2=full-spoke constant
%    sampfac  Sampling factor (>1 over-, <1 under-sampled)      (1.5 1)
%             1=frequency encoding direction (speed in k-space)
%             2=spokes spacing@kmax=sqrt(opnex/pi) in 3dradial
%               if sampfac(2)>10 -> #spokes=sampfac(2)
%   rewinder  Gradient ramp down method:                        (0)
%             0=fast, 1=return to centre of k-space, 2=crushing
%             2->[2 area-factor(default=1)]
%      order  Spokes ordering:                                  (1)
%             equidistant:1=subsequent,2=random; 3=golden angle; 
%             4=halton randomized spiral
%             if 3->order(2)=type (input to golden_angle.m)     ([... 1])
%                   order(4)=0->tangential correct; 1->equidistant dcf
%             if dim==4->order(3)=ordering in z-PE
%do_rot_file  Write single static spoke into gradient waveform file and 
%             rotate via vap_phiXX.fdl and vap_thetaXX.fdl      (true)
%        dda  Dummy scans to discard                            ([0 1])
%      flips  flip angles for interleaves                       ([])
%      freqs  frequencies for interleaves                       ([])
%
%          k  k-Space trajectory (normalised to +-0.5)
%             #acq-pts #interleaves #dim
%        dcf  Density compensation function (analytical)
%             assumes equal spokes distribution; some error for 3D mtx<32
%          t  Time                                    [s]
%       grad  Gradient trajectory                     [T/m/s]
%        ind  Index for acquired data
%
% Literature: Handbook of MRI Pulse Sequence; Bernstein, King, Zhou
%             mrm65_1091_2011_konstandin.pdf,mrm62_1565_2009_nagel.pdf
% 11/2022 Rolf Schulte
if (nargin<1), help(mfilename); return; end


%% default parameters
gdt = 4e-6;                  % gradient update time [s]
if ~exist('mtx','var'),      mtx = []; end     % matrix size
if isempty(mtx),             mtx = 32; end
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
if ~exist('sampfac','var'),  sampfac = []; end
if isempty(sampfac),         sampfac = [1.5 1]; end
if ~exist('rewinder','var'), rewinder = []; end
if isempty(rewinder),        rewinder = 0; end
if length(rewinder)==1,      rewinder(1,2) = 1; end
rewinder = rewinder(:).';
if ~exist('order','var'),    order = []; end
if isempty(order),           order = 1; end
if length(order)<2,          order = [order 1]; end
if length(order)<3,          order = [order 1]; end
if length(order)<4,          order = [order 1]; end
if ~exist('do_rot_file','var'), do_rot_file = []; end
if isempty(do_rot_file),     do_rot_file = true; end
if ~exist('dda','var'),      dda = []; end
if isempty(dda),             dda = 0; end
if length(dda)<2,            dda = [dda 1]; end
switch dim
    case 2,     save_reduced = false;
    case {3,4}, save_reduced = do_rot_file;  % save only single view
    otherwise,  error('dim(=%g)~=2,3 or 4',dim);
end

if ~exist('flips','var'),      flips = []; end   
if isempty(flips),             var_flips = false; else, var_flips = true; end
if ~exist('freqs','var'),      freqs = []; end  
if isempty(freqs),             var_freqs = false;  else, var_freqs = true; end
if ~exist('bonus_fids','var'),      bonus_fids = []; end  
if isempty(bonus_fids),             bonus_fids = 0; end


%% checks
if length(fov)<1,      error('length(fov)<1'); end
if any(fov<1),         warning('fov(=%g)<1',fov(1)); end
if length(mtx)<1,      error('length(mtx)<1'); end
if gmax<1,             warning('gmax(=%g)<1',gmax); input(''); end
if length(sampfac)~=2, error('length(sampfac)~=2'); end
if abs(gdt-4d-6)>1d-10, warning('gdt(=%g[us])~=4[us]',gdt*1d6); end
if kdt<1.9,             warning('kdt(=%g)<1.9[us]',kdt); end
if abs(round(kdt/2)*2-kdt)>1d-10
    warning('kdt(=%g) must be multiple of 2[us]',kdt); 
end

if var_flips && var_freqs
    if length(flips)~=length(freqs)
        error('flip and freq interleaves do not match');
    end
    intlvs = length(flips);
    if length(bonus_fids)==1 && length(bonus_fids)~=intlvs
        bonus_fids = repmat(bonus_fids,1,intlvs);
    end
end
if sum(bonus_fids(:))>0 && do_rot_file
    error('bonus fids not implemented with rot files');
end
if sum(bonus_fids(:))>0 && dim == 4
    error('bonus fids not implemented for z-PE');
end


%% misc parameters
RFS_MAX_NSCANS = 16382;                % maximum number of scans
fov = fov*1d-3;                        % [mm] -> [m]
kdt = kdt*1d-6;                        % [us] -> [s]
bw = 1/kdt;                            % full (+-) acq bandwidth [Hz]
if bw<5d3, warning('bw(=%g)<5000[Hz]',bw); end
gmax = gmax*1d-3;                      % [mT/m] -> [T/m]
res = fov./mtx;                        % resolution         [m]
kmax = 1/(2*res(1));
GAM = abs(gyrogamma(nucleus))/2/pi;    % gyromagnetic ratio [Hz/T]
nk_off = 2;                            % skip first acq points (artefact)
ng_off = nk_off*kdt/gdt;               % #grad pts offset for start of traj
if abs(ceil(ng_off)-ng_off)>1d-10
    warning('ng_off(=%g) not even; rounding',ng_off);
end
ng_off = round(ng_off);
fprintf('gdt=%g[us]; kdt=%g[us]; ',gdt*1d6,kdt*1d6);
fprintf('nk_off=%g; ng_off=%g\n',nk_off,ng_off);


%% readout gradient strength of single spoke
g_read = bw/(fov(1)*GAM)/sampfac(1);   % [T/m]
if g_read>gmax
    warning('g_read(=%g[mT/m])>gmax(=%g[mT/m]; setting g_read=gmax',...
        g_read*1d3,gmax*1d3);
    fprintf('Better to increase kdt (less oversampling)\n');
    g_read = gmax;
end
fprintf('g_read=%g[mT/m]\n',g_read*1d3);


%% construct single spoke
switch typ
    case 0    % centre-out constant grad strength
        ng_ramp1 = ceil(g_read/smax/gdt)+1;       % #ramp-up points
        ng_plat = ceil(kmax/(GAM*g_read*gdt)-ng_ramp1/2+2*kdt/gdt); 
        % #plateau pts
        if ng_plat<0
            warning('ng_plat(=%g)<0; seeting to 0',ng_plat);
            ng_plat = 0;
        end
        fprintf('ng_plat=%g; ng_ramp1=%g\n',ng_plat,ng_ramp1);
        gs = [zeros(ng_off,1) ; ...
            linspace(0,g_read,ng_ramp1).' ; ...
            g_read*ones(ng_plat,1)];
        ng_pre = 0;
    case 1    % centre-out density-adapted radial
        gs = sub_da_spoke(g_read,smax,kmax,gdt,GAM,dim,ng_off);
        ng_pre = 0;
    case 2    % full spoke constant grad strength, no ramp sampling
        ng_plat = 2*ceil(kmax/(GAM*g_read*gdt)); % #plateau pts
        ga_pre = -g_read*ng_plat*gdt/2;          % area prewinder
        [G_pre,ng_pre] = gradient_lobe(ga_pre,g_read,gmax,smax,gdt,true);
        fprintf('ng_plat=%d; ng_pre=%d\n',ng_plat,ng_pre);
        
        gs = [zeros(ng_off,1) ; ...
            G_pre(1,end:-1:1).' ; ...
            g_read*ones(ng_plat,1)];
        
    otherwise
        error('unknown typ(=%d)',typ);
end
switch rewinder(1,1)
    case 0    % regular fast ramp down
        ng_ramp2 = ceil(gs(end,1)/smax/gdt)+1;   % #ramp-down points
        gs = [gs ; linspace(gs(end,1),0,ng_ramp2).'];
    case 1    % return to centre of k-space and ramp down
        ga_rew = -sum(gs,1)*gdt;                 % area rewinder
        [G_rew,ng_ramp2]  = gradient_lobe(ga_rew,gs(end,1),gmax,smax,gdt,true);
        gs = [gs ; G_rew.'];
    case 2    % crushing
        ga_rew = rewinder(1,2)*sum(gs,1)*gdt;    % area crushing
        [G_rew,ng_ramp2]  = gradient_lobe(ga_rew,gs(end,1),gmax,smax,gdt,true);
        gs = [gs ; G_rew.'];
end
if isodd(size(gs,1))
    gs = [gs ; 0];                               % add 0 for even #grad pts
    ng_ramp2 = ng_ramp2+1;
end


%% construct 3D radial through rotation matrices
switch dim
    case {2,4}
        if sampfac(2)>10
            Ns = round(sampfac(2));
        else
            Ns = ceil(pi*kmax*fov(1)*sampfac(2))*2;
            if typ==2, Ns = round(Ns/2); end
            % full spoke must use 360� for artefact reduction
        end
        if order(1)==3       % golden angle
            fprintf('Number of spokes (fully sampled): Ns=%g\n',Ns);
            l = 1;
            while false       % rounding up to next Fibonacci number
                if Ns<fibonacci(l)
                    Ns = fibonacci(l);
                    break;
                else
                    l = l+1;
                end
            end
            % fprintf('\trounding up to next Fibonacci number: Ns=%g\n',Ns);
            ga = golden_angle(order(2));
            phi = (0:Ns-1)*ga;
            phi = mod(phi,2*pi);
        else
            phi = (0:Ns-1)/Ns*2*pi;
            if (typ==2) && ~isodd(Ns) 
                 phi(1,(Ns/2+1):Ns) = phi(1,(Ns/2+1):Ns) + pi/Ns;
                 % full spoke: use 360� + add half angle step for >180�
            end
        end
    case 3
        if order(1)==3       % golden angle
            % Fyrdahl, MAGMA 2020
            ga1 = sub_fibonacci_series(1,100)/...
                sub_fibonacci_series(order(2),102);
            ga2 = sub_fibonacci_series(1,100)/...
                sub_fibonacci_series(order(2),101);
            fprintf('Golden angle-%d: ga1=%g; ga2=%g\n',order(2),ga1,ga2);

            if sampfac(2)>10
                Ns = round(sampfac(2));
            else
                if typ==2   % full spoke
                    Ns = ceil(pi*(kmax*fov(1)*sampfac(2))^2)*2;
                else        % centre-out spoke
                    Ns = ceil(2*pi*(kmax*fov(1)*sampfac(2))^2)*2;
                end
            end
            if typ==2   % full spoke
                theta = acos(mod(ga1*(0:(Ns-1)),1)); % probably mod(x,2)
            else        % centre-out spoke
                theta = acos(mod(ga1*(0:(Ns-1)),2)-1); 
            end
            phi = 2*pi*mod(ga2*(0:(Ns-1)),1);
        elseif order(1)==4       % randomized halton spiral
            if sampfac(2)>10
                Ns = round(sampfac(2));
            else
                if typ==2   % full spoke
                    Ns = ceil(pi*(kmax*fov(1)*sampfac(2))^2)*2;
                else        % centre-out spoke
                    Ns = ceil(2*pi*(kmax*fov(1)*sampfac(2))^2)*2;
                end
            end
            %calc spiral order
            M_PI = 3.14159265358979323846;
            PrevAngle = 0;
            for proj = 0:Ns-1
                currZ = -1.0 + 2.0 * proj/Ns;
                arch_z(proj+1) = currZ;
                polar_angle(proj+1) = acos(currZ);
                if (proj == 0)
                    arch_azi(proj+1) = 0;
                else
                    arch_azi(proj+1) = mod(PrevAngle + 3.6/(sqrt(Ns*(1-currZ*currZ))), 2.0*M_PI);
                end
                PrevAngle = arch_azi(proj+1);
            end
            
            %calc halton order
            p1 = 2;
            %p2 = 3;
            for proj = 0:Ns-1
                z = haltonnumber(proj+1,p1) * 2 - 1;
                %phi = 2 * M_PI * haltonnumber(proj+1,p2);
                halt_polar(proj+1) = acos(z);
                %halt_azi(proj+1) = phi; %not used here
            end
            
            %sort spiral via halton
            [~,sort_index] = sort(halt_polar);
            phi = arch_azi(:,sort_index);            
            theta = polar_angle(:,sort_index); 
        else
            Ntheta = ceil(pi*kmax*fov(1)*sampfac(2));
            phi = []; theta = [];
            if typ==2   % full spoke
                th_arr = ((0:Ntheta/2)/Ntheta*pi);
                Ns = 2*pi*(kmax*fov(1)*sampfac(2))^2;
            else        % centre-out spoke
                th_arr = ((0:Ntheta)/Ntheta*pi);
                Ns = 4*pi*(kmax*fov(1)*sampfac(2))^2;
            end
            for th = th_arr
                Nphi = ceil(2*pi*kmax*fov(1)*sampfac(2)*sin(th));
                if Nphi<1, Nphi = 1; end
                phi = [phi (0:Nphi-1)/Nphi*2*pi];
                theta = [theta repmat(th,[1 Nphi])];
            end
        end
    otherwise, error('dim(=%g)~=2,3 or 4',dim);
end


%% add z-phase encoding (stack-of-stars)
if dim==4
    if length(fov)<2, error('length(fov)(=%d)<2',length(fov)); end
    if length(mtx)<2, error('length(mtx)(=%d)<2',length(mtx)); end
    if res(2)<res(1)
        warning('res(2)(=%g)<res(1)(=%g)',res(2),res(1));
    end
    agz = ndgrid((-mtx(2)/2:mtx(2)/2-1)/mtx(2)*2);
    if order(3)==1      % subsequent
        agz = repelem(agz.',1,length(phi));
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
    ga_pha = mtx(2)/fov(2)/GAM/2;
    sder = res(1)/res(2)/sqrt(2);
    [G_pha,ng_pha] = gradient_lobe(ga_pha,0,gmax,sder*smax,gdt,true);
    if ng_pha>ng_pre, warning('ng_pha(=%g)>ng_pre(=%g)',ng_pha,ng_pre); end
end


%% checks
if any(isnan(phi(:))), error('phi isnan'); end
if any(isinf(phi(:))), error('phi isinf'); end


%% ordering
ng1 = length(gs);
n2 = length(phi);
if order(1)==2
    ii = randperm(n2);
    phi = phi(ii);
    if dim==3, theta = theta(ii); end
end
if dim==4, n2 = n2*mtx(2); end

%% interleaves
if intlvs>1
    phi = repelem(phi,intlvs); %repeat each angle for each interleave
    if exist('theta','var')
        theta = repelem(theta,intlvs);
    end
end
if sum(bonus_fids(:))>1 % if bonus spec added, want additional (meaningless) phi/theta
    phi_spec = zeros(1,sum(bonus_fids(:)));
    phi = [phi phi_spec];
    if exist('theta','var')
        theta_spec = zeros(1,sum(bonus_fids(:)));
        theta = [theta theta_spec];
    end
end
if dim~=4, n2 = length(phi); end%update with interleaves and bonus spec

%% rotations
if n2>RFS_MAX_NSCANS*16
    warning('#interleaves(=%g)>%g*16 (fidall limit)',...
        n2,RFS_MAX_NSCANS);
    input('press key to continue (might also crash matlab)');
end
if do_rot_file
    switch dim
        case 2, grad = gs;        % first waveform pointing in x
        case 3
            grad = zeros(ng1,1,3);
            grad(:,1,3) = gs;     % first waveform pointing in z
        case 4
            grad = zeros(ng1,1,3);
            grad(:,1,1) = gs;               % first waveform pointing in x
            grad(1:ng_pha,1,3) = G_pha;      % phase encoding along z
            if rewinder(1)==1
                ii = (1:ng_pha)+ng1-ng_pha;
                grad(ii,:,3) = -G_pha;
            end
    end
else
    grad = zeros(ng1,n2,dim-(dim==4));
    switch dim
        case 2
            grad(:,:,1) =  gs*cos(phi);
            grad(:,:,2) = -gs*sin(phi);
        case 3
            grad(:,:,1) = -gs*(sin(theta).*cos(phi));
            grad(:,:,2) =  gs*(sin(theta).*sin(phi));
            grad(:,:,3) =  gs*(cos(theta));
        case 4
            grad(:,:,1) =  repmat(gs*cos(phi),[1 mtx(2)]);
            grad(:,:,2) = -repmat(gs*sin(phi),[1 mtx(2)]);
            grad(1:ng_pha,:,3) = bsxfun(@times,G_pha.',agz);
            if rewinder(1)==1
                ii = (1:ng_pha)+ng1-ng_pha;
                grad(ii,:,3) = -bsxfun(@times,G_pha.',agz);
            end
        otherwise, error('dim(=%g) must be 2,3 or 4',dim);
    end
end

%% bonus spec
if sum(bonus_fids(:))>1 % if bonus spec added, zero gradients for them
    grad(:,end-sum(bonus_fids(:))+1:end,:) = 0;
end

%% print info + plotting
fprintf('size(grad)=(pts/intlv intlv groups)=(%g %g %g)\n',...
    ng1,n2,size(grad,3));
if dim==3, fprintf('theoretical number of spokes=%g\n',Ns); end
dur = ng1*gdt*1d3;
fprintf('grad duration=%g[ms]\n',dur);
if do_rot_file
    figure(101); clf;
    if dim==3, subplot(2,1,1); end
    plot(phi);
    axis tight;
    ylabel('phi'); 
    grid on;
    if dim==3
        subplot(2,1,2);
        plot(theta);
        axis tight;
        ylabel('theta');
        grid on;
    end
else
    y_str = {'Gx','Gy','Gz'};
    figure(101); clf;
    for l3=1:(dim-(dim==4))
        subplot(dim-(dim==4),1,l3); 
        if ~exist('ng_ramp1','var')
            ng_ramp1 = ceil(g_read/smax/gdt)+1;  % #ramp points
        end
        plot(grad(ng_off+ng_ramp1+1,:,l3)); 
        axis tight;
        ylabel(y_str(l3));
    end
end


%% calculate k-space of single spoke
ks = cumsum(gs,1)*gdt*GAM*res(1);      % k-space normalised to +/-0.5; gdt
nks0 = length(ks);                     % #k-pts relative to gdt
nks1 = floor(gdt/kdt*nks0);            % #k-pts relative to kdt
% resampling ks from gdt to kdt
ks = interp1((0:(nks0-1)).'*gdt,ks,(0:(nks1-1)).'*kdt,'spline');
nk_pre = round(ng_pre*gdt/kdt);        % #samples prewinder (typ~=2 -> 0)


% indexing k-space
if max(ks)<0.5
    warning('max(ks)(=%g)<0.5',max(ks));
    if rewinder(1,1)==1      % exclude rewinder
        ii = (nk_off+nk_pre):(nks1-floor(ng_ramp2*gdt/kdt));
    else                % include ramp down
        ii = (nk_off+nk_pre):nks1;
    end
else
    tmp = find(ks>0.5);
    ii = (nk_off+nk_pre):tmp(1);
end
ks = ks(ii,1);

fprintf('max(ks)=%g\n',max(ks));
if ks(end-1,1)>0.5, warning('ks(end-1,1)(=%g)>0.5',ks(end-1,1)); end
if max(ks)<0.5,     warning('max(ks)(=%g)<0.5',max(ks)); end
nks = size(ks,1);                      % #k-pts
ts = (0:(nks-1))*kdt;                  % time of single spoke [s]


%% calculate index for reco
n_ind = ceil((nks+nk_off+nk_pre)/32)*32; % #acq pts rounded to multiple of 32
if n_ind < 64, n_ind = 64; end         % fidal miminum=64
inds = false(1,n_ind);
inds(1,(nk_off+nk_pre-1)+(1:nks)) = true;
ind = repmat(inds,[n2 1]);
t = repmat(ts,[n2 1]);
fprintf('acq: BW=%g[kHz]; #pts=%g; duration=%g[ms]\n',...
    1d-3/kdt,n_ind,n_ind*kdt*1d3);
fprintf('time of spoke: ts=%g[ms]\n',ts(1,end)*1d3);

%% assemble 2D/3D k-space via rotating spoke
k = zeros(nks,n2,dim-(dim==4));
kz = [];
switch dim
    case 2
        k(:,:,1) =  ks*cos(phi);
        k(:,:,2) = -ks*sin(phi);
        
        dcfs = pi/n2*(ks(2:end,1).^2-ks(1:end-1,1).^2)*mtx(1)^2;
        theta = [];
        % bonus spec
        if sum(bonus_fids(:))>1 % if bonus spec added, zero gradients for them
            k(:,end-sum(bonus_fids(:))+1:end,:) = 0;
        end
    case 3
        k(:,:,1) = -ks*(sin(theta).*cos(phi));
        k(:,:,2) =  ks*(sin(theta).*sin(phi));
        k(:,:,3) =  ks*(cos(theta));

        dcfs = 4/3*pi/n2*(ks(2:end).^3-ks(1:end-1).^3)*mtx(1)^3;

        figure(103); clf
        plot3(k(end,:,1),k(end,:,2),k(end,:,3),'.')
    case 4
        k(:,:,1) =  repmat(ks*cos(phi),[1 mtx(2)]);
        k(:,:,2) = -repmat(ks*sin(phi),[1 mtx(2)]);
        kz = agz/2;
        k(:,:,3) = repmat(kz,[nks,1]);
        
        dcfs = pi/length(phi)*(ks(2:end,1).^2-ks(1:end-1,1).^2)*mtx(1)^2;
        theta = [];
end


%% calculate density-compensation function
dcfs = abs(dcfs);
dcfs = (dcfs([1:(nks-1),(nks-1)],1)+dcfs([1,1:(nks-1)],1))/2;   % spoke dir
% dcfs = interp1(linspace(-1,1,nks-1),dcfs,linspace(-1,1,nks),'spline');
if (typ==2)&&(order(1)~=3), dcfs = dcfs/2; end   % correct scaling for n2
% if (typ==2), dcfs = dcfs/2; end   % correct scaling for n2

if order(1)==3          % golden angle
    dcft = sub_get_dcft(phi,theta,dda(1),dda(2));  % distance on circle/area on surface
    if order(4)==1      % averaging out distance; required for QTI
        dcft = repmat(mean(dcft),[1 length(dcft)]);
    end
    if dim==4, dcft = repmat(dcft,[1 mtx(2)]); end
    dcf = dcfs*dcft;
    % if dim==4, dcf = repmat(dcf,[1 mtx(2)]); end
    % bonus spec
    if sum(bonus_fids(:))>1 % if bonus spec added, not for recon
        dcf(:,end-sum(bonus_fids(:))+1:end) = 0;
    end
    dcf = reshape(dcf,[1 nks*n2]);
else
    if dda(1)>0, warning('dda(1)(=%g)>0 && not Golden Angle -> ignoring',dda(1)); end
    dcf = repmat(dcfs,[1 n2]);
    % bonus spec
    if sum(bonus_fids(:))>1 % if bonus spec added, not for recon
        dcf(:,end-sum(bonus_fids(:))+1:end) = 0;
    end
    dcf = reshape(dcf,[1 nks*n2]);
    dcft = [];
end
k = reshape(k,[1 nks*n2 size(k,3)]);


%% plotting
figure(102); clf;
subplot(2,1,1); 
ttg = (((1:ng1)-ng_off-ng_pre-1))*gdt*1d3;
% ii = false(size(ttg));
ii = (ttg>=0) & (ttg<=(ts(1,end)*1d3));
plot(ttg,gs*1d3,'b-',ttg(1,ii),gs(ii,1)*1d3,'bx');
if dim==4
    hold on; plot(ttg,grad(:,1,3)*1d3,'m-'); hold off
end
ylabel('grad [mT/m]'); 
grid on
ttax = [-ng_off-ng_pre,(ng1-ng_off-ng_pre)]*gdt*1d3;
axis([ttax(1) ttax(2) min(grad(:))*1.1d3 max(grad(:))*1.1d3]);

subplot(2,1,2);
plot(ts*1d3,ks,'b-x',ts*1d3,dcfs,'r-');
xlabel('time [ms]');
ylabel('k - dcf');
legend({'k','dcf'});
grid on
axis([ttax(1) ttax(2) min(ks) max([max(dcfs),0.55])]);


%% check gmax + smax
gmax_act = max(abs(grad(:)));
smax_act = max(max(max(abs(diff(grad,[],1)))))/gdt;
if gmax_act>(gmax+1d-10)
    warning('gmax: actual(=%g)>nominal(=%g) [mT/m]',gmax_act*1d3,gmax*1d3);
end
if smax_act>smax
    warning('smax: actual(=%g)>nominal(=%g) [T/m/s]',smax_act,smax);
end


%% misc checks
if sum(ind(:))~=size(k,2)
    warning('sum(ind(:))(=%g)~=size(k,2)(=%g)',sum(ind(:)),size(k,2));
end
if max(dcf)>1
    warning('max(dcf)(=%g)>1',max(dcf));
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
        switch typ
            case 0, typ_str = '_coca'; 
            case 1, typ_str = '_coda'; 
            case 2, typ_str = '_fsca';
        end
        switch rewinder(1,1)
            case 0, rew_str = '';
            case 1, rew_str = ['_rew' num2str(rewinder(1,1))];
            case 2, rew_str = sprintf('_rew%d_%.2f',rewinder(1,1:2));
        end
        if dim==4
            fov_str = sprintf('%g_%g',round(fov(1)*1d3),round(fov(2)*1d3));
            mtx_str = sprintf('%d_%d',mtx(1),mtx(2));
        else
            fov_str = num2str(round(fov(1)*1d3));
            mtx_str = num2str(mtx(1));
        end
        
        fname = sprintf('interleaved_radial%gD_%s_fov%s_mtx%s_intlv%g_kdt%g_gmax%g_smax%g_dur%g%s%s', ...
            dim,nuc,fov_str,mtx_str,n2,kdt*1d6,...
            round(gmax_act*1d3),round(smax_act),round(dur*10)/10,...
            typ_str,rew_str);
        fname = strrep(fname,'.','p');
        switch order(1)
            case 2, fname = [fname '_RND'];
            case 3, fname = [fname '_GA' num2str(order(2))];
            case 4, fname = [fname '_HAL'];
        end
        if ((dim==4) && (order(3)>1))
            fname = [fname '_zPE' num2str(order(3))];
        end
    else
        fname = [];
    end
end


%% saving waveform
if ~isempty(fname)
    fov1H = fov(1)/abs(gyrogamma('1h'))*abs(gyrogamma(nucleus));
    desc = sprintf('design_interleaved_radial(%g,%d,%d,fname,%g,%g,%s,%d,%d,[%g,%g],%d,%d,%d)\n',...
        fov(1)*1d3,mtx(1),kdt*1d6, gmax*1d3,smax,nuc,dim,typ,sampfac(1),...
        sampfac(2),rewinder(1,1),order(1,1),do_rot_file);

    phi = phi*180/pi;                  % convert [rad] -> [deg]
    theta = theta*180/pi;              % convert [rad] -> [deg]
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
    out = write_ak_wav([fname,'.wav'],grad,bw,fov1H,desc,n_ind);
    if isfield(out,'wave'), out = rmfield(out,'wave'); end
    if isfield(out,'grad'), out = rmfield(out,'grad'); end
    if save_reduced
        nviews = n2;
        tmp = zeros(1,nks,3);
        if dim==3
            tmp(1,:,3) = ks.';
        else
            tmp(1,:,1) = ks.';
        end
        ks = tmp;
        save(fname,'out','ks','ts','inds','dcfs','dcft','nviews','fov',...
            'fov1H','mtx','gmax','smax','nucleus','n_ind','npix',...
            'phi','theta','kz','flips','freqs');
    else
        dcf = single(dcf); k = single(k); t = single(t);
        save(fname,'out','k','t','ind','fov','fov1H','mtx','gmax','smax',...
            'nucleus','dcf','n_ind','npix','phi','theta','flips','freqs');
    end
    

    if do_rot_file
        write_fdl([fname '_phi.fdl'],mod(phi,360),'phi');
        if dim==3
            write_fdl([fname '_theta.fdl'],mod(theta,360),'theta');
        end
        if dim==4
            write_fdl([fname '_agz.fdl'],agz,'agrad');
        end
    end
    if var_flips
        flips_all = repmat(flips,1,Ns); % repeat flips for all trajs
        for i=1:intlvs
            flips_all = [flips_all, repmat(flips(i),1,bonus_fids(i))];
        end
        write_fdl([fname '_flip.fdl'],flips_all,'flip');
    end
    if var_freqs
        freqs_all = repmat(freqs,1,Ns); % repeat flips for all trajs
        for i=1:intlvs
            freqs_all = [freqs_all, repmat(freqs(i),1,bonus_fids(i))];
        end
        write_fdl([fname '_freq.fdl'],freqs_all,'freq');
    end
end


%% output arguments
if nargout==1
    k = fname;
end


end      % design_interleaved_radial.m


%% sub-functions
function gs = sub_da_spoke(g_read,smax,kmax,gdt,GAM,dim,ng_off)
% slow down outside to improve SNR
% 3D=mrm62_1565_2009_nagel.pdf
% 2D=mrm65_1091_2011_konstandin.pdf
ng_ramp1 = ceil(g_read/smax/gdt)+1;       % #ramp points
switch dim
    case {2,4}, k0 = GAM*g_read^2/smax;      % Eq9(2D)
    case 3,     k0 = 2*GAM*g_read^2/smax;    % Eq9(3D)
end
if k0>kmax, warning('k0(=%g)>kmax(=%g)',k0,kmax); end
ng0 = ceil(k0/(gdt*GAM*g_read)-0.5*ng_ramp1);
if ng0<0, warning('ng0(=%g)<0'); end
switch dim
    case {2,4} % 2D
        ngc = ceil((kmax^2-k0^2)/(2*GAM*k0*g_read*gdt)); % Eq4(2D)
        gda = k0*g_read./sqrt(2*GAM*k0*g_read*...
            (0:(ngc-1))*gdt+k0^2); % Eq5(2D)
    case 3     % 3D
        ngc = ceil((kmax^3-k0^3)/(3*GAM*k0^2*g_read*gdt)); % Eq6(3D)
        gda = k0^2*g_read*(3*GAM*k0^2*g_read*...
            (0:(ngc-1))*gdt+k0^3).^(-2/3); % Eq7(3D)
end
if ngc>5d4, warning('#waveform/spoke too large: ngc(=%g)>5d4',ngc); end
if ngc>2d5, error('#waveform/spoke too large: ngc(=%g)>2d5',ngc); end
fprintf('Density-adapted radial:\n');
fprintf('ng_ramp_up=%g; ng0=%g; ngc=%g\n',ng_ramp1,ng0,ngc);
fprintf('kmax=%g; k0=%g\n',kmax,k0);

gs = [zeros(ng_off,1) ; ...
    linspace(0,g_read,ng_ramp1).' ; ...
    g_read*ones(ng0,1) ; ...
    gda.'];
ng_tmp = ng_off+ng_ramp1+ng0+ngc;
kact = sum(gs(1:ng_tmp,1))*gdt*GAM;
if kact<kmax,     warning('kact(=%g)<kmax(=%g)',kact,kmax); end
if kact>1.1*kmax, warning('kact(=%g)>1.1*kmax(=%g)',kact,kmax); end

end      % sub_da_spoke


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

%% halton number calculation
function result = haltonnumber(index, base)
    result = 0;
    f = 1.0;
    i = index;
    while i>0
        f = f/base;
        result = result + f*mod(i,base);
        i = fix(i/base);
    end
end

%% density-compensation function along tangential direction
function dcft = sub_get_dcft(phi,theta,dda,nseg)

nr = 5;
if dda>0
    if nseg>1
        nphi = size(phi,2);
        nn = nphi/nseg;
        if abs(nn-round(nn))>1d-10
            warning('size(phi,2)(=%g) not multiple of nseg(=%g)',nphi,nseg);
        end
        iida = [];
        for l=1:nseg
            iida = [iida ((dda+1):(nphi/nseg))+(l-1)*nphi/nseg];
        end
    else
        iida = (dda+1):size(phi,2);
    end
    phi   = phi(1,iida);
    if ~isempty(theta), theta = theta(1,iida); end
end

n2 = size(phi,2);
dcft = zeros(1,n2);

if isempty(theta)       % 2D: distance on circle
    [sphi,ii] = sort(mod(phi,2*pi));
    if all(phi<1.01*pi) % full spoke
        sphi = [sphi,sphi(1,1)+pi];
    else
        sphi = [sphi,sphi(1,1)+2*pi];
    end
    sphi = diff(sphi);
    sphi = (sphi+sphi(1,[n2,1:n2-1]))/2;
    dcft = zeros(1,n2);
    dcft(1,ii) = sphi*n2/(2*pi);
else                    % 3D: area of 3D radial on surface
    % approximation via radius + area of circle
    x = -(sin(theta).*cos(phi));
    y =  (sin(theta).*sin(phi));
    z =  (cos(theta));
    x = single(x); y = single(y); z = single(z);
    rr_thresh = 1/sqrt(n2);
    ni = zeros(1,n2);
    for ll=1:n2
        rr = (x-x(1,ll)).^2+(y-y(1,ll)).^2+(z-z(1,ll)).^2;
        rr2 = rr(1,rr<rr_thresh);      % thresholding to improve speed
        ni(1,ll) = size(rr2,2);
        if (ni(1,ll)-1)<nr
            warning('ni(=%d)-1<nr(=%d): increase rr_thresh',...
                ni(1,ll),nr);
            rr_thresh = rr_thresh*1.25;
            rr2 = rr;
        end
        rr = sort(rr2);
        dcft(1,ll) = mean(sqrt(rr(1,2:(nr+1))),2).^2;
    end
    fprintf('rr_thresh=%g; n2=%g; nr=%g; ni=%g+-%g\n',...
        rr_thresh,n2,nr,mean(ni),std(ni));
    dcft = dcft*n2*0.065;
end

if dda>0
    if nseg>1
        tmp = zeros(1,nphi);
        tmp(iida) = dcft;
        dcft = tmp;
    else
        dcft = [zeros(1,dda) dcft];
    end
end

end      % sub_get_dcft
