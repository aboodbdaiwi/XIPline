function [ksca,dcf,t,grad,out,fname] = design_mrsi(fov,mtx,scheme,ktraj,bw,...
    npts,gmax,smax,fname,samp,nucleus,npre,verb,crush)
%DESIGN_MRSI  Design trajectory for MRSI with arbitrary sampling patterns
% [k,dcf,t,grad,out,fname] = design_mrsi(fov,mtx,scheme,ktraj,bw,npts,...
%                       gmax,smax,fname,samp,nucleus,npre,verb,crusher)
%                                                         unit     default
%     fov  Field of view                                [mm]
%          anisotropic FOV: reference: fov(1)=fov in freq encode dir
%     mtx  Matrix size
%  scheme  Sampling scheme                                       (1)
%          1=uniform rectangular
%          2=uniform circular
%          3=accumulation-weighted Hanning
%          4=density-weighted Hanning
%   ktraj  Spatial k-space trajectory                            (1)
%          1=consecutively
%          2=random
%          3=in-out
%      bw  Spectral bandwidth (full)                    [Hz]     (5000)
%    npts  Spectral acquisition points                           (1024)
%    gmax  Maximum gradient amplitude                   [mT/m]   (33)
%    smax  Maximum slew rate                            [T/M/s]  (120)
%   fname  filename (if given->write_ak_wav) (true=generate)     false
%    samp  Sampling factor: 
%          scheme=3: #averages at k=0                            (6)
%          scheme=4: >1=over, <1=under-sampling; 
%                    directions=[radial,circumferential]         ([1,1])
% nucleus  Nucleus                                               '1H'
%    npre  #acq points before phase encoding gradient            (0)
%    verb  1 -> print output + plot gradients & k-space          (1)
%          2 -> also run analyse_trajectory.m
%   crush  Crusher: 0=off                                        (0)
%                   1=same dir as phase encode;crush(2)=factor kmax
%                   2=refocusing (balanced)
%
%       k  k-space trajectory        [1,#exc,direction]
%     dcf  density compensation function
%       t  time
%    grad  gradient trajectory       [#pts,#exc,direction]
%          directions: x+y+z
%     out  output structure of write_ak_wav.m
%
%  Literature: mrm50_1266_2003_greiser.pdf
%  12/2020 Rolf Schulte
if (nargin<1), help(mfilename); return; end


%% input
fov = 1d-3*fov(:).';       % [m]
if isempty(fov),           error('isempty(fov)'); end
if length(fov)>3,          error('length(fov)>3'); end
if max(fov)>fov(1),        warning('max(fov)(=%g)>fov(1)(=%g)',max(fov),fov(1)); end
mtx = mtx(:).';
if isempty(mtx),           error('isempty(mtx)'); end
if length(mtx)>3,          error('length(mtx)>3'); end
if any(size(fov)~=size(mtx)), error('any(size(fov)~=size(mtx))'); end
if any(rem(mtx,1)),        error('mtx not natural number'); end
if ~exist('scheme','var'), scheme = []; end
if isempty(scheme),        scheme = 1; end
if ~exist('ktraj','var'),  ktraj = []; end
if isempty(ktraj),         ktraj = 1; end
if ~exist('bw','var'),     bw = []; end
if isempty(bw),            bw = 5000; end
if ~exist('npts','var'),   npts = []; end
if isempty(npts),          npts = 1024; end
if ~exist('gmax','var'),   gmax = []; end
if isempty(gmax),          gmax = 33; end
gmax = gmax*1d-3;          % [T/m]
if ~exist('smax','var'),   smax = []; end
if isempty(smax),          smax = 120; end
if ~exist('fname','var'),  fname = []; end
if ~exist('samp','var'),   samp = []; end
if isempty(samp)
    if scheme==3
        samp = 6;
    else
        samp = [1 1];
    end
end
if scheme~=3 && length(samp)==1, samp = [samp 1]; end
nak0 = samp(1);              % averages at k=0 for accumulation-weighting
if ~exist('nucleus','var'),  nucleus = []; end
if isempty(nucleus),         nucleus = '1H'; end
if ~exist('npre','var'),     npre = []; end
if isempty(npre),            npre = 0; end
if ~exist('verb','var'),     verb = []; end
if isempty(verb),            verb = 1; end
if ~exist('crush','var'),    crush = []; end
if isempty(crush),           crush = 0; end


%% misc parameters
MAX_FILTER_POINTS = 32768;   % max #points (cv1)
gdt = 4e-6;                  % gradient update time [s]
dim = size(mtx,2);           % spatial dimensions
res = fov./mtx;              % voxel-size [m]
apotune = [];                % fine tuning apodisation

alpha = [1.61  1.71  1.78];
beta = [1.25 1.47 1.73];
beta = 0.93*beta;  % fudge factor: better match nominal and real resolution


%% Cartesian grid
for l=1:dim, x{l} = (-mtx(l)/2:mtx(l)/2-1); end

xx = zeros(1,prod(mtx),dim);
switch dim
    case 1, xx = x{l};
    case 2
        [x1,x2] = ndgrid(x{1},x{2});
        xx(1,:,1) = x1(:).';
        xx(1,:,2) = x2(:).';
    case 3
        [x1,x2,x3] = ndgrid(x{1},x{2},x{3});
        xx(1,:,1) = x1(:).';
        xx(1,:,2) = x2(:).';
        xx(1,:,3) = x3(:).';
    otherwise, error('dim(=%g)~=1,2 or 3');
end


%% spatial k-space
switch scheme    
    case {1,2}     % uniform sampling: 1=rectangular, 2=circular
        k = bsxfun(@rdivide,xx,permute(mtx,[1 3 2]));
        if scheme==2
            rad = sqrt(sum(k.^2,3));        % radius in k-space
            ii = rad<=0.501;                % index inside circle
            k = k(:,ii,:);
            grd_drt = 1;
        else
            grd_drt = sqrt(dim);
        end
        dcf = ones(1,size(k,2));
        apo = dcf;
        
    case 3         % accumulation-weighted Hanning
        nexc_theo = prod(mtx)*nak0/beta(dim);
        wk = (1+cos(2*pi*x{1}.'/mtx(1)/alpha(dim)))/2;
        for l=2:dim
            xrshp = reshape(x{l}/mtx(l),[ones(1,l-1) mtx(l)]);
            wk = bsxfun(@times,wk,(1+cos(2*pi*xrshp/alpha(dim)))/2);
        end
        % wk = beta(dim)*wk;
        aa = wk(:).';
        wk = round(wk(:).'*nak0);
        
        kk = bsxfun(@rdivide,xx,permute(mtx,[1 3 2]));
        k = []; apo = []; dcf = [];
        rad = sqrt(sum(kk.^2,3));      % radius in k-space
        for l=1:nak0
            ii = ((rad<=0.501)&(wk>=l));
            k = [k kk(:,ii,:)];
            apo = [apo aa(1,ii)];
            dcf = [dcf 1./wk(1,ii)];
        end
        grd_drt = 1;
        apotune = nak0*apo.*dcf;       % fine tuning apodisation
        
    case 4         % density-weighted Hanning
        [k,dcf] = sub_dw(mtx,dim,samp,alpha(dim),beta(dim),verb);
        nexc_theo = prod(mtx)*2*beta(dim);
        grd_drt = 1;
        apo = (1+cos(2*pi*k/alpha(dim)))/2;
        
    case 5         % Golden angle - Hanning
        zhalf = false;
        if dim==3
            if ((mtx(1)==mtx(2)) && (mtx(1)==2*mtx(3)))
                zhalf = true;
            end
        end
        if any(diff(mtx)~=0) && ~zhalf
            warning('density-weighting sub-optimal for varying mtx');
        end
        Ns = round(prod(mtx)*2*beta(dim)*samp(1));
        if dim==3, Ns = round(Ns/4); end
        x = (0:(Ns-1));
        rad = (x/Ns).^(1/dim)/2;
        rad = rad./sqrt(1+cos(pi*x/alpha(dim)/Ns));
        switch dim
            case 2
                % ga = pi*(3-sqrt(5));           % Golden angle [rad]
                ga = pi/((sqrt(5)+1)/2);         % Golden angle [rad]
                kk = rad.*exp(1i*ga*x)*samp(2);
                k = zeros(1,Ns,2);
                k(1,:,1) = real(kk);
                k(1,:,2) = imag(kk);
                dcf = prod(mtx)*voronoi_area(k);
            case 3
                ga = [0.4656 0.6823];
                phi = 2*pi*mod(x*ga(2),1);
                if zhalf
                    theta = 2*acos(mod(x*ga(1),1));
                else
                    theta = acos(mod(x*ga(1),1));
                    % mirror to lower hemisphere
                    theta = [theta ; theta+pi]; theta = theta(:).';
                    phi = repmat(phi,[2 1]); phi = phi(:).';
                    rad = repmat(rad,[2 1]); rad = rad(:).';
                end
                k = bsxfun(@times,rad,...
                    [-sin(theta).*cos(phi) ; sin(theta).*sin(phi) ; cos(theta)]);
                k = permute(k,[3 2 1])*samp(2);
                dcf = calc_sdc3(k,mtx(1));
            otherwise
                error('dim(=%d)~=2 or 3',dim);
        end
        apo = 1./dcf;
        grd_drt = 1;
        
    case 6         % 3D density-weighted Hanning: 2D circle + 1D line
        if dim~=3, error('dim(=%d)~=3',dim); end
        if length(samp)<3, samp = [samp(:).' 1]; end
        [k2d,dcf2d] = sub_dw(mtx(1,1:2),2,samp(2:3),alpha(2),beta(2),verb);
        [k1d,dcf1d] = sub_dw(mtx(1,3),1,samp(1),alpha(1),beta(1),verb);
        if verb>0
            figure(13); clf
            subplot(2,1,1); 
            plot(k1d,dcf1d,'bx',0.5,1,'r*'); grid on
            axis([1.05*min([k1d -0.5]) 1.05*max([k1d 0.5]) 0 1.05*max([dcf1d 1])])
            subplot(2,1,2);
            k2dabs = sqrt(sum(k2d.^2,3));
            plot(k2dabs,dcf2d,'bx',0.5,1,'r*'); grid on
            axis([0 1.05*max(k2dabs) 0 1.05*max([dcf2d(:);1])]);
        end
        n2d = size(k2d,2); n1d = size(k1d,2);
        k = zeros(1,n1d*n2d,3);
        k(1,:,1:2) = repmat(k2d,[1 n1d 1]);
        k1d = repmat(k1d,[n2d 1]); k1d = k1d(:).';
        k(1,:,3) = k1d;
        dcf = dcf2d.'*dcf1d; dcf = dcf(:).';
        
        nexc_theo = prod(mtx)*2*beta(dim);
        grd_drt = sqrt(2);
        apo = (1+cos(2*pi*k/alpha(dim)))/2;
        
    case 7         % 3D density-weighted Hanning: 2D circle + 1D line
        if dim~=3, error('dim(=%d)~=3',dim); end
        if length(samp)<3, samp = [samp(:).' 1]; end
        k = zeros(1,0,3); dcf = zeros(1,0);
        [k1d,dcf1d] = sub_dw(mtx(1,3),1,samp(1),alpha(1),beta(1),verb);
        k1max = 1.1*max(abs(k1d));
        for l=1:size(k1d,2)
            ss = sqrt(k1max.^2-abs(k1d(1,l)).^2);
            [k2d,dcf2d] = sub_dw(mtx(1,1:2),2,[2*ss*samp(2) samp(3)],alpha(2),beta(2),verb);
            ii = size(k,2)+(1:size(k2d,2));
            fprintf('%d %g %g\n',l,ss,length(ii))
            k(1,ii,1:2) = k2d*ss*2;
            k(1,ii,3) = k1d(1,l);
            dcf(1,ii) = dcf2d*dcf1d(l);
        end

        if verb>0
            figure(13); clf
            subplot(2,1,1); 
            plot(k1d,dcf1d,'bx',0.5,1,'r*'); grid on
            subplot(2,1,2); 
            plot(sqrt(sum(k2d.^2,3)),dcf2d,'bx',0.5,1,'r*'); grid on
        end
        
        nexc_theo = prod(mtx)*2*beta(dim);
        grd_drt = sqrt(2);
        apo = (1+cos(2*pi*k/alpha(dim)))/2;
 
    otherwise, error('scheme(=%g) unkown',scheme);
end
nexc = size(k,2);
switch dim
    case 1, nexc_str = '';
    case 2, nexc_str = sprintf('; Cart. circle=%g',pi/4*prod(mtx));
    case 3, nexc_str = sprintf('; Cart. sphere=%g',pi/6*prod(mtx));
end
if exist('nexc_theo','var')
    nexc_str = sprintf('%s; theo = %g',nexc_str,nexc_theo); 
end
fprintf('nexc: actual=%g; Cartesian=%g%s\n',nexc,prod(mtx),nexc_str);
fprintf('mean(dcf)=%g; max(dcf)=%g; std(dcf)=%g\n',...
    mean(dcf),max(dcf),std(dcf));
if mean(dcf)>1, fprintf('mean(dcf)(=%g)>1\n',mean(dcf)); end


%% spatial k-space trajectory
switch ktraj
    case 1         % consecutive: nothing to be done
    case 2         % randomise
        ii = randperm(size(k,2));
        k = k(:,ii,:);
        dcf = dcf(:,ii,:);
    case 3
        [~,ii] = sort(sqrt(sum(k.^2,3)));
        k = k(:,ii,:);
        dcf = dcf(:,ii,:);
    otherwise, warning('ktraj(=%g) unknown',ktraj);
end


%% scale k-space for anisotropic FOV
if any(abs(diff(res))>1d-10)
    fprintf('Scaling k-space for anisotropic resolution\n');
    scale = shiftdim(res(1)./res,-1);
    ksca = bsxfun(@times,k,scale);
else
    ksca = k;
end
kabs = sqrt(sum(ksca.^2,3));
kmax = max(kabs);
dcf05 = max(dcf.*(kabs<0.51));
if dcf05>1
    tmp_str = 'undersampling';
else
    tmp_str = 'oversampling';
end
fprintf('max(dcf) at kabs=0.5: %g -> %s\n',dcf05,tmp_str);
fprintf('kmax = %g\n',kmax);


%% phase encoding gradient
% gmax_phas = gmax/grd_drt/beta(dim)*1.1;      % derate more because kmax>0.5
% smax_phas = smax/grd_drt/beta(dim)*1.1;
gmax_phas = gmax/grd_drt/kmax/2;
smax_phas = smax/grd_drt/kmax/2;
ga_phas = pi/min(res)/gyrogamma(nucleus);
[G_phas,ng_phas] = gradient_lobe(ga_phas,0,gmax_phas,smax_phas,gdt,true);
G_phas = G_phas*min(res)/res(1);
ng_pre = ceil(npre/gdt/bw);


%% crusher gradient
switch crush(1)
    case 0
        ng_crush = 0;
        ng_acq = 0;
    case 1
        if length(crush)==1, crush(2) = 1.5; end
        [G_crush,ng_crush] = gradient_lobe(crush(2)*ga_phas,0,gmax/grd_drt,...
            smax/grd_drt,gdt,true);

        ng_acq = ceil(npts/bw/gdt);
    case 2
        G_crush = -G_phas;
        ng_crush = ng_phas;
        ng_acq = ceil(npts/bw/gdt);
    otherwise
        ng_crush = 0;
        ng_acq = 0;
        warning('crush(1)=%d unknown; skipping',crush(1));
end


%% full gradient waveforms
ng_tot = ceil((ng_phas+ng_pre+ng_acq+ng_crush)/2)*2;
if ng_tot>MAX_FILTER_POINTS
    warning('ng_tot(=%d)>MAX_FILTER_POINTS(=%d)',ng_tot,MAX_FILTER_POINTS); 
end
gg = zeros(ng_tot,1);
gg(ng_pre+(1:ng_phas),1) = 2*G_phas.';
% grad = bsxfun(@times,gg,k);
grad = bsxfun(@times,gg,ksca);
switch crush(1)
    case 1
        kn = bsxfun(@rdivide,k,sqrt(sum(k.^2,3)));     % normalised
        kn(isinf(kn)) = 0.5; kn(isnan(kn)) = 0.5;
        ii = (1:ng_crush)+ng_tot-ng_crush;
        grad(ii,:,:) = bsxfun(@times,G_crush.',kn);
    case 2
        ii = (1:ng_crush)+ng_tot-ng_crush;
        % grad(ii,:,:) = bsxfun(@times,2*G_crush.',k);
        grad(ii,:,:) = bsxfun(@times,2*G_crush.',ksca);
end


%% indexing
read_start = npre+ceil((ng_phas+1)*gdt*bw);  % start of spectral readout
nind = ceil((read_start+npts)/32)*32;       % round-up to multiple of 32
ind = false(1,nind);                        % initialise index list
ind(1,read_start+(1:npts)) = true;          % index list for MRSI
t = (0:nind-1)/bw;                          % time [s]
t = t(1,ind);
fprintf('Rounding #acq to multiple of 32; last %d samples not considered\n',...
    nind-(read_start+npts));
fprintf('\tnind = %d (from %d)\n',nind,read_start+npts);


%% print info + plot grad & k-space
gmax_act = max(abs(grad(:)));
smax_act = abs(diff(grad,[],1)/gdt);
smax_act = max(smax_act(:));
if gmax_act>gmax, warning('gmax_act(=%g)>gmax(=%g)',gmax_act*1d3,gmax*1d3); end
if smax_act>smax, warning('smax_act(=%g)>smax(=%g)',smax_act,smax); end
if verb>0
    fprintf('max grad = %g [mT/m]; specified gmax = %g [mT/m]\n',...
        gmax_act*1d3,gmax*1d3);
    fprintf('max slew = %g [mT/m]; specified smax = %g [mT/m]\n',...
        smax_act,smax);
    fprintf('Duration: phase encoding = %g [ms]; total = %g [ms]\n',...
        ng_phas*gdt*1d3,nind/bw*1d3);
    fprintf('size(ind) = %g %g\n',size(ind,1),size(ind,2));

    % plotting
    figure(11); clf;
    switch dim
        case 1, plot(k,zeros(size(k)),'x');
        case 2
            % plot(k(1,:,1),k(1,:,2),'x');
            plot(ksca(1,:,1),ksca(1,:,2),'x');
            axis square
            axis(1.1*max(ksca(:))*[-1 1 -1 1]); grid on
        case 3
            subplot(2,2,1);
            plot3(k(1,:,1),k(1,:,2),k(1,:,3),'x'); grid on;
            subplot(2,2,2);
            plot(k(1,:,1),k(1,:,2),'x'); grid on;
            subplot(2,2,3);
            plot(k(1,:,1),k(1,:,3),'x'); grid on;
            subplot(2,2,4);
            plot(k(1,:,2),k(1,:,3),'x'); grid on;
    end
    figure(12); clf
    title('density compensation function');
    subplot(2,1,1);
    plot(dcf,'b-x');
    axis([1 size(dcf,2) 0 max([1.05*max(dcf) 1])]);
    grid on
    xlabel('encoding'); ylabel('dcf');
    subplot(2,1,2);
    plot(kabs,dcf,'bx',0.5,1,'r*');
    axis([0 max(kabs) 0 max([1.05*max(dcf) 1])]);
    grid on
    xlabel('|k|'); ylabel('dcf');
    drawnow;
end


%% run analyse_trajectory
if verb>1
    wf.k = k; wf.mtx = mtx; wf.dcf = dcf; wf.fov = fov; wf.t = t;
    wf.scheme = scheme;
    analyse_trajectory(wf,0,'',4);
end


%% generate filename
if (islogical(fname) || isnumeric(fname))
    if fname
        if any(diff(fov))
            fov_str = sprintf('%g_%g',round(fov(1)*1d3),round(fov(2)*1d3));
            if dim>2
                fov_str = sprintf('%s_%g',fov_str,round(fov(3)*1d3)); 
            end
        else
            fov_str = sprintf('%g',round(fov(1)*1d3));
        end
        if any(diff(mtx))
            mtx_str = sprintf('%d_%d',mtx(1),mtx(2));
            if dim>2
                mtx_str = sprintf('%s_%g',mtx_str,mtx(3)); 
            end
        else
            mtx_str = sprintf('%d',mtx(1));
        end
        if isnumeric(nucleus), error('Please enter string for nucleus'); end
        fname = sprintf('mrsi%gD_%s_scheme%g_fov%s_mtx%s_bw%g_npts%g_nexc%g_gmax%g_smax%g',...
            dim,nucleus,scheme,fov_str,mtx_str,bw,npts,...
            nexc,round(gmax_act*1d3),round(smax_act));
    else
        fname = [];
    end
end


%% echo-time
if verb>0
    TE = read_start/bw;
    fprintf('TE += %g [ms]\n',TE*1d3);
end


%% OUTPUT
if ~isempty(fname)
    fprintf('Saving trajectory to\n   %s\n',fname);
    fov1H = fov(1)/abs(gyrogamma('1h'))*abs(gyrogamma(nucleus));
    n_ind = size(ind,2); npix = max(mtx);
    
    desc = sprintf('MRSI trajectories; dim=%g, fov=%s, mtx=%s, nexc=%g\n',...
        dim,fov_str,mtx_str,nexc);
    desc = sprintf('%sgmax=%g, smax=%g, nucleus=%s\n',...
        desc,gmax,smax,nucleus);
    desc = sprintf('%sgdt=%g, fov1H=%g, bw=%g, n_ind=%g\n',...
        desc,gdt,fov1H,bw,n_ind);
    
    out = write_ak_wav([fname,'.wav'],grad,bw,fov1H,desc,n_ind);
    save(fname,'k','dcf','t','out',...
        'fov','mtx','scheme','bw','npts','gmax','smax','fname',...
         'nucleus','npre','nexc','gdt','dim',...
         'ind','fov1H','dcf','t','desc','n_ind','npix',...
         'gmax_act','smax_act','apo','apotune');
else
    out = [];
end


end    % main function design_mrsi.m


%% sub-functions
function [k,dcf] = sub_dw(mtx,dim,samp,alpha,beta,verb)
if any(diff(mtx)~=0)
    warning('density-weighting sub-optimal for varying mtx');
end
wk_scale = prod(mtx)^(1/dim)*samp(1);


kk = linspace(0,0.5,10000);
wk = (1+cos(2*pi*kk/alpha));             % Hanning filter; Eq 4
wk = wk/sum(wk)*wk_scale;
ii = [(dim>1) diff(rem(cumsum(wk),1))<0];% index to wk==natural number
rad = beta*kk(1,ii);                     % lists of radii
if size(rad,2)<2, warning('rad'); disp(rad); end

switch dim
    case 1      % 1D density-weighted
        k = [-rad(1,end:-1:1) 0 rad];
        dcf = diff(k)*mtx(1);
        dcf = interp1(linspace(-1,1,size(k,2)-1),dcf,...
            linspace(-1,1,size(k,2)));
    case 2      % 2D density-weighted
        k = [0,0];
        dcf = rad(1,2)^2*prod(mtx);
        for lk=2:size(rad,2)
            Nphi = round(2*pi*rad(1,lk)/(rad(1,lk)-rad(1,lk-1))*samp(2));
            % Nphi = ceil(Nphi/2)*2;
            % b =  (res(1)/res(2));
            % if b>1, warning('b(=%g)>1',b); end
            % tmp = cumsum(sqrt(1-(1-b^2)*sin(linspace(0,pi,Nphi/2))));
            % tmp = tmp/max(tmp);
            % phi = pi*[tmp,1+tmp].';
            phi = 2*pi*(0:Nphi-1).'/Nphi;
            if verb>3
                fprintf('lk=%g; Nphi=%g; rad=%g\n',lk,Nphi,rad(1,l));
            end
            
            if ~isodd(lk), phi = phi+pi/Nphi; end
            % phi = phi + randn(1,1)/Nphi;
            k = [k; rad(1,lk)*[cos(phi),-sin(phi)]];
            ddcf = (rad(1,lk)-rad(1,lk-1))^2*prod(mtx);
            dcf = [dcf repmat(ddcf,[1 Nphi])];
        end
        dcf = dcf/samp(2);
        k = permute(k,[3 1 2]);
    case 3      % 3D density-weighted
        k = [0,0,0];
        dcf = rad(1,2)^3*prod(mtx);
        for lk=2:size(rad,2)
            d = (rad(1,lk)-rad(1,lk-1))/rad(1,lk);
            d = d/samp(2);
            N = 4*pi/d.^2;
            Ntheta = round(pi/d);
            Dtheta = pi/Ntheta;
            Dphi = d.^2/Dtheta;
            
            if verb>3
                fprintf('lk=%g; d=%g; N=%g; Ntheta=%g; rad=%g; Dtheta=%g; Dphi=%g\n',...
                    lk,d,N,Ntheta,rad(1,lk),Dtheta,Dphi);
            end
            
            for lth=(0:Ntheta-1)
                theta = pi*(lth+0.5)/Ntheta;
                %x = (lth+0.5)/Ntheta;
                %theta = acos(1-2*x);
                
                Nphi = round(2*pi*sin(theta)/Dphi);
                
                phi = 2*pi*(0:Nphi-1).'/Nphi;
                
                knew = rad(1,lk)*[-sin(theta)*cos(phi),...
                    sin(theta)*sin(phi),...
                    repmat(cos(theta),[Nphi 1])];
                
                k = [k; knew];
                ddcf = (rad(1,lk)-rad(1,lk-1))^3*prod(mtx);
                dcf = [dcf repmat(ddcf,[1 Nphi])];
            end
        end
        dcf = dcf/samp(2).^2;
        k = permute(k,[3 1 2]);
end

end      % sub_dw

