function [k,dcf,t,grad,ind] = design_corings(fov,mtx,bw_spec,...
    n_spec,scheme,gmax,smax,fname,sampfac,gdt,nucleus,verb)
%DESIGN_CORINGS  Design trajectory for concentric rings MRSI
% [k,dcf,t,grad,ind] = design_corings(fov,mtx,bw_spec,...
%    n_spec,scheme,gmax,smax,fname,sampfac,gdt,nucleus,verb)
%                                                         unit     default
%       fov  Field of view                                [mm]
%       mtx  Matrix size
%   bw_spec  Spectral bandwidth (full)                    [Hz]     (1000)
%    n_spec  Spectral acquisition points                           (256)
%    scheme  Sampling scheme                                       (1)
%            1=uniform
%            2=density-weighted Hanning
%      gmax  Maximum gradient amplitude                   [mT/m]   (33)
%      smax  Maximum slew rate                            [T/M/s]  (120)
%     fname  filename (if given->write_ak_wav) (true=generate)     false
%   sampfac  Sampling factor: >1=over, <1=under-sampling           (1 1 1)
%            1=circumferential; 2=radial; 3=z-direction
%       gdt  Gradient update time                         [us]     (8)
%   nucleus  Nucleus                                               '1H'
%      verb  1&2->print output; 2=plotting                         (2)
%
%         k  k-space trajectory        [1,#exc,direction]
%       dcf  density compensation function
%         t  time
%      grad  gradient trajectory       [#pts,#exc,direction]
%            directions: x+y+z
%       ind  Data selection index
%
%  12/2020 Rolf Schulte
if (nargin<1), help(mfilename); return; end


%% input
if ~exist('scheme','var'),   scheme = []; end
if isempty(scheme),          scheme = 1; end
if ~exist('bw_spec','var'),  bw_spec = []; end
if isempty(bw_spec),         bw_spec = 1000; end
if ~exist('n_spec','var'),   n_spec = []; end
if isempty(n_spec),          n_spec = 256; end
if ~exist('gmax','var'),     gmax = []; end
if isempty(gmax),            gmax = 33; end
if ~exist('smax','var'),     smax = []; end
if isempty(smax),            smax = 120; end
if ~exist('fname','var'),    fname = []; end
if ~exist('sampfac','var'),  sampfac = []; end
if isempty(sampfac),         sampfac = [1 1];  end
sampfac = sampfac(:).';
switch length(sampfac)
    case 1, sampfac = [sampfac 1 1];
    case 2, sampfac = [sampfac 1];
    case 3
    otherwise, warning('length(sampfac)(=%d)~=1,2 or 3',length(sampfac));
end
if ~exist('gdt','var'),      gdt = []; end
if isempty(gdt),             gdt = 8; end
if ~exist('nucleus','var'),  nucleus = []; end
if isempty(nucleus),         nucleus = '1H'; end
if ~exist('verb','var'),     verb = []; end
if isempty(verb),            verb = 2; end


%% checks + convert units
do_intlv = true;        % allow spectral interleaves
n_intlv = 1;            % #spectral interleaves
fov = 1d-3*fov(:).';    % [m]
if length(fov)>3,       error('length(fov)>3'); end
if length(fov)==1,      fov = [fov fov]; end
mtx = mtx(:).';
if isempty(mtx),        error('isempty(mtx)'); end
if length(mtx)>3,       error('length(mtx)>3'); end
if length(mtx)==1,      mtx = [mtx mtx]; end
if any(size(fov)~=size(mtx)), error('any(size(fov)~=size(mtx))'); end
if any(rem(mtx,1)),     error('mtx not natural number'); end
if mtx(1)~=mtx(2),      error('mtx(1)~=mtx(2)'); end
if fov(1)~=fov(2),      error('fov(1)~=fov(2)'); end
gmax = gmax*1d-3;       % [T/m]
if gdt<4,               error('gdt(=%g)<4[us]',gdt); end
if abs(rem(gdt,4))>1d-10, error('gdt(=%g) not multiple of 4[us]',gdt); end
gdt = gdt*1d-6;         % gradient update time [s]
dt = 2d-6;              % trajectory generation dwell time [s]
                        % resampling later to kdt and gdt
MAX_FILTER_POINTS = 32768;  % max #points (cv1)


%% misc parameters
drtsmx = 1/sqrt(2);                         % smax derate factor grad lobes
bw_max = 500d3;                             % max sampling BW [Hz]
ng_max = 32766;                             % max #grad-pts per excitation
GAM = abs(gyrogamma(nucleus))/2/pi;         % gyromagnetic ratio [Hz/T]
dim = size(mtx,2);                          % spatial dimensions
res = fov./mtx;                             % voxel-size [m]
alpha = [1.61 1.71 1.78];                   % Hanning constants: 1D,2D,3D


%% calculate possible BW and duration of single ring
nk_circ = ceil(pi*mtx(1)*sampfac(1));       % #k-pts per ring (samp w kdt)
bw_spat = bw_spec*nk_circ;                  % desired spatial BW
if verb>0
    fprintf('Desired: nk_circ=%d; bw_spec=%g[Hz]; bw_spat=%g[kHz]\n',...
        nk_circ,bw_spec,bw_spat*1d-3);
end
bw_spat_max = gmax*fov(1)*GAM;              % max possible spatial BW
% bw_spat_max = min([bw_spat_max bw_max]);

if bw_spat>bw_max
    warning('bw_spat(=%g)>bw_max(=%g); capping bw_spat',bw_spat,bw_max);
    bw_spat = bw_max;
else
    if bw_spat>bw_spat_max
        if do_intlv
            fprintf('Attention: bw_spat(=%g)>bw_spat_max(=%g)\n',...
                bw_spat,bw_spat_max);
            n_intlv = ceil(bw_spat/bw_spat_max);
            fprintf('\tusing n_intlv=%d spectral interleaves\n',n_intlv);
            bw_spat = bw_spec*nk_circ/n_intlv;
        else
            warning('bw_spat(=%g)>bw_spat_max(=%g); dereasing bw_spec',...
                bw_spat,bw_spat_max);
            bw_spec = bw_spat_max/nk_circ;  % derate spectral BW
            bw_spat = bw_spat_max;          % derate spatial BW
        end
    end
    bw_spat = 0.5d6/floor(0.5d6/bw_spat);   % round up to multiple of 2us
    nk_circ = floor(bw_spat/bw_spec)*n_intlv;
    bw_spat = bw_spec*nk_circ/n_intlv;
    bw_spat = 0.5d6/floor(0.5d6/bw_spat);   % round up to multiple of 2us
end
bw_spec = bw_spat*n_intlv/nk_circ;          % actual spectral BW [Hz]
% n_spec = ceil(n_spec/n_intlv)*n_intlv;      % round-up to multiple
if verb>0
    fprintf('Actual:  nk_circ=%d; bw_spec=%g[Hz]; bw_spat=%g[kHz]\n',...
        nk_circ,bw_spec,bw_spat*1d-3);
end


%% calculate oversampling factors for dt->kdt-gdt
kdt = 1/bw_spat;                            % sampling dwell time [s]
dki = kdt/dt;                               % oversampling factor dt->kdt
if abs(dki-round(dki))>1d-10, warning('dki(=%g) not natural number',dki); end
dki = round(dki);
dgi = gdt/dt;                               % oversampling factor dt->gdt
if abs(dgi-round(dgi))>1d-10, warning('dgi(=%g) not natural number',dgi); end
dgi = round(dgi);


%% calculate readout gradient strength + checks
g_read = pi*mtx(1)*bw_spec/(fov(1)*GAM*n_intlv); % readout grad [T/m]
if g_read>gmax, error('g_read(=%g)>gmax(=%g)',g_read,gmax); end
ntmp = nk_circ*kdt/gdt*n_spec/n_intlv;
if ntmp>ng_max
    error('#grad points(=%d) exceeding fidall limit(=%d); decrease gdt',...
        ntmp,ng_max);
end
if verb>0
    fprintf('g_read=%.4g[mT/m]; dt=%g[us]; kdt=%g[us]; gdt=%g[us]\n',...
        g_read*1d3,dt*1d6,kdt*1d6,gdt*1d6);
end


%% spatial k-space: radii + (for 3D) phase encodes
switch scheme    
    case 1     % uniform sampling
        n_ring = ceil(sampfac(2)*mtx(1)/2)+1;    % #rings
        rad = linspace(0,0.5,n_ring);            % position of rings
        % pe = (-mtx(3)/2:(mtx(3)/2-1))/mtx(3);  % phase encoding (for 3D)
        apo = [];                                % apodisation function
        if dim>2
            n_zpe = mtx(3);
            kz = (-n_zpe/2:(n_zpe/2-1))/n_zpe; % z positions
            kz = [kz(1,(floor(n_zpe/2)+1):end) ...
                kz(1,floor(n_zpe/2):-1:1)];    % start from centre
            apo3 = ones(1,n_zpe);
        else
            n_zpe = 1;
        end
    case 2    % density-weighted Hanning
        % calculate k-space positions in radial direction
        wk_scale = mtx(1)*sampfac(2);
        kk = linspace(0,0.5,10000);
        wk = (1+cos(2*pi*kk/alpha(1)));          % Hanning filter; Eq 4
        wk = wk/sum(wk)*wk_scale;
        ii = [(dim>1) diff(rem(cumsum(wk),1))<0];% index to wk==natural number
        rad = kk(1,ii);                          % lists of radii
        if max(rad)<0.49
            warning('max(rad)(=%g)<0.49: normalising to 0.5', max(rad)); 
            rad = rad/max(rad)*0.5;
        end
        n_ring = size(rad,2);                    % #rings
        apo = (1+cos(pi*(0:(n_ring-1))/n_ring/alpha(1)))/2; % apodisation

        if dim>2
            % calculate k-space positions in z direction
            wk_scale = mtx(3)*sampfac(3)/2;
            kk = linspace(0,0.5,10000);
            wk = (1+cos(2*pi*kk/alpha(1)));      % Hanning filter; Eq 4
            wk = wk/sum(wk)*wk_scale;
            ii = [(dim>1) diff(rem(cumsum(wk),1))<0];
            kz = kk(1,ii);                     % lists of radii
            if max(kz)<0.49
                warning('max(kz)(=%g)<0.49: normalising to 0.5', max(kz));
                kz = kz/max(kz)*0.5;
            end
            kz = [kz -kz(1,2:end)];        % mirror
            n_zpe = size(kz,2);                % #z-phase encodes
            apo3 = (1+cos(pi*(0:(n_zpe-1))/n_zpe/alpha(1)))/2;
        else
            n_zpe = 1;
        end
    otherwise, error('scheme(=%g) unkown',scheme);
end
if max(rad)>0.501, warning('max(rad)(=%g)>0.501',max(rad)); end
if max(rad)<0.49,  warning('max(rad)(=%g)<0.49', max(rad)); end


%% phase encoding gradient
if dim>2
    ga_phas = 1/(2*res(3)*GAM);
    % [G_phas,n_phas] = gradient_lobe_old(ga_phas,dt,gmax,drtsmx*smax,0,true);
    [G_phas,n_phas] = gradient_lobe(ga_phas,0,gmax,drtsmx*smax,dt,true);
else
    n_phas = 1;
end


%% ramp-up gradient
ga_des1 = -nk_circ*kdt/2/pi*g_read;
% [g_ramp1,ngr1] = gradient_lobe_old(ga_des1,dt,g_read,drtsmx*smax,0,'',true);
% [g_ramp2,ngr2] = gradient_lobe_old(0,dt,g_read,drtsmx*smax,g_read,'',true);
[g_ramp1,ngr1] = gradient_lobe(ga_des1,0,g_read,drtsmx*smax,dt);
[g_ramp2,ngr2] = gradient_lobe(0,[0 g_read],g_read,drtsmx*smax,dt);

n_ramp = ceil(max([ngr1 ngr2 n_phas])/dki)*dki;
g_up = zeros(2,n_ramp);
g_up(1,(1:ngr1)+n_ramp-ngr1) = g_ramp1;
g_up(2,(1:ngr2)+n_ramp-ngr2) = g_ramp2;


%% ramp-down gradient
n_down = ceil(g_read/smax/dt)+1;
g_down = g_read*[zeros(1,n_down) ; linspace(1,0,n_down)];


%% one full ring
n_circ = nk_circ*dki;
x = 2*pi*((1:n_circ)-0.5)/n_circ;
g_ring = [g_up , ...
    g_read*repmat([sin(x);cos(x)],[1 ceil(n_spec/n_intlv)]) , ...
    g_down];
g_ring = [g_ring , zeros(2,8-rem(size(g_ring,2),4))];


%% assemble gradient trajectory + resample to gdt
% size of grad: #pts/exc,#exc,dim
grad = bsxfun(@times,2*rad,permute(g_ring(:,1:dgi:end),[2 3 1]));
if isodd(size(grad,1)), grad = [grad ; zeros(1,size(grad,2),dim)]; end
if n_intlv>1
    % adding spectral interleaves -> rotating waveforms
    gtmp = zeros(size(grad,1),n_ring*n_intlv,dim);
    for li=1:n_intlv
        phi = 2*pi*(li-1)/n_intlv;
        ii = (0:(n_ring-1))*n_intlv+li;
        gtmp(:,ii,1) = grad(:,:,1)*cos(phi)-grad(:,:,2)*sin(phi);
        gtmp(:,ii,2) = grad(:,:,1)*sin(phi)+grad(:,:,2)*cos(phi);
    end
    grad = gtmp;
    clear gtmp
end
if dim>2
    rr = 2*repmat(kz,[n_ring 1]); rr = rr(:).';
    gtmp = zeros(size(grad,1),size(grad,2),dim);
    gtmp(:,:,1:2) = grad;
    gtmp = repmat(gtmp,[1 n_zpe 1]);
    gtmp(1:ceil(n_phas/dgi),:,3) = G_phas(1,1:dgi:end).'*rr;
    grad = gtmp;
    clear gtmp
end


%% calculate k-space trajectory
nk_read = ceil(nk_circ*n_spec/n_intlv);
k_ring = cumsum(g_ring,2)*dt*GAM*res(1);
if verb>2
    fprintf('max(k_ring)=%g %g; min(k_ring)=%g %g\n',...
        max(k_ring,[],2),min(k_ring,[],2));
end
if max(abs(k_ring(:)))>0.501
    warning('max(abs(k_ring(:)))(=%g)>0.501',max(abs(k_ring(:))));
end
if max(abs(k_ring(:)))<0.49
    warning('max(abs(k_ring(:)))(=%g)<0.49',max(abs(k_ring(:))));
end

ii = n_ramp+(0:(nk_read-1))*dki;
kk = bsxfun(@times,2*rad,permute(k_ring(:,ii),[2 3 1]));
if n_intlv>1
    ktmp = zeros(size(kk).*[1 n_intlv 1]);
    for li=1:n_intlv
        phi = 2*pi*(li-1)/n_intlv;
        ii = (0:(n_ring-1))*n_intlv+li;
        ktmp(:,ii,1) = kk(:,:,1)*cos(phi)-kk(:,:,2)*sin(phi);
        ktmp(:,ii,2) = kk(:,:,1)*sin(phi)+kk(:,:,2)*cos(phi);
    end
    kk = ktmp;
    clear ktmp
end
k = reshape(kk,[1 nk_read*n_ring*n_intlv 2]);    % save 2D k-space only

t = (0:(nk_read-1))*kdt;
t = repmat(t,[1 n_ring*n_intlv]);
t1 = (0:(nk_circ-1))*kdt;


%% density compensation function
k1 = reshape(kk((1:(nk_circ/n_intlv)),:,:),[1 nk_circ*n_ring 2]);
dcf1 = voronoi_area(k1*mtx(1));
if ~isempty(apo)   % re-apply apodisation from ring spacing
    dcf1 = dcf1.*reshape(repmat(apo,[nk_circ 1]),[1 n_ring*nk_circ]);
end
dcf = repmat(dcf1,[n_spec 1])/n_spec;
dcf = dcf(:).';


%% index list for acquired data
nexc = n_ring*n_intlv*n_zpe;
nk_ramp = n_ramp/dki;
n_ind = ceil((nk_ramp+nk_read)/32)*32;      % round-up to multiple of 32
if dim>2
    ind = false(1,n_ind);
else
    ind = false(nexc,n_ind);
end
ii = nk_ramp-1+(1:nk_read);
ind(:,ii) = true;


%% print info
gmax_act = max(abs(grad(:)));
smax_act = abs(diff(grad,1)/gdt);
smax_act = max(smax_act(:));
if verb>0
    fprintf('gmax: specified=%.4g; actual=%.4g[mT/m]\n',...
        gmax*1d3,gmax_act*1d3);
    if gmax_act>gmax, warning('gmax_act(=%g)>gmax(=%g)',gmax_act,gmax); end
    fprintf('smax: specified=%.4g; actual=%.4g[T/m/s]\n',...
        smax,smax_act);
    if smax_act>smax, warning('smax_act(=%g)>smax(=%g)',smax_act,smax); end
    fprintf('size(ind) = %g %g\n',size(ind,1),size(ind,2));
end
if size(dcf,2)~=size(k,2)
    warning('size(dcf,2)(=%d))~=size(k,2)(=%d)',...
        size(dcf,2),size(k,2)); 
end
if size(t,2)~=size(k,2)
    warning('(size(t,2)(=%d)~=size(k,2)(=%d)',size(t,2),size(k,2)); 
end
if verb>0
    fprintf('n_ring=%d; n_intlv=%d; n_zpe=%d; nexc=%d\n',...
        n_ring,n_intlv,n_zpe,nexc);
    fprintf('n_ind=%g\n',n_ind);
end
if n_ind>MAX_FILTER_POINTS
    warning('n_ind(=%g)>MAX_FILTER_POINTS(=%g)',n_ind,MAX_FILTER_POINTS);
end


%% plot grad, k-space, dcf
if verb>1
    figure(100); clf
    tg = (0:size(grad,1)-1)*gdt*1d3;
    plot(tg,1d3*squeeze(grad(:,end,:)));
    grid on
    xlabel('time [ms]'); ylabel('grad [mT/m]');
    
    figure(101); clf
    plot(k1(1,:,1),k1(1,:,2),'x-');
    axis square
    axis(0.55*[-1 1 -1 1]); grid on
    
    figure(102); clf
    plot3(t*1d3,k(1,:,1),k(1,:,2),'x');
    
    figure(103); clf
    plot(dcf1,'r-x');
    axis([1 size(dcf1,2) 0 max(dcf1)]);
    grid on
    drawnow;
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
        
        fname = sprintf('corings%dD_%s_fov%g_mtx%g_nexc%d_bw%g_ns%g_nr%d_s%d_gmax%g_smax%g_dur%g', ...
            dim,nuc,round(fov(1)*1d3),mtx(1),nexc,round(bw_spec),n_spec,n_ring,scheme,...
            round(gmax_act*1d3),round(smax_act),...
            round(size(grad,1)*gdt*1d3));
        fname = strrep(fname,'.','p');
    else
        fname = [];
    end
end


%% echo-time
if verb>0
    TE = (nk_ramp-1)*kdt;
    fprintf('TE+=%g[ms]; t_readout=%g[ms]\n',...
        TE*1d3,round(size(grad,1)*gdt*1d4)/10);
end


%% OUTPUT
if ~isempty(fname)
    fprintf('Saving trajectory to\n   %s\n',fname);
    fov1H = max(fov)/abs(gyrogamma('1h'))*abs(gyrogamma(nucleus));
    npix = max(mtx); bw = bw_spat; npts = n_ind;
    esp_dir = 0;   % all axes
    esp_time = round(n_intlv/bw_spec*1d6);
    
    desc = sprintf('design_corings(%g,%g,%g,%g,%g,%g,%g)',...
        fov(1)*1d3,mtx(1),bw_spec,n_spec,scheme,gmax*1d3,smax);
    desc = sprintf('%s,fname,[%g %g],%g,%s,%g)',...
        desc,sampfac,gdt,nuc,verb);
     
    out = write_ak_wav([fname,'.wav'],grad,bw_spat,fov1H,desc,n_ind,...
        esp_time,esp_dir,gdt*1d6);
    save(fname,'k','dcf','t','out',...
        'fov','mtx','scheme','bw','npts','gmax','smax','fname',...
        'nucleus','nexc','gdt','dim',...
        'ind','fov1H','t','desc','n_ind','npix',...
        'gmax_act','smax_act','n_spec','bw_spec','nk_circ','n_ring',...
        'n_intlv','n_zpe','kz',...
        'k1','dcf1','t1','apo','esp_dir','esp_time');
end


end    % main function design_corings.m
