function [grad,ind] = design_cart(fov,mtx,bw,gmax,smax,fname,blncd,...
    nucleus,sampfac,dda,pe_sort,gdt2)
%DESIGN_CART  Design Cartesian readout (1D,2D or 3D)
% [grad,ind] = design_cart(fov,mtx,bw,gmax,smax,fname,blncd,nucleus,...
%    sampfac,dda,pe_sort,gdt2)
%                                                       [unit]    (default)
%        fov  Field of view                             [mm]
%        mtx  Matrix size                               [1]
%         bw  Full sampling bandwidth (+-)              [Hz]      ('max')
%       gmax  Maximum gradient strength                 [mT/m]    (33)
%       smax  Maximum slewrate                          [T/m/s]   (120)
%****  fov,mtx,gmax,smax: 1D=[freq]; 2D=[freq phase]; 3D=[freq,phase,phase]
%             dimensions=max(length(fov,mtx,gmax,smax))
%      fname  Filename (string, or logical)                       (false)
%      blncd  Balance gradient moments (for SSFP)                 (true)
%    nucleus  Nucleus                                             ('1H')
%    sampfac  Sampling factor (>1 -> oversampling along freq)     (2)
%        dda  Dummy scans (acquired)                              (0)
%    pe_sort  Sorting of phase encoding:                          (0)
%             0=sequential, 1=centre-out
%       gdt2  Gradient update time for storing          [us]      (4)
%             must be multiple of 4
%
%       grad  Gradient trajectory                       [T/m/s]
%             (points,1,direction) 
%             direction: 1=x=freq, 2=y=phase1, 3=z=phase2
%        ind  Index for acquired data
%
%  2/2021 Rolf Schulte
if (nargin<1), help(mfilename); return; end


%% input default + checking
fov = fov(:).'; mtx = mtx(:).';
if any(size(fov)~=size(mtx)), error('size(fov)~=size(mtx)'); end
if ~exist('bw','var'),     bw = []; end
if isempty(bw),            bw = 'max'; end
if ~exist('gmax','var'),   gmax = []; end
if isempty(gmax),          gmax = 33; end
if length(gmax)<2,         gmax = [gmax gmax]; end
if length(gmax)<3,         gmax = [gmax gmax(2)]; end
gmax = gmax(:).';
if any(gmax<1), warning('gmax(=%g)<1',gmax(1)); input(''); end
if ~exist('smax','var'),   smax = []; end
if isempty(smax),          smax = 120; end
if length(smax)<2,         smax = [smax smax]; end
if length(smax)<3,         smax = [smax smax(2)]; end
smax = smax(:).';
if ~exist('fname','var'),  fname = []; end
if ~exist('blncd','var'),  blncd = ''; end
if isempty(blncd),         blncd = true; end
if ~exist('nucleus','var'),nucleus = ''; end
if isempty(nucleus),       nucleus = '1H'; end
if ~exist('sampfac','var'),sampfac = []; end
if isempty(sampfac),       sampfac = 2;    end
if length(sampfac)~=1, error('length(sampfac)(=%g)~=1',length(sampfac)); end
if ~exist('dda','var'),    dda = []; end
if isempty(dda),           dda = 0; end
if ~exist('pe_sort','var'),pe_sort = []; end
if isempty(pe_sort),       pe_sort = 0; end
if ~exist('gdt2','var'),   gdt2 = []; end
if isempty(gdt2),          gdt2 = 4; end
if gdt2<4,                 error('gdt2(=%g)<4',gdt2); end
if abs(mod(gdt2,4))>1d-10, error('gdt2 (=%g) not multiple of 4',gdt2); end

gdt = 4e-6;        % gradient update time [s]
bw_max = 500d3;    % maximum sampling BW  [Hz]


%% misc parameters
dim = length(mtx);
mtx = [mtx , ones(1,3-dim)];           % add 1's for grad pre-allocation
GAM = abs(gyrogamma(nucleus))/2/pi;    % gyromagnetic ratio [Hz/T]
gmax = gmax*1d-3;                      % [mT/m] -> [T/m]
fov = fov*1d-3;                        % [mm] -> [m]


%% bandwidth and kdt
if ischar(bw)
    if ~strcmpi(bw,'max')
        warning('bw(=%s)~=''max''; setting to max anyways',bw); 
    end
    bw = gmax(1)*fov(1)*GAM*sampfac(1);
    if bw>bw_max, bw = bw_max; end
end
kdt = 1/bw;                            % round to kdt multiple of gdt
kdt = gdt*ceil(kdt/gdt-1d-10);
bw = 1/kdt;

if bw<5d3,    warning('bw (=%g [Hz]) < 5000',bw); end
if bw>bw_max, warning('bw(=%g) > bw_max(=%g)',bw,bw_max); end
dki = ceil(kdt/gdt-1d-10);             % round up to multiple of gdt
if abs(dki-kdt/gdt)>1d-10
    error('1/BW(=%g[us]) not multiple of gdt(=%g[us])',kdt*1d6,gdt*1d6); 
end
n_goff = 1;                            % #grad pts offset for start of traj


%% strength of single readout gradient
g_read = bw/(fov(1)*GAM)/sampfac(1);   % [T/m]
if g_read>gmax(1)
    warning('g_read(=%g[mT/m])>gmax(=%g[mT/m]; setting g_read=gmax',...
        g_read*1d3,gmax(1)*1d3);
    fprintf('Better to increase kdt (less oversampling)\n');
    g_read = gmax(1);
end


%% frequency encoding gradient along Gx
% parameters
n_plat = round(mtx(1)*dki*sampfac(1)); % #plateau pts for readout gradient
n_ramp = ceil(g_read/smax(1)/gdt);   % #ramp pts
ga_rew = -g_read*n_plat*gdt/2;    % gradient area for (p)rewinder
[G_rew,n_rew]  = gradient_lobe(ga_rew,g_read,gmax(1),smax(1),gdt,true);

if blncd
    % gradient moment nulled
    n_grd = n_plat+2*n_rew;
    gfreq = [G_rew(1,end:-1:1) g_read*ones(1,n_plat) G_rew]; 
else
    % normal ramp down
    n_grd = n_plat+n_rew+n_ramp+1;
    gfreq = [G_rew(1,end:-1:1) g_read*ones(1,n_plat+1) ...
        g_read*((n_ramp-1):-1:0)/n_ramp]; 
end


%% add phase-encoding gradients
if dim>1           % 2D MRI - Gy
    [x1,x2] = ndgrid((-mtx(2)/2:mtx(2)/2-1)/mtx(2)*2,...
        (-mtx(3)/2:mtx(3)/2-1)/mtx(3)*2);
    if pe_sort==1
        rr = x1.^2+x2.^2;
        [~,itmp] = sort(rr(:));
        ind_pe = zeros(mtx(2),mtx(3),2);
        [ind_pe(:,:,1),ind_pe(:,:,2)] = ind2sub([mtx(2),mtx(3)],...
            reshape(itmp,[mtx(2) mtx(3)]));
        x1tmp = x1; x2tmp = x2;
        for l1=1:mtx(2)
            for l2=1:mtx(3)
                x1(l1,l2) = x1tmp(ind_pe(l1,l2,1),ind_pe(l1,l2,2));
                x2(l1,l2) = x2tmp(ind_pe(l1,l2,1),ind_pe(l1,l2,2));
            end
        end
        clear x1tmp x2tmp
    end
    ga_pha1 = mtx(2)/fov(2)/GAM/2;
    [G_pha1,n_pha1] = gradient_lobe(ga_pha1,0,gmax(2),smax(2),gdt,true);
    gpha1 = bsxfun(@times,G_pha1.',x1(:).');
    if n_pha1>n_rew, warning('n_pha1(=%g)>n_rew(=%g)',n_pha1,n_rew); end
end
if dim>2           % 3D MRI - Gz
    ga_pha2 = mtx(3)/fov(3)/GAM/2;
    [G_pha2,n_pha2] = gradient_lobe(ga_pha2,0,gmax(3),smax(3),gdt,true);
    gpha2 = bsxfun(@times,G_pha2.',x2(:).');
    if n_pha2>n_rew, warning('n_pha2(=%g)>n_rew(=%g)',n_pha2,n_rew); end
end


%% assemble full gradient matrix
grad = zeros(n_grd,mtx(2)*mtx(3),dim);
grad(:,:,1) = repmat(gfreq.',[1 mtx(2)*mtx(3)]);
if dim>1
    ii1 = (1:n_pha1);
    ii2 = (1:n_pha1)+n_grd-n_pha1;
    grad(ii1,:,2) =  gpha1;
    if blncd, grad(ii2,:,2) = -gpha1; end
end
if dim>2
    ii1 = (1:n_pha2);
    ii2 = (1:n_pha2)+n_grd-n_pha2;
    grad(ii1,:,3) =  gpha2;
    if blncd, grad(ii2,:,3) = -gpha2; end
end


%% index list 
n_ind = ceil(n_grd/dki);          % #acquired points
if n_ind<64
    fprintf('Warning: n_ind(=%g)<64; setting to 64\n',n_ind);
    n_ind = 64;
end
n_ind = 32*ceil(n_ind/32);        % rounding to multiple of 32 for fidall
ind_grad = false(n_grd,1);        % gradient gdt
ii = n_goff+n_rew+(1:n_plat)-1;   % index
ind_grad(ii,1) = true;
ind = false(n_ind,1);             % acq index list (kdt)
ind(1:ceil(n_grd/dki),:) = ind_grad(1:dki:n_grd,1);


%% k-space trajectory (for recon_spiral.m and recon_grid3d.m)
kk = cumsum(grad,1)*gdt*abs(gyrogamma(nucleus));
k = kk(1:dki:end,:,:);
t = (0:size(k,1)-1)/bw;
k = k(ind,:,:);
k = reshape(k,[1 size(k,1)*size(k,2) size(k,3)]);
k = bsxfun(@rdivide,k,2*pi*reshape(mtx(1,1:dim)./fov,[1 1 dim]));

t = t(ind).';
ind = repmat(ind,[1 mtx(2)*mtx(3)]).';
t = repmat(t,[1 mtx(2)*mtx(3)]);


%% add dummy scans
if dda>0
    ind = [false(dda,size(ind,2)) ; ind];
    grad = [zeros(size(grad,1),dda,dim) , grad];
    nexc = mtx(2)*mtx(3) + dda;
else
    nexc = mtx(2)*mtx(3);
end


%% print info
fprintf('BW=%g [kHz] (full)=+-%g [kHz] (imaging)\n',bw*1d-3,bw*0.5d-3);
fprintf('g_read=%g [mT/m]; n_plat=%g; n_rew=%g; n_ramp=%g\n',...
    g_read*1d3,n_plat,n_rew,n_ramp);
fprintf('gdt=%g [us]; kdt=%g [us]; dki=%g; ',gdt*1d6,kdt*1d6,dki);
fprintf('n_grd=%g; n_ind=%g; sum(ind(:,1)=%g\n',n_grd,n_ind,sum(ind(:,1)));
if length(gfreq)~=n_grd
    warning('length(gfreq)(=%g) ~= ng(=%g)',length(gfreq),n_grd); 
end


%% print info + plotting
fprintf('#points = %g\n',n_grd);
fprintf('size(grad) = %g %g %g\n',size(grad,1),size(grad,2),size(grad,3));
dur = n_grd*gdt*1d3;
fprintf('grad duration = %g [ms]\n',dur);
fprintf('acq. duration = %g [ms]\n',n_ind*kdt*1d3);


%% plotting
clf;
tg = (0:n_grd-1)*gdt*1d3;
nsbp = 3 + (dim>1);
subplot(nsbp,1,1);
plot(tg,gfreq*1d3,'b',tg(ind_grad),gfreq(ind_grad)*1d3,'rx',...
    [0 dur],gmax(1)*[1d3 1d3],'m--',[0 dur],-gmax(1)*[1d3 1d3],'m--');
ylabel('G-freq [mT/m]'); grid on; axis([0 dur gmax(1)*1.1d3*[-1 1]]);
subplot(nsbp,1,2); 
csgf = cumsum(gfreq)*gdt;
plot(tg,csgf,'b'); 
ylabel('sum grad'); grid on; 
axis([0 dur -1.1*max(abs(csgf)) 1.1*max(abs(csgf))]);
subplot(nsbp,1,3); 
plot(tg(1:end-1),diff(gfreq)/gdt,'b',...
    [0 dur],smax(1)*[1 1],'m--',[0 dur],-smax(1)*[1 1],'m--');
ylabel('slewrate [T/m/s]'); grid on; axis([0 dur -1.1*smax(1) 1.1*smax(1)]);
xlabel('time [ms]');
if dim>1
    subplot(nsbp,1,4);
    [~,ii] = max(reshape(grad(:,:,2),[size(grad,1)*size(grad,2) 1]));
    [~,i2] = ind2sub([size(grad,1) size(grad,2)],ii);
    plot(tg,grad(:,i2,2)*1d3,'b',...
        tg(ind_grad),grad(ind_grad,1,2)*1d3,'rx',...
        [0 dur],gmax(2)*[1d3 1d3],'m--',[0 dur],-gmax(2)*[1d3 1d3],'m--');
    ylabel('G-pe [mT/m]'); grid on; axis([0 dur -gmax(2)*1.1d3 gmax(2)*1.1d3]);
end


%% check gmax + smax
gmax_act = max(abs(grad(:)));
sr = diff(grad,[],1)/gdt;
smax_act = max(abs(sr(:)));
if gmax_act>max(gmax)
    warning('gmax: actual(=%g)>nominal(=%g) [mT/m]',gmax_act*1d3,max(gmax)*1d3);
end
if smax_act>max(smax)
    warning('smax: actual(=%g)>nominal(=%g) [T/m/s]',smax_act,max(smax));
end


%% generate filename
if (islogical(fname) || isnumeric(fname))
    if fname
        nuc = nucleus; if isnumeric(nuc), nuc = num2str(nuc); end
        if isnumeric(nucleus), error('Please enter string for nucleus'); end
        if any(diff(fov(1,1:dim))~=0)
            switch dim
                case 2, fov_str = sprintf('%gx%g',round(fov(1:2)*1d3));
                case 3, fov_str = sprintf('%gx%gx%g',round(fov(1:3)*1d3));
            end
        else
            fov_str = sprintf('%g',round(max(fov)*1d3));
        end
        if any(diff(mtx(1,1:dim))~=0)
            switch dim
                case 2, mtx_str = sprintf('%gx%g',mtx(1:2));
                case 3, mtx_str = sprintf('%gx%gx%g',mtx(1:3));
            end
        else
            mtx_str = sprintf('%g',max(mtx));
        end
        
        dur_str = sprintf('%.3g',dur);
        dur_str = regexprep(dur_str,'\.','p');   % replace '.' by 'p'
        fname = sprintf('cart%gD_%s_fov%s_mtx%s_nexc%d_kdt%g_gmax%g_smax%g_dur%s', ...
            dim,nuc,fov_str,mtx_str,nexc,kdt*1d6,round(gmax_act*1d3),...
            round(smax_act),dur_str);
        if blncd, fname = [fname '_bal']; end
        if pe_sort==1, fname = [fname '_co']; end
    else
        fname = [];
    end
end


%% saving waveform
if ~isempty(fname)
    fov1H = max(fov)/abs(gyrogamma('1h'))*abs(gyrogamma(nucleus));
    desc = ['Cartesian readout; design_cart.m\nfname=' fname];
    npix = max(mtx);
    fprintf('Saving gradient + mat file to\n\t%s\n',fname);
    if gdt2>4
        dgi = round(gdt2/(gdt*1d6));
        grad2 = grad(1:dgi:end,:,:);
        if any(any(abs(grad2(end,:,:))>1d-10))
            grad2 = [grad2;zeros(1,size(grad2,2),dim)];
        end
        if isodd(size(grad2,1))
            grad2 = [grad2;zeros(1,size(grad2,2),dim)];
        end
    else
        grad2 = grad;
    end
    
    out = write_ak_wav([fname,'.wav'],grad2,bw,fov1H,desc,n_ind,[],[],gdt2);
    if dim>2       % reduce storage space for 3D
        k = single(k); t = single(t);
    end
    dcf = ones(1,prod(mtx)*sampfac(1),class(k))/sampfac(1);
    if ~exist('ind_pe','var'), ind_pe = []; end
    save(fname,'out','ind','fov','fov1H','mtx','gmax','smax','nucleus',...
        'npix','n_ind','t','k','dcf','sampfac','dim','dda','ind_pe');
end

end      % design_cart.m
