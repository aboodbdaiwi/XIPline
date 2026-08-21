function [k,dcf,t,grad,ind,out]=design_spiral_bni(fov,mtx,tacq,vardens,fname,gmax,smax,...
    samp,nucleus,balanced,calc_dcf,plt,gdt2,order,do_rot_file,bonus_spec)
%DESIGN_SPIRAL_BNI  Design a 2D spiral trajectory from BNI code
%[k,dcf,t,grad,ind,out] = design_spiral_bni(fov,mtx,tacq,vardens,fname,gmax,smax,...
%                 samp,nucleus,balanced,calc_dcf,plt,gdt,order,do_rot_file)
%
%        fov  Field of view                                  [mm]
%        mtx  #pixels (Cartesian resolution after gridding)
%       tacq  Duration of single spiral arm                  [ms]       (5)
%    vardens  Variable density (ustype,us_0,us_1,us_r)            (0,1,1,1)
%      fname  Filename: if given->write_ak_wav 
%             if logical true: generate name                        (false)
%       gmax  Max gradient amp                               [mT/m]    (33)
%       smax  Max slew rate                                  [T/m/s]  (120)
%       samp  Sampling factor: >1 over, <1 under (freq,arms)      ([1.2,1])
%    nucleus  Nucleus                                                ('1H')
%   balanced  Balancing gradient area                               (false)
%   calc_dcf  Calculate density compensation function                (true)
%        plt  Plotting                                               (true)
%        gdt  Gradient update time                           [us]       (4)
%      order  Trajectory order:                                         (1)
%             0=consecutive,1=random,2=ll:order(2):end,3=ga with order(n)
%do_rot_file  Write single static spiral into gradient 
%             waveform file & rotate via vap_phiXX.fdl              (false)
% bonus_spec  Number of bonus spect to acquire                          (0)
%
%          k  k-space trajectory [1,#samples*#intlv,2] - [-0.5..0.5]
%        dcf  Density compensation function [1,#samples*#intlv]
%          t  Time               [#intlv,#samples]
%       grad  Gradient waveform [T/m] with dt=4us
%             [#pts/intlv,#intlv,3]
%        ind  Index (2nd dim) for k-space points on spiral (excl ramps,etc)
%        out  Output structure of wrt_wavs
%
% Litrature: mrm66_39_2011_pipe.pdf
% 12/2022  Rolf Schulte
if (nargin<1), help(mfilename); return; end


%% fixed parameters
fufa_gmax = 0.99;
fufa_smax = 0.99;
nk_off = 2;                     % skip first acq points (artefact)
RFS_MAX_NSCANS = 16382;         % maximum number of scans
gdt = 4;                        % gradient update time 4us (compiled into bnispiralgen.c)
GR = (1+sqrt(5))/2;             % golden ratio


%% input parameters
if ~exist('tacq','var'),        tacq = []; end   % duration of single arm
if isempty(tacq),               tacq = 5; end
if ~exist('fname','var'),       fname = []; end  % for waveform file
if ~exist('vardens','var'),     vardens = [0,1,1,1]; end
if ~exist('gmax','var'),        gmax = []; end
if isempty(gmax),               gmax = 33; end   % [mT/m]
if ~exist('smax','var'),        smax = []; end
if isempty(smax),               smax = 120; end
if ~exist('samp','var'),        samp = []; end
if isempty(samp),               samp = [1.2 1]; end  % [freq arms]
if ~exist('nucleus','var'),     nucleus = []; end
if isempty(nucleus),            nucleus = '1H'; end
if ~exist('balanced','var'),    balanced = []; end
if isempty(balanced),           balanced = false; end
if ~exist('calc_dcf','var'),    calc_dcf = []; end
if isempty(calc_dcf),           calc_dcf = true; end
if ~exist('plt','var'),         plt = []; end
if isempty(plt),                plt = true; end
if ~exist('gdt2','var'),        gdt2 = []; end
if isempty(gdt2),               gdt2 = 4; end   % downsampling gradients
if gdt2<4, error('gdt2(=%g)<4',gdt2); end
if ~exist('order','var'),       order = []; end
if isempty(order),              order = 0; end
if length(order)<2 & order(1)==2,             order = [order 13]; end
if length(order)<2 & order(1)==3,             order = [order 1]; end
if ~exist('do_rot_file','var'), do_rot_file = []; end
if isempty(do_rot_file),        do_rot_file = false; end
if ~exist('bonus_spec','var'),       bonus_spec = []; end
if isempty(bonus_spec),              bonus_spec = 0; end
if bonus_spec<0 && do_rot_file
    warning('bonus spec requires saving all gradients, switching to that method...')
    do_rot_file = false;
end


%% convert to standard SI units
if gmax<1, warning('gmax(=%g)<1',gmax); input(''); end
if fov<1,  warning('fov(=%g)<1',fov);   input(''); end
fov = fov*1d-3;           % [mm] -> [m]
gdt = gdt*1d-6;           % [us] -> [s]
gdt2 = gdt2*1d-6;
gmax = gmax*1d-3;         % [mT/m] -> [T/m]


%% design parameters
rgamma = abs(gyrogamma('1h')/gyrogamma(nucleus));
fov_mns = fov/rgamma;
fprintf('gamma ratio=%g\n',rgamma);
fprintf('scaled FOV (fov_mns) = %g [m]\n',fov_mns);


%% design 2D spiral
fprintf('Calculating 2D spiral\n');
arms = 1; 
goup = false;
godown = false;
for l=1:1000
    [gspir,nn] = bnispiralgen(fov_mns,mtx,0,arms,gmax*1d3*fufa_gmax,...
        smax*fufa_smax,vardens,balanced+1);
    if nn(1)*4d-3>tacq
        if godown, break; end
        arms = arms+1;
        goup = true;
    else
        if goup, break; end
        arms = arms-1;
        godown = true;
    end
    if arms==0
        arms = 1;
        break;
    end
    if l==1000
        warning('tact=%g[ms] did not reached tacq=%g[ms] after 1000 iterations\n',...
            nn(1)*4d-3,tacq);
    end
end

gspir = gspir*1d-3;                                   % [mT/m] -> [T/m]
kspir = gdt*cumsum(gspir,1);                          % [s*T/m]
kspir = kspir(:,1) + 1i*kspir(:,2);                   % convert to complex


%% calculate required sampling BW
kdt1 = 2d6*pi/(max(sqrt(sum(gspir.^2,2)))*gyrogamma('1h')*fov_mns*samp(1));
fprintf('Full sampling freq-dir: kdt = %g [us]\n',kdt1);
if kdt1>4
    kdt = floor(kdt1/4)*4;
else
    kdt = 2;
end
if kdt<4
    kdt = 2;
    if kdt1<kdt, warning('Undersampling will occur; reduce gmax'); end
end
fprintf('Setting kdt = %g [us]\n',kdt);

kdt = kdt*1d-6;           % [us] -> [s]
dki = kdt/gdt;
if dki<0.5, warning('kdt: dki(=%g)<0.5',dki); end
if (abs(dki-round(dki))>1d-10)
    if abs(dki-0.5)>1d-10
        warning('kdt(=%g[us]) not multiple of gdt(=%g[us])',kdt,gdt);
    end
else
    dki = round(dki);
end
ng_off = nk_off*dki;

%% account for over/undersampling arms
init_arms=arms;
arms = ceil(init_arms*samp(2));


%% rotate trajectory
phi = [];
if order(1)==3 %ga https://onlinelibrary.wiley.com/doi/10.1002/jmri.28187    
    if order(2) > arms
        warning('Tiny golden angle %g (%gdeg) requires a minimum of %g arms for uniform coverage',order(2),mod((360/(GR+order(2)-1)),360),order(2));
        order(2) = arms;
        warning('Set tiny golden angle to %g arms (%gdeg)\n',arms,mod((360/(GR+order(2)-1)),360));
    end
    tyga_deg = mod((360/(GR+order(2)-1)),360);    % note order(2)=1 gives 222.5deg rotation while order(2)=2 give the more commonly used angle of 137.51
end
if arms>1
    if do_rot_file
        phi = mod(-(0:arms-1)/arms*360,360);
        if order(1)==3             
            phi = mod(tyga_deg*(0:arms-1),360);
        end
    else
        ksp = bsxfun(@times,kspir,exp(1i*(0:arms-1)/arms*2*pi));
        if order(1)==3
            phi = mod(tyga_deg*(0:arms-1),360);
            [x,y]=pol2cart(deg2rad(phi),1);
            phi = complex(x,y);
            ksp = bsxfun(@times,kspir,phi);
        end
    end
end


%% re-order trajectory
orderstr = '';
switch order(1)
    case 0
    case 1         % random trajectory order
        ind_order = randperm(arms);
        ksp = ksp(:,ind_order,:);
        orderstr = '_R';
    case 2         % mix every order(2)-th spiral
        if do_rot_file, warning('do_rot_file && order(1)==2: not implemented'); end
        ind_order = [];
        for ll=1:order(2)
            ind_order = [ind_order,ll:order(2):arms];
        end
        if length(ind_order)~=arms
            error('length(ind_order)(=%d)~=arms(=%d)',...
                length(ind_order),arms); 
        end
        if any(diff(sort(ind_order))==0)
            error('any(diff(sort(ind_order))==0)');
        end
        ksp = ksp(:,ind_order,:);
        orderstr = sprintf('_M%d',order(2));
    case 3
        %order taken into account above
        orderstr = sprintf('_GA%d',order(2));
    otherwise
        warning('order(=%d) unknown',order(1));
end

%% bonus spec
if bonus_spec>0
    arms = arms+bonus_spec;
    ksp_zeros = zeros(nn(2),bonus_spec);
    ksp = [ksp, ksp_zeros];
    ind_order(1,end-bonus_spec:end)=(arms-bonus_spec):arms;
end

%% calculate back to gradient trajectory
grad = [zeros(1+ng_off,arms); diff(ksp,1,1)/gdt];
if mod(size(grad,1),2)       % adding 0 if odd
    grad = [grad ; zeros(1,arms)];
end
grad = cat(3, real(grad), imag(grad), zeros(size(grad))); % convert from r/i to x/y/z


%% k-space trajectory
kfac = fov_mns/mtx*gyrogamma('1H')/2/pi;
if dki>1
    kk = kfac*ksp(1:dki:nn(1),:,:); 
else
    kk = kfac*ksp(1:nn(1),:,:);
end
if abs(dki-0.5)<1d-10
    kk = interp1((1:nn(1)),kk,(1:0.5:nn(1)+0.5),'spline');
end
nk = size(kk,1);
k = cat(3, reshape(real(kk),[1, nk*arms,1]), reshape(imag(kk),[1, nk*arms,1]), zeros(1, nk*arms,1));
kk = reshape(k,nk,arms,3);


%% calculate density compensation
if calc_dcf
    fprintf('Calculating density compensation function\n');
    tmp_k = k(:,1:((arms-bonus_spec)*nk),:); % don't include spec
    dcf = calc_sdc3(tmp_k,mtx,8);
    dcf = [dcf, zeros(1,bonus_spec*nk)]; % add spec as 0s
    dcf = reshape(dcf,[nk arms]);
    kkabs = sqrt(sum(kk(:,1,:).^2,3));      % intlv index does not matter
    ii = (kkabs>0.4);%&(kkabs<0.45);
    dcf_thresh = mean(dcf(ii,:),1);
    % dcf_thresh = sort(dcf(kkabs<0.45,:),1,'descend');
    % ii = (1:round(0.2*nk));
    % dcf_thresh = mean(dcf_thresh(ii,:),1);
    %     for li=1:nintlv
    %         dcf(dcf(:,li)>dcf_thresh(1,li),li) = dcf_thresh(1,li);
    %     end
    for li=1:arms
        ii = (dcf(:,li)>dcf_thresh(1,li)) & (kkabs>0.45);
        dcf(ii,li) = dcf_thresh(1,li);
    end    
    dcf = reshape(dcf,[1 nk*arms]);
else
    dcf = [];
end


%% calculate time and index list
acq_pts = ceil((length(gspir)+ng_off)/dki/32)*32; % round to multiple of 32
if acq_pts>16384
    warning('#sampling pts/exc (=%g)>16384: exceeds fidall limit',acq_pts);
end
inds = false(1,acq_pts);               % index for acq data (single arm)
inds(1,(nk_off+(1:nk))) = true;
ind = repmat(inds,[arms 1]);         % index acq all interleaves
ts = (0:nk-1)*kdt;                     % time of single trajectory acq [s]
t = repmat(ts,[arms 1]);             % time of all trajectories [s]

%% convert back to 2D
k=k(:,:,1:2);
grad=grad(:,:,1:2);

    
%% print info about waveform
gmax_act = max(max(sqrt(sum(grad.^2,3))));
smax_act = max(max(sqrt(sum(diff(grad,[],1).^2,3))))/gdt;
fprintf('actual gmax = %g [mT/m]\n',gmax_act*1d3);
if (gmax_act>gmax), warning('gmax exceeded'); end
fprintf('actual smax = %g [T/m/s]\n',smax_act);
if (smax_act>smax), warning('smax exceeded'); end

fprintf('Sequence details\n');
fprintf('Acq BW = %g [kHz] (full)\n',1d-3/kdt);
fprintf('gdt = %g [us]; kdt = %g [us]\n',gdt*1d6,kdt*1d6);
fprintf('g_pts = %dx%dx%d; k_pts = %d; acq_pts/exc = %g\n',...
    size(grad,1),size(grad,2),size(grad,3),size(k,2),acq_pts);
fprintf('t_seq = %g [ms]; arms = %d\n',nn(2)*gdt*1d3,arms);


%% checks
if size(k,1)~=1,   error('size(k,1)(=%g)~=1',  size(k,1)); end
if calc_dcf
    if size(dcf,1)~=1, error('size(dcf,1)(=%g)~=1',size(dcf,1)); end
    if size(dcf,2)~=size(k,2)
        warning('size(dcf,2)~=size(k,2): interpolating dcf');
        dcf = interp1(linspace(0,1,size(dcf,2)),dcf,linspace(0,1,size(k,2)));
    end
end
if size(k)~=sum(ind(:)), error('size(k)~=sum(ind(:))'); end


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
        fname = sprintf('spiralbni_%s_fov%g_mtx%g_nexc%g_bonus%g_vd%gstart%gstop%gr%g_kdt%g_gmax%g_smax%g_dur%g%s%s', ...
            nuc,round(fov*1d3),mtx,arms,bonus_spec,vardens(1),vardens(2),vardens(3),vardens(4),kdt*1d6,...
            round(gmax_act*1d3),round(smax_act),...
            round(nn(2)*gdt*1d4)/10,balstr,orderstr);
        fname = strrep(fname,'.','p');
        if abs(gdt2*1d6-4)>1d-10, fname = [fname '_gdt' num2str(gdt2*1d6)]; end
    else
        fname = [];
    end
end


%% downsampling grad
if gdt2>gdt
    fprintf('Downsampling gradient waveforms from %g to %g [us]\n',...
        gdt*1d6,gdt2*1d6);
    if abs(gdt2/gdt-floor(gdt2/gdt))>1d-10
        error('gdt2(=%g) not multiple of gdt(=%g)',gdt2,gdt);
    end
    dgi = round(gdt2/gdt);
    grad = grad(1:dgi:end,:,:);
    if mod(size(grad,1),2)       % adding 0 if odd
        ng_add = 3;
    else
        ng_add = 2;
    end
    grad = [grad ; zeros(ng_add,nintlv,3)];

end


%% export waveforms
if ~isempty(fname)
    desc = sprintf('design_spiral_bni(%g,%g,%g,[%g,%g],[%g,%g,%g,%g],fname,%g,%g,%g,%s,%g,%g,%g,%g)\n',...
        fov*1d3,mtx,tacq,samp(1),samp(2),vardens(1),vardens(2),vardens(3),vardens(4),...
        arms,gmax*1d3,smax,nucleus,balanced,calc_dcf,plt,gdt2*1d6);
    fprintf('Saving gradient + mat file to\n\t%s\n',fname);
    fprintf('%s\n',desc);
    if ~calc_dcf, warning('dcf missing'); end
    
    if do_rot_file
        out = write_ak_wav([fname '.wav'],grad(:,1,:),1/kdt,fov_mns,desc,acq_pts,...
            [],[],gdt2*1d6);
        write_fdl([fname '_phi.fdl'],  phi,  'phi');
    else
        out = write_ak_wav([fname '.wav'],grad,1/kdt,fov_mns,desc,acq_pts,...
            [],[],gdt2*1d6);
    end

    if isfield(out,'wave'), out = rmfield(out,'wave'); end
    if isfield(out,'grad'), out = rmfield(out,'grad'); end
    
    nviews = arms;
    dcf_range = max(dcf);
    npix = mtx;
    
    dcf = single(dcf);
    k = single(k);

    
    save(fname,'out','k','dcf','dcf_range','ts','inds','ind',...
        'fov','mtx','gdt','kdt','gmax','smax','nucleus',...
        'ng_off','rgamma','fov_mns','arms','nviews','npix');
else
    out.gmax = gmax;
    out.smax = smax;
    out.gdt = gdt;
    out.kdt = kdt;
    out.grad = grad;
    out.bw = 1/kdt;
    out.fov = fov;
end


%% plotting
if plt
    fprintf('Plotting\n');
    if arms>500
        ii2 = round(linspace(1,arms,500));
    else
        ii2 = (1:arms);
    end
    figure(100); clf
    yl = {'Gx [mT/m]','Gy [mT/m]','Gz [mT/m]'};
    for l=1:2
        subplot(2,1,l);
        plot((1:size(grad,1))*gdt2*1d3,grad(:,ii2,l)*1d3);
        xlabel('time [ms]');
        ylabel(yl{l})
        grid on
    end
    
    figure(101); clf
    plot(kk(:,ii2,1),kk(:,ii2,2),'.','MarkerSize',1);
    grid on;
    title('k-space');
    
    if calc_dcf
        figure(102); clf
        dcftmp = reshape(dcf,[nk arms]);
        plot(ts*1d3,dcftmp(:,ii2));
        grid on;
        title('dcf');
        xlabel('time [ms]');
        ylabel('dcf');
    end
end

if nargout==1
    k = fname;
end

end      % design_spiral_bni.m
