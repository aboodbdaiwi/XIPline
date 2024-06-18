function [k,dcf,t,ind,out,grad]=design_spiral(fov,npix,arms,ksamp,...
    fname,gmax,smax,nucleus,acq_round2n,do_rot_file,balanced)
%DESIGN_SPIRAL Design a spiral with delayed acq for fast CSI
%[k,dcf,t,ind,out,grad]=design_spiral(fov,npix,arms,ksamp,...
%    fname,gmax,smax,nucleus,acq_round2n,do_rot_file,balanced)
%
%        fov  Field of view                                  [mm]
%       npix  #pixels (Cartesian resolution after gridding)
%       arms  #spatial interleaves                    (1)
%      ksamp  k-space sampling time                   (16)   [us]
%      fname  Filename: if given->write_ak_wav 
%             if logical true: generate name          (false)
%       gmax  Max gradient amp                        (33)   [mT/m]
%       smax  Max slew rate                           (120)  [T/m/s]
%    nucleus  Nucleus                                 ('1H')
%acq_round2n  Round up #acq pts to to 2^n             (true)
%do_rot_file  Write single static spiral into gradient 
%             waveform file & rotate via vap_phiXX.fdl (false)
%   balanced  Balancing gradient area                 (true)
%
%          k  k-space trajectory  [-0.5..0.5]
%        dcf  Density compensation function (calculated with vornoi_area)
%          t  Time  (nexc x recon-pts)
%        ind  Index (2nd dim) for k-space points on spiral (excl ramps, etc)
%        out  Output structure of wrt_wavs
%       grad  Gradient waveform [T/m] with dt=4us 
%
% 10/2020 Rolf Schulte
if (nargin<1), help(mfilename); return; end


%% fixed parameters
fufa_gmax = 0.99;
fufa_smax = 0.99;
g_offset = 16;               % #grad pts offset for start of trajectory
gsamp = 4e-6;                % gradient update time [s]


%% input parameters
if ~exist('arms','var'),        arms = []; end    % #spatial interleaves
if isempty(arms),               arms = 1; end
if ~exist('ksamp','var'),       ksamp = []; end
if isempty(ksamp),              ksamp = 16; end   % [us] k-space dwell time
if ksamp<2, warning('ksamp(=%g)<2',ksamp); input(''); end
if ~exist('fname','var'),       fname = []; end   % for waveform file
if ~exist('gmax','var'),        gmax = []; end
if isempty(gmax),               gmax = 33; end   % [mT/m]
if ~exist('smax','var'),        smax = []; end
if isempty(smax),               smax = 120; end  % [T/m/s]
if ~exist('nucleus','var'),     nucleus = []; end
if isempty(nucleus),            nucleus = '1H'; end
if ~exist('acq_round2n','var'), acq_round2n = []; end
if isempty(acq_round2n),        acq_round2n = true; end
if ~exist('do_rot_file','var'), do_rot_file = []; end
if isempty(do_rot_file),        do_rot_file = false; end
if ~exist('balanced','var'),    balanced = []; end
if isempty(balanced),           balanced = true; end


%% checks
if gmax<1, warning('gmax(=%g)<1',gmax); input(''); end
if fov<1,  warning('fov(=%g)<1',fov);   input(''); end
dki = ksamp/gsamp*1d-6;
if dki<0.5
    error('ksamp: dki<0.5');
end
if (abs(dki-round(dki))>1d-10)
    if abs(dki-0.5)>1d-10
        error('ksamp (=%g us) not multiple of gsamp (=%g us)',...
            ksamp*1d6,gsamp*1d6);
    end
else
    dki = round(dki);
end
k_offset = g_offset/dki;


%% convert to standard SI units
fov = fov*1d-3;           % [mm] -> [m]
ksamp = ksamp*1d-6;       % [us] -> [s]
gmax = gmax*1d-3;         % [mT/m] -> [T/m]
res = fov/npix;           % [m] resolution


%% design spiral
rgamma = abs(gyrogamma('1h')/gyrogamma(nucleus));
fov_mns = fov/rgamma;
res_mns = res/rgamma;
fprintf('gamma ratio=%g\n',rgamma);
fprintf('scaled FOV (fov_mns) = %g [m]\n',fov_mns); 
fprintf('scaled res (res_mns) = %g [mm]\n',res_mns*1d3);
gmax_nyquist = 2*pi/(gyrogamma('1h')*ksamp*fov_mns);
fprintf('gmax_nyquist = %g [mT/m]\n',gmax_nyquist*1d3);
if (gmax_nyquist<gmax)
    fprintf('Attention: approaching sampling BW limited regime\n');
    if (ksamp>gsamp)
        fprintf('!!! Undersampling will occur: reduce ksamp !!!\n'); 
        input('press key to continue');
    end
end
pause(1);
[~,g1,~,t1] = vds(fufa_smax*smax*1e2,fufa_gmax*gmax*1e2,gsamp, ...
   arms,[fov_mns*1e2,0],1/(2*res_mns*1e2));


%% calculate single interleave
gspir = 1e-2*g1;                % convert to SI unit [T/m]
kspir = gsamp*cumsum(gspir);    % [s*T/m]


%% rewinders
if balanced
    [grewx,ngx] = gradient_lobe(-real(kspir(1,end)),real(gspir(1,end)),...
        fufa_gmax*gmax/sqrt(2),fufa_smax*smax/sqrt(2),gsamp);
    [grewy,ngy] = gradient_lobe(-imag(kspir(1,end)),imag(gspir(1,end)),...
        fufa_gmax*gmax/sqrt(2),fufa_smax*smax/sqrt(2),gsamp);
    
    grew = zeros(1,max([ngx,ngy])+1);
    grew(1:ngx) = 1*grewx;
    grew(1:ngy) = grew(1:ngy)+1i*grewy;
    
else
    % nrew = ceil(gmax/smax/gsamp/fufa_gmax);
    nrew = ceil(abs(gspir(1,end))/smax/gdt/fufa_gmax)+1;
    grew = linspace(real(gspir(1,end)),0,nrew) + ...
        1i*linspace(imag(gspir(1,end)),0,nrew);
end


%% construct full (not yet rotated) gradient waveform
gg = [zeros(1,g_offset),gspir,grew].';
ng = size(gg,1);
nk = ceil(length(kspir)/dki);


%% calculate time and index list
if acq_round2n
    acq_pts = 2^ceil(log2(ng/dki));
else
    acq_pts = ceil(length(gspir)/dki/2)*2;
end
if acq_pts>16384
    warning('#sampling pts/exc (=%g)>16384: exceeds fidall limit',nk);
end
ind = false(arms,acq_pts);             % index for acquired data 
ind(:,(k_offset+(1:nk))) = true;
t = repmat((0:nk-1)*ksamp,[arms 1]);   % time of acq pts [s]


%% rotate trajectory
phi = [];
if arms>1
    if do_rot_file
        phi = mod(-(0:arms-1)/arms*360,360);
    else
        gg = bsxfun(@times,gg,exp(1i*(0:arms-1)/arms*2*pi));
    end
end

%% k-space trajectory
if dki>1, k = kspir(1:dki:length(kspir))/2/pi*res_mns*gyrogamma('1H'); 
else,     k = kspir/2/pi*res_mns*gyrogamma('1H');
end
if abs(dki-0.5)<1d-10
    k = interp1((1:length(k)),k,(1:0.5:length(k)+0.5),'spline');
end
if arms>1
    % tmp = [];
    % for l=1:arms, tmp = [tmp , k*exp(1i*(l-1)/arms*2*pi)]; end
    % k = tmp;
    k = bsxfun(@times,k,exp(1i*(0:arms-1)/arms*2*pi).').';
    k = k(:).';
end


%% calculate density compensation function (dcf)
% approximate dcf via analytical radial expression to determine threshold
n2 = arms; mtx = npix; 
tmp = pi/n2*(k(1,2:end).*conj(k(1,2:end))-k(1,1:end-1).*conj(k(1,1:end-1)))*mtx^2;
dcf_thresh = max(tmp);

dcf = voronoi_area(k*npix,'','',dcf_thresh);


%% print info about waveform
gmax_act = max(max(abs(gg)));
smax_act = max(max(abs(diff(gg,[],1))))/gsamp;
fprintf('actual gmax  = %g [mT/m]\n',gmax_act*1d3);
if (gmax_act>gmax), warning('gmax exceeded'); end
fprintf('actual smax  = %g [T/m/s]\n',smax_act);
if (smax_act>smax), warning('smax exceeded'); end

fprintf('Sequence details\n');
fprintf('Acq BW = %g [kHz] (full)\n',1d-3/ksamp);
fprintf('gsamp = %g [us]; ksamp = %g [us]\n',gsamp*1d6,ksamp*1d6);
fprintf('g_pts = %gx%g; k_pts = %g; acq_pts/exc = %g\n',...
    size(gg,1),size(gg,2),size(k,2),acq_pts);
t_arm = t1(1,end); t_rew = (length(grew)+1)*gsamp;
fprintf('t_arm = %g [ms]; t_rew = %g [ms]\n',t_arm*1d3,t_rew*1d3);
fprintf('t_seq = %g [ms]\n',size(gg,1)*gsamp*1d3);


%% checks
if any(size(dcf)~=size(k))
    warning('size(dcf)~=size(k): interpolating dcf'); 
    dcf = interp1(linspace(0,1,size(dcf,2)),dcf,linspace(0,1,size(k,2)));
end
if size(k)~=sum(ind(:)),     error('size(k)~=sum(ind(:))'); end
if prod(size(t))~=size(k,2), error('prod(size(t))~=size(k,2)'); end


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
        fname = sprintf('spiral_%s_fov%g_mtx%g_arms%g_kdt%g_gmax%g_smax%g_dur%g%s', ...
            nuc,round(fov*1d3),npix,arms,ksamp*1d6,round(gmax_act*1d3),...
            round(smax_act),round(size(gg,1)*gsamp*1d4)/10,balstr);
        fname = strrep(fname,'.','p');
        if do_rot_file, fname = [fname '_rofi']; end
    else
        fname = [];
    end
end


%% export waveforms
if ~isempty(fname)
    desc = sprintf('design_spiral(%g,%g,%g,%g,fname,%g,%g,%s,%g,%g,%g)\n',...
        fov*1d3,npix,arms,ksamp*1d6, gmax*1d3,smax,nuc,acq_round2n,...
        do_rot_file,balanced);
    fprintf('Saving gradient + mat file to\n\t%s\n',fname);
    fprintf('%s\n',desc);
    out = write_ak_wav([fname '.wav'],gg,1/ksamp,fov_mns,desc);
    save(fname,'out','k','dcf','t','ind','fov',...
        'npix','gsamp','ksamp','gmax','smax','nucleus',...
        't_rew','t_arm','g_offset','gmax_nyquist','rgamma','fov_mns','phi');
    if do_rot_file, write_fdl([fname '_phi.fdl'],phi,'phi'); end
else
    out.gmax = gmax;
    out.smax = smax;
    out.gdt = gsamp;
    out.kdt = ksamp;
    out.grad = gg;
    out.bw = 1/ksamp;
    out.fov = fov;
end


%% reshape + output grad
grad = zeros(ng,size(gg,2),2);
grad(1:ng,:,1) = real(gg);
grad(1:ng,:,2) = imag(gg);


%% plotting
if true
    clf
    yl = {'Gx [mT/m]','Gy [mT/m]'};
    for l=1:2
        subplot(2,1,l);
        plot((1:size(grad,1))*gsamp*1d3,grad(:,:,1)*1d3);
        xlabel('time [ms]');
        ylabel(yl{l})
        grid on
    end
end


end  % main function design_spiral.m
