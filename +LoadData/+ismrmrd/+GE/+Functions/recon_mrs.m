function [spec,hz,ddd,time,outp] = recon_mrs(d,h,zf,lb,cohave,do_ecc,...
    ws_bounds,plt,verb,lcm_basis)
%RECON_MRS  Reconstruction of (single voxel) MRS data 
%  Coil-combination, eddy current compensation, coherent averaging, HSVD
%  water subtraction, writing out to LCModel
%
%[spec,hz,ddd,time,outp]=recon_mrs(d,h,zf,lb,cohave,do_ecc,ws_bounds,plt,verb,lcm_basis)
%
%         d  Raw p-file data [nf,ns,1,1,1,nc]    (or pfile fname)
%            nf=#frames; ns=#spectral pts; nc=#coils
%         h  p-file header                       (or empty)
%        zf  zero-filling                                (default=[])
%        lb  line-broading (exponential Gaussian) [Hz]   (default=[])
%    cohave  coherent averaging of time series           (default=4)
%    do_ecc  phase deconvolution eddy current correction (default=true)
%            requires reference (not water suppressed) frames
% ws_bounds  bounds for HSVD water suppression           (default=[-50 50])
%       plt  plotting                                    (default=true)
%      verb  verbose                                     (default=true)
% lcm_basis  filename of LCModel basis -> save files     (default=[])
%
%      spec  Spectrum        [ns,nf2]  (spec(:,1)=MRS; spec(:,2)=reference)
%        hz  Frequency axis  [Hz]
%       ddd  FID: coil-combined and averaged
%      time  Time axis       [s]
%      outp  misc output structure
%
% 9/2015 Rolf Schulte
if nargin<1, help(mfilename); return; end


%% input/default parameters
if ~exist('zf','var'),       zf = []; end
if ~exist('lb','var'),       lb = []; end
if ~isempty(lb), if length(lb)==1, lb(2) = 0; end; end
if ~exist('cohave','var'),   cohave = []; end
if isempty(cohave),          cohave = 4;  end
if ~exist('do_ecc','var'),   do_ecc = []; end
if ~exist('ws_bounds','var'),ws_bounds = []; end
if isempty(ws_bounds),       ws_bounds = true; end
if islogical(ws_bounds)
    if ws_bounds, ws_bounds = 50*[-1 1]; 
    else, ws_bounds = []; end
end
if ~exist('plt','var'),      plt = []; end
if isempty(plt),             plt = true; end
if ~exist('verb','var'),     verb = []; end
if isempty(verb),            verb = true; end
if ~exist('lcm_basis','var'), lcm_basis = []; end


%% reading in data, if pfile name is given
if ~isnumeric(d)
    if exist(d,'file'), [d,h] = read_p(d);
    else, warning('strange input d/file not existing');
    end
end


%% miscelaneous parameters
csi_dims = h.rdb_hdr.csi_dims;         % #CSI dims (must be 0)
nf   = h.rdb_hdr.nframes;              % acquired frames (size(d,1)
ns   = h.rdb_hdr.frame_size;           % acq pts (user1)
bw   = h.rdb_hdr.spectral_width;       % BW [Hz] (user0)
nref = h.rdb_hdr.user19;               % #reference frames (w/o WS)
nspec = h.image.user4/h.image.nex;     % #spectral frames  (w/ WS)
nc = size(d,6);                        % #coils
time = (0:ns-1)/bw;                    % time axis
if isempty(zf), zf = ns; end           % zero filling; default=off

if isempty(do_ecc)
    if (nref==0), do_ecc = false;
    else, do_ecc = true; end
end
if do_ecc && (nref==0)
    warning('No reference data; switching off eddy current compensation');
    do_ecc = false;
end
signal_bound = -100;                   % cut-off freq for SNR calc [Hz]


%% print info
if verb
    fprintf('TE=%g [ms]; TR=%g [ms]; nacqs(cv4)=%g; BW(cv0)=%g [Hz]; ns(cv1)=%g\n',...
        h.image.te*1d-3,h.image.tr*1d-3,h.image.user4,h.image.user0,ns);
    fprintf('nf=%g, nc=%g, nref=%g\n',nf,nc,nref);
    fprintf('size(d) = %s\n',num2str(size(d)));
end


%% data checks
nrcvr = (h.rdb_hdr.dab(2)-h.rdb_hdr.dab(1))+1;   % #receiver channels
if nf~=size(d,1),  warning('nf(=%g) ~= size(d,1)(=%g)',nf,size(d,1)); end
if ns~=size(d,2),  warning('ns(=%g) ~= size(d,2)(=%g)',ns,size(d,2)); end
if nc~=nrcvr,      warning('nc(=%g) ~= nrcvr(=%g)',nc,nrcvr); end
if csi_dims~=0,    warning('csi_dim(=%g)~=0: use recon_csi.m for CSI data',csi_dim); end
if (size(d,3)~=1), warning('size(d,3)(=%g)~=1',size(d,3)); end 
if (size(d,4)~=1), warning('size(d,4)(=%g)~=1',size(d,4)); end 
if (size(d,5)~=1), warning('size(d,5)(=%g)~=1',size(d,5)); end 
if ((nref+nspec)~=nf), warning('(nref(%g)+nspec(%g)~=nf(=%g)',nref,nspec,nf); end

% unchop data (in case opnex=1)
if ischop(h) && (nf>1), d(2:2:nf,:) = -d(2:2:nf,:); end


%% coil combination
if nref>0, ind1 = (1:nref);
else,      ind1 = (1:nf); end
dd = mrs_coil_combine(d,ind1,[],[],verb);      % signal-weighted
% dd = mrs_coil_combine(d,ind1,[],100,verb);     % SNR-weighted
% dd = csi_coil_combine(d,ind1);                   % signal-weighted svds
warning('check coil combination');

%% eddy current correction - phase deconvolution
if do_ecc
    dd_ref = mean(dd(1:nref,:),1);
    if verb, fprintf('Eddy current correction by phase deconvolution\n'); end
    dd(((nref+1):nf),:) = ecc(dd(((nref+1):nf),:),dd_ref,[],plt);
end


%% averaging time series
if cohave          % phase + frequency alignment
    dd = mrs_coave(dd,bw,cohave,nref,verb);
end
if nref>0          % actual averaging; 1=reference MRS; 2=ws MRS
    ddd(1,:) = mean(dd(((nref+1):nf),:),1);
    ddd(2,:) = mean(dd((1:nref),:),1);
else
    ddd = mean(dd,1);
end


%% water subtraction via HSVD; only for water suppressed frames
if ~isempty(ws_bounds)
    if ischar(ws_bounds)
        if verb, fprintf('Water subtraction\n'); end 
        ddd(1,:) = h2osupGE(ddd(1,:),ddd(2,:));
    else
        if verb
            fprintf('HSVD water suppression: bounds=(%g %g) [Hz]\n',...
                ws_bounds(1),ws_bounds(2));
        end
        [~,~,~,ddd(1,:)] = ...
            hsvd(ddd(1,:).',25,ws_bounds,bw,false);
    end
end


%% line broadening
if ~isempty(lb)
    fprintf('Apodising: exponential lb=%g [Hz]; Gaussian lb=%g [Hz]\n',lb(1),lb(2));
    apo = lb_fun(ns,bw,lb);
    ddd = ddd.*repmat(apo,[size(ddd,1) 1]);
end


%% spectral reconstruction
hz  = -(-zf/2:zf/2-1)*bw/zf;             % frequency axis
spec = ifftshift(fft(ddd,zf,2),2);


%% SNR
[snr,signal,noise,f_max,ind_noise] = mrs_snr(ddd(1,:),spec(1,:),bw,signal_bound,[],false);


%% line-width
[lw,out] = fwhm(spec,bw,[],false,false);


%% print info
if verb
    vol = h.rdb_hdr.roilenx*h.rdb_hdr.roileny*h.rdb_hdr.roilenz*1d-3;
    acq_dur = h.image.tr*1d-6*h.image.user4;
    fprintf('Voxel size = (%g, %g, %g) [mm^3] = %g [ml]\n',...
        h.rdb_hdr.roilenx,h.rdb_hdr.roileny,h.rdb_hdr.roilenz,vol);
    fprintf('Acqusition time = %g [s]\n',acq_dur);
    fprintf('SNR = %g/%g = %.3g\n',signal,noise,snr);
    fprintf('SNR/acq_dur/vol = %.3g [ml*min]^-1\n',snr/vol/acq_dur*60);
    fprintf('linewidth (FWHM) = %.3g [Hz]',lw(1));
    if nref>0, fprintf('; ref = %.3g [Hz]',lw(2)); end
    fprintf('\n');
    fprintf('**************************************************\n');
    outp.label = ...
        sprintf('P-file\tdate\tTE\tTR\tcv4\tBW\tns\tnc\tdur\tvol\tSNR\tunitSNR\tLW\n');
    outp.str = ...
        sprintf('P%.5d.7\t%s\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%.3g\t%.3g\t%.3g\n',...
        h.image.rawrunnum,h.rdb_hdr.scan_date,...
        h.image.te*1d-3,h.image.tr*1d-3,h.image.user4,h.image.user0,ns,nc,...
        acq_dur,vol,snr,snr/vol/acq_dur*60,lw(1));
    fprintf('%s',outp.label);
    fprintf('%s',outp.str);    
end


%% plotting
if plt
    clf
    if nref>0
        sca = max(abs(spec),[],2);
        ss = abs(spec)./repmat(sca,[1 zf]);
        plot(hz,ss(1,:),'b',hz,ss(2,:),'r',...
            hz(ind_noise),ss(1,ind_noise),'cx',...
            f_max(1),1,'cs',...
            out.freqs(1,:),abs(ss(1,out.ind_max(1,:))),'c*',...
            [out.fp(1,:) out.fm(1,:)],...
            repmat(abs(ss(1,out.ind_max(1,:)))/2,[2 1]),'c*-',...
            out.freqs(2,:),abs(ss(2,out.ind_max(2,:))),'m*',...
            [out.fp(2,:) out.fm(2,:)],...
            repmat(abs(ss(2,out.ind_max(2,:)))/2,[2 1]),'m*-',...
            'MarkerSize',8,'LineWidth',1,'MarkerFaceColor','r');
        legend('spectrum','reference');
    else
        ss = abs(spec)./signal;
        plot(hz,ss(1,:),'b',...
            hz(ind_noise),ss(1,ind_noise),'cx',...
            f_max(1),1,'cs',...
            out.freqs(1,:),abs(ss(1,out.ind_max(1,:))),'c*',...
            [out.fp(1,:) out.fm(1,:)],...
            repmat(abs(ss(1,out.ind_max(1,:)))/2,[2 1]),'c*-',...
            'MarkerSize',8,'LineWidth',1,'MarkerFaceColor','r');
        legend('spectrum');
    end
    set(gca,'XDir','reverse');
    axis([min(hz) max(hz) 0 1]);
end


%% write out to LCModel
if ~isempty(lcm_basis)
    if islogical(lcm_basis)
        if ~lcm_basis, return; end
        if h.image.te*1d-3~=35
            warning('TE(=%g)~=35[ms] (TE of default basis)',h.image.te*1d-3); 
        end
        lcm_basis = '/home/schulte/data/basis/lcmodel_simulated/press35/gamma_basis_TE35_B0_30_TE1_18_BW5000_N4096.basis';
    end
    fprintf('Writing out raw data and control file for LCModel\n'); 
    if ~isempty(lb)
        warning('You should use raw data for LCModel; remove line-broadening');
        fprintf('   lb exp=%g; Gaussian=%g [Hz]\n',lb(1),lb(2));
    end
    if (zf~=ns)
        warning('You should not apply zero-filling (ns=%g; zf=%g)',ns,zf);
    end
    if ((nref>0) && ~do_ecc)
        write_lcmodel(ddd(1,:),h,lcm_basis,ddd(2,:));
    else
        write_lcmodel(ddd(1,:),h,lcm_basis);
    end
end

end      % recon_mrs.m
