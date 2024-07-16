function [bb,bbabs,k,dd] = recon_zte(d,h,mtx_reco,delay,lbt,fname,f0,...
    coil_arr,plt,do_ori,overrecofa,do_single,max_waspi,do_phi0corr,...
    do_tecorr,save_coils,do_homodyne,export_dcm)
%RECON_ZTE Reconstruct ZTE 3D radial data from 3dradial psd
%[bb,bbabs,h,k] = recon_zte(d,h,mtx_reco,delay,lbt,fname,f0,coil_arr,...
%        plt,do_ori,overrecofa,do_single,max_waspi,do_phi0corr,...
%        do_tecorr,save_coils,do_homodyne,export_dcm);
%                                                               (default)
%          d  Raw (p-file) data  (or pfile fname)
%          h  Header from p-file (or empty)
%   mtx_reco  Reconstruction matrix resolution (Fourier interp) ([])
%      delay  Gradient-acquisition delay [us]                   ([])
%             two values possible for burst_mode>0 
%             []->automatic calib via Block;ISMRM11;#2816
%        lbt  Gaussian linebroadening time [Hz]                 (0)
%      fname  Print <fname>.png and save reco as <fname>.mat    ([]->no)
%             also export dicom if template.dcm is found
%         f0  Change in centre frequency                   [Hz] ([])
%   coil_arr  Reconstruct only specific coil elements           ([]->all)
%        plt  Plotting: 0=off;1+2+5+6=imagesc_row
%                       3+4+5+6=imagesc_ind3d
%                       2+4+6=invisible (for saving)            (3)
%     do_ori  Orient reconstructed data                         (true)
% overrecofa  Reconstruction factor beyond FOV                  (1)
%  do_single  Reconstruct/store data in single precision        (true)
%  max_waspi  Maximise WASPI/crop highres (for slow TR switches)
%             (default=false for specnuc=1; true for specnuc>1)
%do_phi0corr  Apply phase correction (centre of kspace) to burst echoes
%             0=off, 2=filtered phase removal, 1=mid k-space    (1)
%  do_tecorr  Correct signal level for T2* decay of inconsistent TE for
%             even acquisitions (burst_mode>0)                  (true)
% save_coils  Save individual coil images                       (false)
%do_homodyne  Homodyne reconstruction for ZTE & burst_half==1   (false)
% export_dcm  Save images in dicom format                       (1)
%
%         bb  Reconstructed data       (mtx,mtx,mtx,#timesteps,#coils);
%      bbabs  RMS coil-combined data   (mtx,mtx,mtx,#timesteps);
%          k  k-space trajectory
%
% 8/2020  Rolf Schulte
if (nargin<1), help(mfilename); return; end

timerVal = tic;               % record execution time


%% input variables
if ~exist('mtx_reco','var'), mtx_reco = []; end
if isempty(mtx_reco),        mtx_reco = 0; end
if ~exist('delay','var'), delay = []; end
if ischar(delay)
    fprintf('isstr(delay) -> applying str2num\n'); 
    delay = str2num(delay); 
end
if ~exist('lbt','var'),   lbt = []; end
if isempty(lbt),          lbt = 0; end
if ~exist('fname','var'), fname = []; end
if ~isempty(fname)
    if ~islogical(fname)
        if ~isempty(regexpi(fname,'\.7$')), fname = fname(1:end-2); end
        if ~isempty(regexpi(fname,'\.h5$')), fname = fname(1:end-3); end
        if ~isempty(regexpi(fname,'\.mat$')), fname = fname(1:end-4); end
    end
end
if ~exist('f0','var'),        f0 = []; end
if isempty(f0),               f0 = 0; end
if length(f0)~=1,             error('length(f0)~=1'); end
if ~exist('coil_arr','var'),  coil_arr = []; end
if ~exist('plt','var'),       plt = []; end
if isempty(plt),              plt = 3; end
if length(plt)~=1,            error('length(plt)~=1'); end
if ~exist('do_ori','var'),    do_ori = []; end
if isempty(do_ori),           do_ori = true; end
if length(do_ori)~=1,         error('length(do_ori)~=1'); end
if ~exist('overrecofa','var'),overrecofa = []; end
if length(overrecofa)>1,      error('length(overrecofa)>1'); end
if ~isempty(overrecofa)
    if ~isnumeric(overrecofa),error('~isnumeric(overrecofa)'); end
    if overrecofa<1,          warning('overrecofa(=%g)<1',overrecofa); end
end
if ~exist('do_single','var'), do_single = []; end
if isempty(do_single),        do_single = true; end
if length(do_single)~=1,      error('length(do_single)~=1'); end
if ~exist('max_waspi','var'), max_waspi = []; end
if length(max_waspi)>1,       error('length(max_waspi)>1'); end
if ~exist('do_phi0corr','var'),do_phi0corr = []; end
if isempty(do_phi0corr),      do_phi0corr = true; end
if length(do_phi0corr)~=1,    error('length(do_phi0corr)~=1'); end
if ~exist('do_tecorr','var'), do_tecorr = []; end
if isempty(do_tecorr),        do_tecorr = true; end
if length(do_tecorr)~=1,      error('length(do_tecorr)~=1'); end
if ~exist('save_coils','var'), save_coils = []; end
if isempty(save_coils),       save_coils = false; end
if length(save_coils)~=1,     error('length(save_coils)~=1'); end
if ~exist('do_homodyne','var'),do_homodyne = []; end
if isempty(do_homodyne),      do_homodyne = false; end
if length(do_homodyne)~=1,    error('length(do_homodyne)~=1'); end
if ~exist('export_dcm','var'),export_dcm = []; end
if isempty(export_dcm),       export_dcm = true; end
if length(export_dcm)~=1,     error('length(export_dcm)~=1'); end


%% reading in data, if pfile name is given
pfname = '';
if ~isnumeric(d)
    if exist(d,'file')
        pfname = d;
        [d,h] = read_p(d,do_single,true);
    else
        warning('strange input d/file not existing');
    end
end
if islogical(fname)
    if fname 
        fname = sprintf('P%05d',h.image.rawrunnum);
    else
        fname = [];
    end
end

%% skip d(1,:..)
nframes = h.rdb_hdr.nframes;
switch (size(d,1)-nframes)
    case 0
    case 1, d = d(2:end,:,:,:,:,:); fprintf('removing d(1,:...\n');
    otherwise
        error('nframes(=%g) ~= size(d,1)(=%g) (+1)',nframes,size(d,1));
end


%% misc variables
xres = h.rdb_hdr.da_xres;
% xtra = h.rdb_hdr.feextra;
nrcv = (h.rdb_hdr.dab(2)-h.rdb_hdr.dab(1))+1;   % #coils
nspokes_low  = h.rdb_hdr.nspokes_lowres;        % #spokes in k-space centre
nspokes_high = h.rdb_hdr.nspokes_highres;       % #spokes outside
nspokes = nspokes_low + nspokes_high;
if isempty(max_waspi)
    max_waspi = false;
    if h.image.specnuc>1, max_waspi = true; end
end

spokesPerSeg = h.image.spokesPerSeg;
% npts = xres-xtra;                    % acquired data points
fov = h.rdb_hdr.fov;                   % field in view [mm]
slthick = h.image.slthick;             % slice thickness = voxel size [mm]
mtx_acq = floor(fov/slthick+0.5);      % nominal matrix size of acquisition
if mtx_reco <= 0, mtx_reco = mtx_acq; end
if ~isempty(overrecofa), mtx_reco = floor(mtx_reco*overrecofa+0.5); end
nechoes = h.rdb_hdr.nechoes;           % enable multiple echoes
bw = h.rdb_hdr.bw*2d3*h.rdb_hdr.oversamplingfactor;   % full acq BW [Hz]
[n1,n2,n3,n4,n5,n6] = size(d);

fprintf('recon_zte input:\n');
if isempty(delay)
    delay_str = '[]';
else
    delay_str = num2str(delay);
end
fprintf('\t\tmtx_reco=%g  delay=%s  lbt=%g  fname=%s  f0=%g',...
    mtx_reco,delay_str,lbt,fname,f0);
fprintf('   coil_arr='); disp(coil_arr);
fprintf('\n\t\tplt=%g  do_ori=%g  overrecofa=%g  do_single=%g\n',...
    plt,do_ori,overrecofa,do_single);
fprintf('\t\tmax_waspi=%g  do_phi0corr=%g  do_tecorr=%g  save_coils=%g\n',...
    max_waspi,do_phi0corr,do_tecorr,save_coils);
fprintf('\t\tdo_homodyne=%g  export_dcm=%g\n',do_homodyne,export_dcm);

if do_ori
    rot90fac = h.data_acq_tab.rotate;
    trnsps = h.data_acq_tab.transpose;
    fprintf('Orienting images: rot90fac=%g; trnsps=%g\n',rot90fac(1),trnsps(1));
else
    rot90fac = 0;
    trnsps = 0;
    fprintf('No image orientation\n');
end
if ~((trnsps(1)==0)||(trnsps(1)==3))
    warning('h.data_acq_tab.transpose (=%g) not 0 or 3:\n',trnsps(1));
    fprintf('\tunknown image orientation\n');
end
if n2~=xres, warning('size(d,2)(=%g)~=xres(=%g)',n2,xres); end
if n3~=1,    warning('size(d,3)(=%g)~=1',n3); end
if n4~=nechoes
    warning('size(d,4)(=%g)~=nechoes(=%g)',n4,nechoes); 
    nechoes = n4;
end
if n6~=nrcv, warning('size(d,6)(=%g)~=nrcv(=%g)',n6,nrcv); end
if nspokes>(n1*n5), warning('nspokes(=%g)>(n1*n5)(=%g)',nspokes,n1*n5); end
if isempty(coil_arr), coil_arr = (1:nrcv); end
coil_arr = sort(floor(coil_arr));
if coil_arr(1)<1, error('coil_arr(1)<1'); end
if coil_arr(end)>nrcv, error('coil_arr(end)>nrcv'); end
if any(diff(coil_arr)==0), error('any(diff(coil_arr)==0)'); end
nc = length(coil_arr);
if (nechoes>1)
    burst_mode = h.rdb_hdr.user12;
else
    burst_mode = 0;
    if isempty(delay)
        fprintf('Delay determination requires burst; setting delay=0\n');
        delay = 0;
    end
end
burst_half = round(h.rdb_hdr.user18);
xloc = h.rdb_hdr.user26;     % off-isocentre in [mm]
yloc = -h.rdb_hdr.user27;
zloc = h.rdb_hdr.user28;
if ((abs(xloc)>0.01) || (abs(yloc)>0.01) || (abs(zloc)>0.01))
    shft = mtx_acq*[xloc yloc zloc]/fov;
    fprintf('Shifting iso centre: [%.4g,%.4g,%.4g] [mm]; ',xloc,yloc,zloc);
    fprintf('[%.3g,%.3g,%.3g] [pixel]\n',shft);
else
    shft = [];
end


%% allocate memory for reconstructed images
if do_single
    fprintf('Reconstructing/storing data as single\n');
    if ~isa(d,'single'), d = single(d); end
end


%% pre-processing of data
% unchop data first
if ischop(h)
    fprintf('Data chopped; unchopping\n');
    d(2:2:nexc,:) = -d(2:2:nexc,:); 
end


%% determine delays
if isempty(delay)
    fprintf('Determining gradient delays\n');
    % Cartesian FFT of spokes
    sig = ifftshift(ifft(fftshift(sqrt(...
        mean(d(:,:,:,2:end,:,:).*conj(d(:,:,:,2:end,:,:)),6)),2),[],2),2);
    
    if do_single
        ssig = complex(zeros(n1*n5,n2,nechoes-1,'single'));
    else
        ssig = complex(zeros(n1*n5,n2,nechoes-1));
    end
    for le=1:(nechoes-1)
        for l5=1:n5
            ssig((1:n1)+(l5-1)*n1,:,le) = sig(:,:,1,le,l5);
        end
    end
    sig = ssig(((1:nspokes_high)+nspokes_low),:,:);
    clear ssig
    
    % index to centre of k-space
    ind2 = h.rdb_hdr.feextra;
    if burst_mode>0
        ind2 = repmat([xres-h.rdb_hdr.feextra-1 ; ind2],[100 1]);
        ind2 = reshape(ind2(1:(nechoes-1),1),[1 1 (nechoes-1)]);
    end
    
    % subtract phase: reference = centre of k-space (and not n2/2)
    cp2 = exp(-2*pi*1i*(bsxfun(@times,(-n2/2:(n2/2-1)),(2*ind2-n2))/n2));
    % calculate auto-correlation (full echo sampling)
    sig = bsxfun(@times,sig.*conj(sig(:,end:-1:1,:)),cp2);

    ng = 10;                                     % #points for phase fitting
    iii = floor(n2/2)+(-ng/2:ng/2-1)+1;          % index list
    phi = unwrap(angle(sig(:,iii,:)),[],2);      % linear phase=shift    
    
    A  = [ones(1,ng) ; (-ng/2:(ng/2-1))/ng];     % linear regression kernel
    pinvA = pinv(A);
    del = zeros(size(sig,1),1,nechoes-1);        % delay in pixels
    for le=1:(nechoes-1)
        tmp = phi(:,:,le)*pinvA;
        del(:,:,le) = tmp(:,2);
    end
    
    for le=1:(nechoes-1)
        ppp = del(:,:,le);
        fprintf('echo=%g: delay = %.3g +- %.3g [pts]\n',...
            le+1,mean(ppp),std(ppp));
    end

    if burst_mode>0
        tmp = del(:,:,2:2:end,:);
        delay(1) = -mean(tmp(:))/bw*1d6;         % convert pixel to [us]
        tmp = del(:,:,1:2:end,:);
        delay(2) = -mean(tmp(:))/bw*1d6;
    else
        delay = -mean(del(:))/bw*1d6;            % convert pixel to [us]
    end
    
    clear sig ind2 ng iii phi A pinvA del ndel
end


%% adjust gradient delay time: shift data by linear phase in f-domain
if any(delay~=0)
    if length(delay)>2
        warning('length(delay)(=%g)>2; using only delay(1)=%g',...
            length(delay),delay(1));
        delay = delay(1);
    end
    if ((length(delay)==2) && (burst_mode==0))
        warning('length(delay)==2 && burst_mode==0; using only delay(1)=%g',...
            delay(1));
        delay = delay(1);
    end
    
    if length(delay)==1
        fprintf('Compensating for gradient delay = %g [us]\n',delay);
        delay = delay*1d-6;
        phafu = fftshift(exp(2*pi*1i*(-(n2/2):(n2/2-1))/n2*bw*delay));
    else
        % shift odd/even acquisitions with different delays
        fprintf('Compensating for gradient delays = [%g,%g] [us]\n',delay);
        delay = delay*1d-6;
        phafu1 = fftshift(exp(2*pi*1i*(-(n2/2):(n2/2-1))/n2*bw*delay(1)));
        phafu2 = fftshift(exp(2*pi*1i*(-(n2/2):(n2/2-1))/n2*bw*delay(2)));
        
        phafu = complex(zeros(1,xres,1,nechoes,'single'));
        phafu(1,:,1,1:2:nechoes) = repmat(phafu1,[1 1 1 ceil(nechoes/2)]);
        phafu(1,:,1,2:2:nechoes) = repmat(phafu2,[1 1 1 floor(nechoes/2)]);
    end
    d = fft(bsxfun(@times,ifft(d,[],2),phafu),[],2);
end


%% pre-allocation of memory
if do_single
    bb = complex(zeros(mtx_reco,mtx_reco,mtx_reco,nechoes,nc,n3,'single'));
else
    bb = complex(zeros(mtx_reco,mtx_reco,mtx_reco,nechoes,nc,n3));
end


%% loop through echoes
for le=1:nechoes
    if nechoes>1, fprintf('**********\necho# %g/%g\n',le,nechoes); end
    isoddle = isodd(le);
    
    
    %% k-space trajectory
    if ((le>1) && (burst_mode>0))
        % mode = 2;     % use all data
        mode = 4;       % discard lowres data
    else
        mode = 0;
    end
    [k,dcf,npts_low,npts_high,t,xtra] = ...
        zte_header2kspace(h,true,mode,burst_half,max_waspi);
    if ~isoddle && (burst_mode>1), k = -k; end
    if true  
        % sagital and coronal z-axis inverted
        % product reco corrects it in dicom writing
        switch h.image.plane
            % Scan Plane (different from opplane)
            % 2=axial, 4=sagittal, 8=coronal, 16=oblique
            case 2
            case {4,8}
                fprintf('Inverting z-kspace');
                k(:,:,3) = -k(:,:,3);
                
            case 16, warning('oblique plane');
            otherwise, warning('h.image.plane(=%g) unknown',h.image.plane);
        end
    end
    if isempty(overrecofa)
        k = k*mtx_acq/mtx_reco;   % scale k-space -> Fourier interpolation
    else
        k = k*mtx_acq/mtx_reco*overrecofa;
    end

    %% Andre's centre of k-space correction
    % subtract phase of lowfreq fft'ed spoke
    if ((do_phi0corr>1) && (mode>0))
        if isoddle
            ind2 = h.rdb_hdr.feextra;
        else
            ind2 = xres-h.rdb_hdr.feextra-1;
        end
        switch do_phi0corr
            case 2
                lbs = 4;
                apo = exp(-(((0:(n2-1))-ind2)/lbs*2).^2);
                cphi = ifft(bsxfun(@times,d(:,:,:,le,:,:),apo),[],2);
                cphi = bsxfun(@times,cphi./abs(cphi),exp(-2*pi*1i*(0:(n2-1))/n2*ind2));
                cphi(isnan(cphi)) = 1;
                d(:,:,:,le,:,:) = fft(bsxfun(@times,ifft(d(:,:,:,le,:,:),[],2),conj(cphi)),[],2);
            case 3
                lbs = 4; ng = ceil(n2/lbs)*2;
                apo = exp(-(((0:(n2-1))-ind2)/lbs*2).^2);
                % apo = zeros(1,n2); apo(1,((-3:3)+ind2)) = 1; apo(1,(([-4,4])+ind2)) = 0.5;
                cphi = ifftshift(ifft(bsxfun(@times,d(:,:,:,le,:,:),apo),[],2),2);
                cphi = bsxfun(@times,cphi,exp(-2*pi*1i*((0:(n2-1))*ind2)/n2));
                iii = floor(n2/2)+(-ng/2:ng/2-1);
                phi = unwrap(angle(cphi(:,iii,:,:,:,:)),[],2);
                
                A = [ones(1,ng) ; (-ng/2:(ng/2-1))/ng];
                AA = [ones(1,n2) ; (-n2/2:(n2/2-1))/n2];
                pinvA = pinv(A);
                phi_fit = zeros(n1,n2,n3,1,n5,n6);
                for l5=1:n5
                    for l6=1:n6
                        pp = phi(:,:,1,1,l5,l6)*pinvA;
                        % pp(:,1) = 0;
                        phi_fit(:,:,1,1,l5,l6) = pp*AA;
                    end
                end
                phi_fit(isnan(phi_fit)) = 0;
                d(:,:,:,le,:,:) = fft(bsxfun(@times,...
                    ifft(d(:,:,:,le,:,:),[],2),...
                    exp(-1i*fftshift(phi_fit,2))),[],2);
            case 4
                sig = sqrt(mean(d(:,:,:,le,:,:).*conj(d(:,:,:,le,:,:)),6));
                sig = ifftshift(ifft(fftshift(sig,2),[],2),2);
                sig = sig.*conj(sig(:,end:-1:1,:,:,:,:,:));
                sig = bsxfun(@times,sig,exp(-2*pi*1i*((-n2/2:(n2/2-1))*(2*ind2-n2))/n2));
                ng = 10;
                iii = floor(n2/2)+(-ng/2:ng/2-1)+1;
                phi = unwrap(angle(sig(:,iii,:,:,:,:)),[],2);
                
                A  = [ones(1,ng) ; (-ng/2:(ng/2-1))/ng];
                pinvA = pinv(A);
                ppp = zeros(n1,1,n3,1,n5);
                for l5=1:n5
                    pp = phi(:,:,1,1,l5)*pinvA;
                    ppp(:,:,1,1,l5) = pp(:,2);
                end
                
                fprintf('delay = %.3g +- %.3g [pts] = %.3g +- %.3g [us]\n',...
                    mean(ppp(:)),std(ppp(:)),...
                    mean(ppp(:))*1d6/bw,std(ppp(:))*1d6/bw);
                
                ppp = mean(ppp(:));
                
                phafu = fftshift(exp(bsxfun(@times,...
                    -2*pi*1i*(-(n2/2):(n2/2-1))/n2,ppp)),2);
                d(:,:,:,le,:,:) = fft(bsxfun(@times,...
                    ifft(d(:,:,:,le,:,:),[],2),phafu),[],2);
                
                
        end
        clear apo cphi phi
    end
    
    
    %% reshaping data (matrix)
    % size(d) = [nframes,frame_size,nphases,nechoes,nslices,nreceivers]
    %           nspokes in nframes and nslices
    if do_single
        dd = complex(zeros(n1*n5,n2,nc,n3,'single'));
    else
        dd = complex(zeros(n1*n5,n2,nc,n3));
    end
    for l3=1:n3
        for l5=1:n5
            dd((1:n1)+(l5-1)*n1,:,:,l3)=d(:,:,l3,le,l5,coil_arr);
        end
    end
    
    
    %% Phase correction of burst echoes using centre of k-space
    % echo separation loop already summed, but only delayed by 1 train
    if (((do_phi0corr==1)||(do_phi0corr>3)) && (mode>0))
        fprintf('Applying phase correction to centre of k-space\n');
        if isoddle
            ind2 = h.rdb_hdr.feextra;
        else
            ind2 = xres-h.rdb_hdr.feextra;
        end
        % ind2 = ind2+(-1:1);
        
        cphi = dd(:,ind2,:,:);
        cphi = sum(bsxfun(@rdivide,cphi,sum(abs(cphi),2)),2);
        cphi(isnan(cphi)) = 1;
        dd = bsxfun(@times,dd,conj(cphi));
    end
    
    
    %% T2* correction for even acquisitions with inconsistent TE
    if (do_tecorr && ~isodd(le) && (burst_mode>0))
        fprintf('Applying inconsistent-TE correction to even acquisition\n');
        ind2 = xres-h.rdb_hdr.feextra;      % index to max of spoke

        % average signal of max of spokes
        % RMS along coil combination
        % ind2 = ind2+(-1:1);     % looks identical
        dtmp = sqrt(mean(abs(dd((nspokes_low+1):nspokes,ind2,:)).^2,3));
        % dtmp = mean(dtmp,2);
        ntmp = floor(nspokes_high/spokesPerSeg);
        dtmp = reshape(dtmp(1:(ntmp*spokesPerSeg),1),[spokesPerSeg ntmp]);
        dtmp = mean(dtmp,2); % data reversed -> signal increasing

        % magnitude correction factor
        te_corr = dtmp(end:-1:1)/mean(dtmp);
        te_corr = repmat(te_corr,[ceil(size(dd,1)/spokesPerSeg) 1]);
        te_corr = te_corr(1:size(dd,1),:);
        
        % correct to data
        dd = bsxfun(@times,dd,te_corr);
    end
    
    
    %% sorting + truncating data
    if (~isoddle && (burst_mode>0))
        ind_acq_low = (npts_low:-1:1);
        ind_acq_high = (npts_high:-1:1);
    else
        ind_acq_low = (1:npts_low) + xtra;
        ind_acq_high = ((xres-npts_high+1):xres);
    end
    
    dd_low = dd(1:nspokes_low,ind_acq_low,:,:);
    dd_high = dd((1:nspokes_high)+nspokes_low,ind_acq_high,:,:);
    clear dd
    
    
    %% reshape data into [nc,nk] matrix
    dd = [reshape(permute(dd_low,[3 2 1 4]),[nc,nspokes_low*npts_low,n3]) , ...
        reshape(permute(dd_high,[3 2 1 4]),[nc,nspokes_high*npts_high,n3])];
    clear dd_low dd_high
    
    
    %% misc checks    
    if size(k,2)~=size(dd,2)
        warning('size(k,2)~=size(dd,2)');
        fprintf('\tsize(k) =  '); disp(size(k));
        fprintf('\tsize(dd) = '); disp(size(dd));
    end
    fprintf('sum(dd(:)==0) = %g\n',sum(sum(sum(sum(dd==0)))));
    
    % size of reshaped data dd:
    %  dim1=#coils; dim2=#indexed kspace data
    % if ((lbt>0) || (f0~=0)), t = t.'; t = t(:).'; end
    if lbt>0, dd = ak_apodise(dd,[],[],t,lbt,false); end    % line-broadening
    if f0~=0, dd = bsxfun(@times,exp(1i*2*pi*t*f0),dd); end % shift recon frequency
    
    if ~isempty(shft)
        % fprintf('Off-isocenter (3D) data modulation shft=(%g,%g,%g)\n',...
        %             shft(1),shft(2),shft(3));
        dd = bsxfun(@times,dd,exp(2*pi*1i*k(1,:,1)*shft(1) + ...
            2*pi*1i*k(1,:,2)*shft(2) + 2*pi*1i*k(1,:,3)*shft(3)));
    end
    
    
    %% actual reconstruction
    for lp=1:n3
        if do_single
            if do_homodyne && (burst_half==1) && (mode==0)
                fprintf('Gridding + homodyne reconstruction (single)\n');
                b1 = gridding3dsingle_homodyne(k,dd(:,:,lp),dcf,...
                    [1 1 1]*mtx_reco,1.125,2.5,npts_low*nspokes_low,3);
            else
                fprintf('Gridding reconstruction (single)\n');
                b1 = gridding3dsingle(k,dd(:,:,lp),dcf,...
                    [1 1 1]*mtx_reco,1.125,2.5);
            end
            % b1 = grid_matlab(dd(lc,:),k,mtx_reco);
            
            if lp==n3, clear dd; end
            
            for lc=1:nc
                if do_ori
                    for l3=1:mtx_reco
                        if trnsps(1)==3
                            bb(:,:,l3,le,lc,lp) = rot90((b1(:,:,l3,lc)).',rot90fac(1));
                        else
                            bb(:,:,l3,le,lc,lp) = rot90((b1(:,:,l3,lc)),rot90fac(1));
                        end
                    end
                else
                    bb(:,:,:,le,lc,lp) = b1(:,:,:,lc);
                end
            end
        else
            fprintf('Gridding reconstruction (double)\n');
            for lc=1:nc
                fprintf('**********\nCoil# %g/%g\n',coil_arr(lc),nrcv);
                b1 = gridding3d(k,dd(lc,:,lp),dcf,[1 1 1]*mtx_reco);
                if do_ori
                    for l3=1:mtx_reco
                        if trnsps(1)==3
                            bb(:,:,l3,le,lc,lp) = rot90(fliplr(b1(:,:,l3)).',rot90fac(1));
                        else
                            bb(:,:,l3,le,lc,lp) = rot90(fliplr(b1(:,:,l3)),rot90fac(1));
                        end
                        bb(:,:,l3,le,lc,lp) = fliplr(bb(:,:,l3,le,lc,lp));
                    end
                else
                    bb(:,:,:,le,lc,lp) = b1;
                end
            end
            if lp==n3, clear dd; end
        end
    end
    
    
    %% plotting slices (imagesc_row)
    for lp=1:n3
        if ( (plt==1) || ((plt==2)&&~isempty(fname)) || ...
                (plt==5) || ((plt==6)&&~isempty(fname)) )
            fprintf('Plotting\n');
            tmp = abs(sqrt(sum(bb(:,:,:,le,:,lp).*conj(bb(:,:,:,le,:,lp)),5)));
            if ((plt==2)||(plt==6))
                fid = figure('Visible','off');
            else
                figure;
            end
            imagesc_row(tmp,[],'',true,true);
            colormap gray;
            figstr = sprintf('P%05d Exam%d Series%d: echo=%d',...
                h.image.rawrunnum,h.exam.ex_no,h.series.se_no,le);
            set(gcf,'name',figstr);
            drawnow
            
            if ~isempty(fname)
                print([fname '_e' num2str(le) '_p' num2str(lp) '.png'],'-dpng','-r600','-painters');
            end
            if ((plt==2)||(plt==6)), close(fid); end
        end
        
        if false
            if do_single, k = single(k); t = single(t); dcf = single(dcf); end
            save([fname '_kspace.mat'],'-mat','k','t','dcf');
        end
    end
end


%% coil combine
% size(bbabs)=(mtx,mtx,mtx,#timesteps);
if nrcv>1
    bbabs = abs(sqrt(sum(bb.*conj(bb),5)));
else
    bbabs = abs(bb);
end


%% plotting 3-planes through centre (imagesc_ind3d)
if ( (plt==3) || ((plt==4)&&~isempty(fname)) || ...
       (plt==5) || ((plt==6)&&~isempty(fname)) )
   for lp=1:n3
       if ((plt==4)||(plt==6))
           fid = figure('Visible','off');
       else
           figure;
       end
       imagesc_ind3d(bbabs(:,:,:,:,1,lp),'','','',false,true,'row');
       colormap gray;
       figstr = sprintf('P%05d Exam%d Series%d',...
           h.image.rawrunnum,h.exam.ex_no,h.series.se_no);
       set(gcf,'name',figstr);
       drawnow
       
       if ~isempty(fname)
           print([fname '_p' num2str(lp) '.png'],'-dpng','-r600','-painters');
       end
       if ((plt==4)||(plt==6)), close(fid); end
   end
end


%% saving result as mat file
if ~isempty(fname)
    % bb = single(bb);   % convert from double (64bit) to single (32bit) precision
    fprintf('Saving file to ''%s''\n',fname);
    if save_coils || (nrcv==1)
        fprintf('\tstoring full complex data\n');
        save([fname '.mat'],'-mat','-v7.3',...
            'bb','h','mtx_reco','mtx_acq','delay','lbt','fname','f0',...
            'do_ori','bw','npts_low','npts_high','nrcv',...
            'nspokes','nspokes_low','nspokes_high','rot90fac','trnsps',...
            'xres','xtra','coil_arr');
    else
        if nechoes>1
            fprintf('Calculating coil-averaged phase bp\n');
            bp = im2phase(bb,1);
            fprintf('\tstoring bbabs and bp\n');
        else
            bp = [];
            fprintf('\tstoring bbabs\n');
        end
        save([fname '.mat'],'-mat','-v7.3',...
            'bbabs','bp','h','mtx_reco','mtx_acq','delay','lbt','fname','f0',...
            'do_ori','bw','npts_low','npts_high','nrcv',...
            'nspokes','nspokes_low','nspokes_high','rot90fac','trnsps',...
            'xres','xtra','coil_arr');
    end
end


%% exporting images as dicom into scanner database
if ~isempty(fname) && export_dcm
    % inp.MRAcquisitionType = '3D';
    switch export_dcm
        case 1
            inp.SeriesNumber = h.exam.ex_numseries + 100;
            inp.SeriesDescription  = sprintf('BURZTE;3D%g;E%g;%s',...
                mtx_reco,nechoes,h.series.se_desc);
            inp.EchoTime = round(max(t)*1d5)*1d-2;         % pw_flat [ms]
            write_dicom(fname,bbabs,h,inp);
        case 2     % Orchestra-based dicom writing
            fprintf('Storing dicom images (using Orchestra)\n');
            inp.description  = sprintf('BURZTE;3D%g;E%g;%s',...
                mtx_reco,nechoes,h.series.se_desc);
            inp.SliceThickness = fov/mtx_reco;
            inp.SpacingBetweenSlices = fov/mtx_reco;
            inp.PixelSpacing = [1 1]*fov/mtx_reco;
            inp.SeriesNumber = h.exam.ex_numseries + 200;
            write_dicom_ox(fname,bbabs,h,inp,pfname);
        case 3
            % export as dicom, if fname present fname_template.dcm is found
            template_fname = [fname '_template.dcm'];
            if ~exist(template_fname,'file')
                fprintf('No template dicom file found\n\t%s\n',...
                    template_fname);
            else
                fprintf('Storing dicom images (using Matlab)\n');
                inp.SeriesNumber = h.exam.ex_numseries + 300;
                inp.SeriesDescription  = sprintf('ZTE;3D%g;E%g;%s',...
                    mtx_reco,nechoes,h.series.se_desc);
                inp.SliceThickness = h.image.slthick;
                inp.RepetitionTime = round(h.image.tr*1d-1)*1d-2;  % [ms]
                inp.EchoTime = round(max(t)*1d5)*1d-2;             % pw_flat [ms]
                write_dicom_template(fname,bbabs,template_fname,inp);
            end
        otherwise
            error('export_dcm(=%d)~=1,2,3',export_dcm);
    end
end

fprintf('recon_zte.m finished: ');
fprintf('total reconstruction time = %.2f s\n',toc(timerVal));
if nargout<1, clear bb; end
    
end  % main function recon_zte.m
