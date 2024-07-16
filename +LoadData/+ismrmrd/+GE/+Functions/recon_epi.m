function [bb,bbabs] = recon_epi(d,h,wfn,mtx,delay,fname,ref_corr,verb)
%RECON_EPI Reconstruct EPI sampled data
% [bb,bbabs] = recon_epi(d,h,wfn,mtx,delay,fname,ref_corr,verb)
%                                                               (default)
%         d  Raw (p-file) data  (or pfile fname)
%         h  Header from p-file (or empty)
%       wfn  Location of waveform .mat file
%       mtx  Reconstruction matrix size (Fourier interpolation) ([])
%            1 value  => isotropic
%            2 values => 1=freq+2=phase
%     delay  Gradient-acquisition delay                    [us] (0)
%     fname  Print <fname>.png and save reco as <fname>.mat     ([]) 
%            also export dicom if template.dcm is found
%  ref_corr  EPI aliasing correction via reference scan         (1)
%            0=off
%            1=linear phase fit + subtraction
%            2=same as 1 with first ref only
%            3=subtract phase of reference scan (w/o coil comb)
%            4=same as 3 with coil combination of reference            
%      verb  Verbose mode                                       (2)
%            0=off, 1=fprintf, 2=fprintf+plotting, 3=add. debugging
%
%        bb  Reconstructed data (mtx,mtx,#exc,#metabolites,#slice,#coils)
%     bbabs  RMS-coil combined images  (mtx,mtx,#exc,#metabolites,#slice)
%
% 6/2020  Rolf Schulte
if (nargin<1), help(mfilename); return; end

thresh = 0.3;                % threshhold for phase from correlation
lb = 1;                      % apodisation for sub_phafit_calibration
show_ref = true;             % reconstruct and show image of reference scan


%% input variables
if ~exist('mtx','var'),      mtx = []; end
if ~exist('delay','var'),    delay = []; end
if isempty(delay),           delay = 0; end
if length(delay)~=1, error('length(delay)(=%g)~=1',length(delay)); end
if ~exist('fname','var'),    fname = []; end
if ~isempty(fname)
    if ~islogical(fname)
        if ~isempty(regexpi(fname,'\.7$')), fname = fname(1:end-2); end
        if ~isempty(regexpi(fname,'\.h5$')), fname = fname(1:end-3); end
        if ~isempty(regexpi(fname,'\.mat$')), fname = fname(1:end-4); end
    end
end
if ~exist('ref_corr','var'), ref_corr = []; end
if isempty(ref_corr),        ref_corr = 1; end
if ~exist('verb','var'),     verb = []; end
if isempty(verb),            verb = 2; end
dbg = false;
if verb>2, dbg = true; end


%% reading in data, if pfile name is given
if ~isnumeric(d)
    if exist(d,'file')
        [d,h] = read_p(d);
    else
        warning('strange input d/file not existing');
    end
end


%% reading waveform file
if isempty(wfn), error('wfn empty'); end
if isstruct(wfn)
    wf = wfn;
else
    if ~isempty(regexpi(wfn,'\.wav$')), wfn = wfn(1:end-4); end
    if isempty(regexpi(wfn,'\.mat$')),  wfn = [wfn '.mat']; end
    if ~exist(wfn,'file'), error('file not found: wfn=%s',wfn); end
    wf = load(wfn);          % load waveform
end


%% check fields required for reconstruction
fina = {'k','dcf','ind','t','npix','mtx','fov','rampsamp','ref','do_burst'};
for l=1:length(fina)
    if ~isfield(wf,fina{l}), error('wf.%s not existing',fina{l}); end
end
if isempty(wf.dcf),      error('wf.dcf is empty'); end
if any(isnan(wf.dcf)),   error('wf.dcf contains NaN'); end
if any(isinf(wf.dcf)),   error('wf.dcf contains inf'); end


%% set reference correction
do_ref = false;
if ref_corr>0
    if (wf.ref(1)>0)
        do_ref = true;
        if verb, fprintf('Applying reference correction (=%g)\n',ref_corr); end
    else
        if verb, fprintf('no reference scan acquired\n'); end
    end
else
    if verb, fprintf('reference correction switched off\n'); end
end
if do_ref
    if length(wf.ref)~=4
        warning('length(wf.ref)(=%g)~=4; adding zeros',length(wf.ref));
        ref = wf.ref(:).';
        wf.ref = [ref zeros(1,4-size(ref,2))];
    end
    if wf.ref(2)<1, warning('wf.ref(2)(=%g) <1',wf.ref(2)); end
end


%% matrix size
mtx_acq = wf.mtx;            % nominal matrix resolution of acquisition
if any(mtx<1), fprintf('mtx<1: setting to default\n'); mtx = []; end
if isempty(mtx)
    mtx = mtx_acq;
else
    switch length(mtx)
        case 1, mtx = [mtx mtx];
        case 2, mtx = mtx(:).';
        otherwise
            warning('length(mtx)(=%g)~=1 or 2',length(mtx));
    end
end
if any(mtx < mtx_acq)
    warning('mtx(=%g %g)<mtx_acq(=%g %g)',mtx,mtx_acq);
end


%% misc variables and checks
fov = h.rdb_hdr.fov*1d-3;    % actual field-of-view [m]
bw = h.rdb_hdr.user0;        % full acquisition bandwidth [Hz]
[nexc,npts,n3,n4,nslices,ncoils] = size(d);
if n3~=1, warning('size(d,3)(=%g)~=1',n3); end
if n4~=1, warning('size(d,4)(=%g)~=1',n4); end

xloc = h.rdb_hdr.user26;     % x off-isocentre shift [mm]
yloc = h.rdb_hdr.user27;     % y off-isocentre shift [mm]
% zloc = h.rdb_hdr.user28;     % z off-isocentre shift [mm]
rot90fac = h.data_acq_tab.rotate;      % image rotate orientation
trnsps = h.data_acq_tab.transpose;     % image translate orientation
if ~((trnsps(1)==0) || (trnsps(1)==3))
    warning('h.data_acq_tab.transpose (=%g) not 0 or 3',trnsps(1));
    fprintf('\tunknown image orientation');
end
if nexc<size(wf.ind,1)
    warning('nexc(=%g)<size(wf.ind,1)(=%g)',nexc,size(wf.ind,1));
end
tmp = nexc/size(wf.ind,1);
if abs(tmp-floor(tmp))>1d-10
    warning('nexc(=%g) not multiple of size(ind,1)(=%g); truncating data',...
        nexc,size(wf.ind,1));
    d = d(1:floor(tmp)*size(wf.ind,1),:,:,:,:,:);
    nexc = size(d,1);
end
nintl = size(wf.ind,1);                % #interleaves
if wf.do_burst
    nintl = wf.nexc;
    if verb>0, fprintf('BURST\n'); end
end

if verb>0
    if wf.rampsamp, fprintf('Ramp sampling; gridding along dim 1\n'); end
    fprintf('nintl = %g\n',nintl);
end


%% pre-processing of data
% off-isocentre shift
if ((abs(xloc)>0.01) || (abs(yloc)>0.01))
    shft = mtx.*[xloc yloc]*1d-3./fov;
else
    shft = [];
end
% shft = []; warning('shft switched off');
if verb>0
    if isempty(shft)
        fprintf('shift = []\n');
    else
        fprintf('shift = (%g,%g)\n',shft);
    end
end
if do_ref
    nrep = floor(nexc/nintl);
    if ref_corr==2, wf.ref(4)=0; end
    if wf.ref(4)>0
        ind_ref = [wf.ref(1):wf.ref(2):nrep ; ...
            wf.ref(4):wf.ref(2):nrep];
        ind_ref = ind_ref(:).';
        if abs(length(ind_ref)-floor(length(ind_ref)))>1d-10
            warning('ind_ref is off');
        end
    else
        ind_ref = wf.ref(1):wf.ref(2):nrep; % index to reference frame(s)        
    end
    nref = length(ind_ref);                 % #ref frames
end


%% selecting indexed data + reshaping
% size(dd)= (nk_phas,nk_freq,nslices,nrep,ncoils)
if verb>0, fprintf('Preprocessing data\n'); end
% dtmp = raw2grid(d,ischop(h),wf.k,shft,[],[],wf.ind,delay*1d-6,bw,false);
if do_ref
    % no shift in phase dir of reference
    dtmp = raw2grid(d,ischop(h),wf.k,shft,[],[],wf.ind,delay*1d-6,bw,false);
    if (nintl>1) && ~wf.do_burst
        ind_ref1 = bsxfun(@plus,(ind_ref-1)*nintl,(1:nintl).');
        ind_ref1 = ind_ref1(:).';
    else
        ind_ref1 = ind_ref;
    end
    dtmp_ref = raw2grid(d(ind_ref1,:,:,:,:,:),ischop(h),wf.k,[shft(1) 0],...
        [],[],wf.ind,delay*1d-6,bw,false);
    dtmp(:,:,ind_ref,:) = dtmp_ref;
else
    dtmp = raw2grid(d,ischop(h),wf.k,shft,[],[],wf.ind,delay*1d-6,bw,false);
end

[~,ns1,nrep,~] = size(dtmp);
nk_phas = wf.mtx(2);                   % #phase encodes
nk_freq = ns1/nk_phas;                 % #freq encoding points
if abs(nk_freq-floor(nk_freq))>1d-10
    error('ns1(=%g) not multiple of nk_phas(=%g)',ns1,nk_freq);
end

dd = complex(zeros(1,ns1,nslices,nrep,ncoils,class(d)));   % initialise
for lr=1:nrep                % reshape
    for lsli=1:nslices
        for lc=1:ncoils
            dd(1,:,lsli,lr,lc) = dtmp(lc,:,lr,lsli);
        end
    end
end
clear dtmp;

dd = reshape(dd,[nk_freq nk_phas nslices nrep ncoils]);    % reshape
% flip every other phase encode (bi-directional epi)
dd(:,2:2:end,:,:,:) = dd(end:-1:1,2:2:end,:,:,:);
if nintl>1
    ddtmp = complex(zeros(size(dd),class(d)));
    for l=1:nintl
        ii2 = (1:nk_phas/nintl)+(l-1)*nk_phas/nintl;
        ii1 = (l:nintl:nk_phas);
        ddtmp(:,ii1,:,:,:,:) = dd(:,ii2,:,:,:,:);
    end
    dd = ddtmp;
    clear ddtmp
end


%% extract + line-broaden reference frames
if do_ref
    ddref = dd(:,:,:,ind_ref,:);            % reference scan (no PE)
    if wf.ref(1)>nrep
        warning('wf.ref(1)(=%g) > nrep(=%g)',wf.ref(1),nrep);
    end
    if wf.ref(4)>nrep
        warning('wf.ref(4)(=%g) > nrep(=%g)',wf.ref(4),nrep);
    end
end
if do_ref && ~show_ref
    if verb>0, fprintf('Removing reference images from data\n'); end
    ind_noref = true(1,nrep);
    ind_noref(1,ind_ref) = false;
    dd = dd(:,:,:,ind_noref,:);
    nrep = nrep-nref;
else
    if verb>0
        fprintf('Reconstructing and displaying reference images\n');
    end
end


%% reconstruction along frequency encoding direction (dim1 of dd)
if verb>0, fprintf('Reconstructing frequency-encoding direction\n'); end
if wf.rampsamp
    k1 = wf.k(1,1:nk_freq,1);               % k-space trajectory along freq
    if mtx(1)~=mtx_acq(1)
        k1 = k1*mtx_acq(1)/mtx(1);          % scale for zero-filling
    end 
    x1 = k1.';
    dcf1 = wf.dcf(1,1:nk_freq);             % density compensation function
    bd = sub_nuift1(dd,mtx(1),k1,dcf1);
    if do_ref
        if lb>0   % apodisation
            fil = exp(-(pi*lb*x1).^2/(4*log(2)));
            ddref = bsxfun(@times,ddref,fil);
        end
        bdref = sub_nuift1(ddref,mtx(1),k1,dcf1);
        if false
            figure(199); clf;
            subplot(2,2,1); imagesc(abs(ddref(:,:,3)))
            subplot(2,2,2); imagesc(angle(ddref(:,:,3)))
            subplot(2,2,3); imagesc(abs(bdref(:,:,3)));
            subplot(2,2,4); imagesc(angle(bdref(:,:,3)));
        end
    end
else
    dd = dd(1:mtx_acq(1),:,:,:,:);          % fft skips symmetric point
    if do_ref, ddref = ddref(1:mtx_acq(1),:,:,:,:); end
    if any(mtx(1)~=mtx_acq(1))              % Fourier interpolation
        dd = truma(dd,false,[mtx(1) mtx_acq(2) nslices nrep ncoils]);
        if do_ref
            ddref = truma(ddref,false,[mtx(1) mtx_acq(2) nslices nref ncoils]);
        end
    end
    x1 = (-mtx(1)/2:mtx(1)/2-1).'/mtx(1);
    bd = fftshift(ifft(ifftshift(dd,1),[],1),1);
    if do_ref
        if lb>0   % apodisation
            fil = exp(-(pi*lb*x1).^2/(4*log(2)));
            ddref = bsxfun(@times,ddref,fil);
        end
        bdref = fftshift(ifft(ifftshift(ddref,1),[],1),1);
    end
end


%% epi echo correction via reference scan
if do_ref
    if verb>0, fprintf('Performing reference correction\n'); end
    switch ref_corr
        case 0     % skip echo separation
        case 1     % linear phase fit + subtraction
            bd = sub_phafit_calibration(bd,bdref,mtx,thresh,3,wf.ref,dbg);
        case 2     % null ref(4); EPI correction with single reference
            bd = sub_phafit_calibration(bd,bdref,mtx,thresh,3,wf.ref,dbg);
        case {3,4}     % subtract phase of reference scan
            bdphase = mean(bdref,4);             % mean over all references
            % weighted coil combination
            if (ncoils>1)&&(ref_corr==4)
                bdphase = bsxfun(@rdivide,bdphase,sum(abs(bdphase),5));
                bdphase = sum(bdphase,5); 
            end
            bdphase = bdphase./abs(bdphase);
            bd = bsxfun(@rdivide,bd,bdphase);
        
        otherwise, warning('ref_corr(=%g) unknown; skipping',ref_corr);
    end
end


%% reconstruction along phase encoding direction (dim2 of dd/bd)
if verb>0, fprintf('Reconstructing phase-encoding direction\n'); end
if any(mtx(2)~=mtx_acq(2))                       % Fourier interpolation
    bd = truma(bd,false,[mtx_acq(1) mtx(2) nslices nrep ncoils]);
end
bb = fftshift(ifft(ifftshift(bd,2),[],2),2);     % actual Fourier trafo


%% image orientation
if verb>0, fprintf('Adjusting image orientation\n'); end
if ~isodd(rot90fac(1))
    ntmp = [mtx(1),mtx(2),nslices,nrep,ncoils];
else
    ntmp = [mtx(2),mtx(1),nslices,nrep,ncoils];
end
bbtmp = complex(zeros(ntmp,class(d)));
for lt=1:nrep
    for lc=1:ncoils
        for l3=1:nslices
            if trnsps(1)==3
                bbtmp(:,:,l3,lt,lc) = rot90(permute(fliplr(bb(:,:,l3,lt,lc)),...
                    [2 1 3]),rot90fac(1));
            else
                bbtmp(:,:,l3,lt,lc) = rot90(fliplr(bb(:,:,l3,lt,lc)),rot90fac(1));
            end
        end
    end
end
bb = bbtmp;
clear bbtmp


%% RMS coil combination
if ncoils>1
    bbabs = sqrt(mean(bb.*conj(bb),5));
else
    bbabs = abs(bb);
end


%% saving result as mat file
if ~isempty(fname)
    if verb>0, fprintf('Saving to file ''%s''\n',fname); end
    bbabs = single(bbabs);   % convert from double (64bit) to single (32bit) precision
    rmflds = {'k','dcf','t','out','ind'};
    for l=1:length(rmflds), wf = rmfield(wf,rmflds{l}); end
    
    save([fname '.mat'],'-mat',...
        'bbabs','h','wf','mtx','mtx_acq','delay','fname',...
        'ref_corr','verb','shft','nintl');
end


%% plotting
if verb>1
    fprintf('Plotting\n');
    figure; set(gcf,'DefaultAxesFontSize',14);
    bbplt = reshape(bbabs,[size(bbabs,1),size(bbabs,2),nslices,nrep]);
    bbplt = permute(bbplt,[1 2 4 3]);
    imagesc_row(bbplt,[],'ind',false,true);
    colormap(gray);
    figstr = sprintf('P%05d Exam%d Series%d',...
        h.image.rawrunnum,h.exam.ex_no,h.series.se_no);
    set(gcf,'name',figstr);
    if ~isempty(fname)
        print(fname,'-dpng','-r600');
    end
end


%% exporting images as dicom into scanner database
if ~isempty(fname)
    fprintf('Storing dicom images\n');
    tmp = sprintf('MNSRP;EPI;SL%g,TS%g',nslices,nrep);
    tmp = [tmp ';' h.series.se_desc];
    inp.SeriesDescription = tmp;
    write_dicom(fname,bbabs,h,inp);
end
if false % ~isempty(fname)
    % export as dicom, if fname present fname_template.dcm is found
    template_fname = [fname '_template.dcm'];
    if ~exist(template_fname,'file')
        fprintf('No template dicom file found\n\t%s\n',template_fname);
    else
        fprintf('Storing dicom images\n');
        inp.SeriesNumber = h.exam.ex_numseries;
        tmp = sprintf('MNSRP;EPI;SL%g,TS%g',nslices,nrep);
        
        tmp = [tmp ';' h.series.se_desc];
        inp.SeriesDescription = tmp;
        inp.MRAcquisitionType = '2D';
        inp.SliceThickness = h.image.slthick;
        inp.RepetitionTime = round(h.image.tr*1d-3);
        inp.EchoTime = round(h.rdb_hdr.user25*1d-2)*1d-1;
        write_dicom_template(fname,bbabs,template_fname,inp);
    end
end

end      % main function recon_epi.m


%% sub-functions
% ********************************************************************

%% 1D non-uniform (slow) Fourier transform along dim1
function b = sub_nuift1(d,n,k,dcf)

nk = size(k,2);                   % #k-space points
nn = size(d);                     % final data matrix size
x1 = 1i*2*pi*((-n/2:(n/2-1))).';  % Fourier kernel
% if nn(2)~=n, x1 = x1/n*nn(2); end % scale or Fourier interpolation
nn(1) = n;
b = complex(zeros(nn,class(d)));


% actual slow Fourier transform
for lk=1:nk
    b = b + bsxfun(@times,d(lk,:,:,:,:,:,:)*dcf(1,lk),exp(k(1,lk)*x1));
end

b = b/n;

end      % sub-function sub_nuift2


%% linear fit to phase difference + subtraction
function bd = sub_phafit_calibration(bd,bdref,mtx,thresh,order,ref,dbg)
%     bd  data (dim1=time,dim2=freq domain)
%         size: [nk_freq nk_phas nslices nrep ncoils]
%  bdref  reference data
%    mtx  Reconstruction matrix size
% thresh  Threshold for masking
%  order  Apply linear phase fit to data: 1=0th,2=1st,3=0th+1st order
%    ref  Reference scans
%    dbg  Debug flag (plotting)

% matrix sizes
[nk_freq,nk_phas,nslices,nrep,ncoils] = size(bd);
nref = size(bdref,4);


%% extract mask (location of object) for fitting
bmag = abs(sqrt(mean(bdref.*conj(bdref),5)));    % RMS coil combination
bmag = mean(mean(bmag,2),4);      % sums up noise, but abs(sum()) not robust
bmag = bsxfun(@minus,bmag,min(bmag,[],1));       % reduce impact of noise
bmag = bsxfun(@rdivide,bmag,max(bmag,[],1));     % normalisation to rows
% mask = regiongrowing1d(bmag,1,thresh);         % calculate mask around max
mask = bmag>thresh;

%% phase calculation
if ncoils>1   % weighted coil combination: normalise first
    bdref = bsxfun(@rdivide,bdref,sum(abs(bdref),5));
end
if false
    % shift half a line to match time of odd echoes
    zf = 4*nk_phas;
    phafu2 = exp(-1i*2*pi*(-zf/2:zf/2-1)/zf/2);
    tmp = fft(fftshift(bsxfun(@times,ifftshift(ifft(bdref(:,2:2:end,:,:,:),zf,2),2),phafu2),2),[],2);
    bdref(:,2:2:end,:,:,:) = tmp(:,1:nk_phas/2,:,:,:);
end

% phase difference even-odd frame
if ref(4)>0
    if nref>2, warning('nref(=%g)>2; not yet implemented',nref); end
    bddiff = bsxfun(@times,bdref(:,:,:,2,:),conj(bdref(:,:,:,1,:)));    
else
    if nref>1, warning('nref(=%g)>1; not yet implemented',nref); end
    bddiff = bsxfun(@times,bdref(:,2:end,:,:,:),conj(bdref(:,1:end-1,:,:,:)));
end
if ncoils>1, bddiff = sum(bddiff,5); end    % actual coil combination
phi = unwrap(angle(bddiff),[],1);           % phase


%% initialise + fit coefficients
phafit = zeros(2,size(bddiff,2),nslices,1); % dim1: 1=0th & 2=1st order phase
x1 = (-mtx(1)/2:mtx(1)/2-1).'/mtx(1);
A  = [ones(mtx(1),1) , x1];                 % lin fit kernel


%% actual phase fit
for lsli=1:nslices
    m1 = mask(:,:,lsli);
    if sum(m1)<8, warning('sum(m1)(=%g)<8) (lsli=%g)\n',sum(m1),lsli); end
    for lref=1 %:nref
        pinvA = pinv(A(m1,:));
        phafit(:,:,lsli,lref) = pinvA*phi(m1,:,lsli,lref);        
    end
end


%% time shift data along phase encodes
if ref(4)>0
    % phafit = phafit/2;
    phafit(2,:,:,:) = phafit(2,:,:,:)/2;
    phafit(1,:,:,:) = cumsum(phafit(1,:,:,:),2);
else
    tmp = zeros(2,nk_phas,nslices,1);
    tmp(:,2:end,:,:) = phafit;
    phafit = cumsum(tmp,2);
end

switch order
    case 0
        warning('doing nothing'); 
        phafu = ones(1,nk_phas,nslices,nref);
    case 1    % 0th order phase correction
        phafu = exp(-1i*phafit(1,:,:,:));
    case 2    % 1st order phase correction
        phafu = exp(-1i*bsxfun(@times,x1,phafit(2,:,:,:)));
    case 3    % 0th + 1st order phase correction
        phafu = exp(-1i*(repmat(phafit(1,:,:,:),[mtx(1) 1 1 1]) + ...
            bsxfun(@times,x1,phafit(2,:,:,:))));
    otherwise, error('order(=%g)~=1,2 or 3');
end


%% actual phase correction of data (in frequency (=spatial) domain in dim1)
if dbg, bd1 = bd; end
bd = bsxfun(@times,bd,phafu);


%% plotting for debugging
if dbg
    figure(102); clf;

    xx = (1:mtx(1));
    for lsli=1:nslices
        m1 = mask(:,:,lsli,1,1);
        subplot(nslices,2,2*lsli);
        plot(xx,phi(:,:,lsli,lref),'-',xx(m1),phi(m1,:,lsli,lref),'x');
        axis tight
        if nslices>1
            set(gca,'xticklabel','');%,'yticklabel','');
            pos = get(gca,'Position');
            % fprintf('%g %g %g\n',pos(2),pos(2)-0.05,(nslices-lsli)/nslices);
            set(gca,'Position',[0.55 pos(2) 0.45 0.8/nslices]);
        end
        
        subplot(nslices,2,2*(lsli-1)+1);
        % plot(xx,abs(bddiff(:,:,lsli,lref)),'-',...
        %     xx(m1),abs(bddiff(m1,:,lsli,lref)),'x');
        plot(xx,bmag(:,:,lsli,lref),'-',...
            xx(m1),bmag(m1,:,lsli,lref),'x');
        
        axis tight
        if nslices>1
            set(gca,'xticklabel','','yticklabel','');
            pos = get(gca,'Position');
            set(gca,'Position',[0 pos(2) 0.5 0.8/nslices]);
        end
    end
    
    figure(104); clf;
    tmp = zeros(nk_freq,nk_phas,nslices,2);
    tmp(:,:,:,1) = bd1(:,:,:,ref(1));
    tmp(:,:,:,2) = bd(:,:,:,ref(1));
    tmp = permute(tmp,[1 2 4 3]);
    imagesc_row(angle(tmp),[-pi pi],'','',true);
    figure(105); imagesc_row(abs(tmp),'','ind','',true);

    figure(106); clf
    subplot(2,1,1);
    plot(squeeze(mod(phafit(1,:,:,1)+pi,2*pi)-pi));
    axis tight
    title('phafit 0-term');
    
    subplot(2,1,2);
    plot(squeeze(phafit(2,:,:,1)));
    axis tight; legend('show');
    title('phafit lin-term');
    
end

end      % sub-function sub_phafit_calibration
