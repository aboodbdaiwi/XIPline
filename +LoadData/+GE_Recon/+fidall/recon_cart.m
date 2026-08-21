function [bb,bbabs] = recon_cart(d,h,wfn,zf,delay,lb,fname,comb_time,grdwrp,plt)
%RECON_CART Reconstruct Cartesian readout acquired with fidall
% [bb,bbabs] = recon_cart(d,h,wfn,mtx,delay,lb,fname,comb_time,grdwrp,plt)
%         d   Raw (p-file) data  (or pfile fname)
%         h   Header from p-file (or empty)
%       wfn   Location of waveform .mat file
%        zf   Reconstruction matrix size: 0=keep; <0->multiplier        (0)
%     delay   Gradient-acquisition delay                        [pts]   (0)
%        lb   Apodisation                                               (0)
%     fname   Print <fname>.png and save reco as <fname>.mat           ([])
%             also export dicom if isdeployed
% comb_time   RMS combine time steps into single image                  (0)
%    grdwrp   Apply gradwarp distortion correction (or gradcoil)     (true)
%       plt   (1) plotting: 0=off,1=on  
%             (2) export as dicom          isdeployed->([1 1]) else ([1 0])
%
%        bb  Reconstructed data        (mtx,mtx,mtx,#timesteps,#coils)
%     bbabs  RMS coil-combined data    (mtx,mtx,mtx,#timesteps)
%
%  2/2025 Rolf Schulte
if (nargin<1), help(mfilename); return; end


%% input variables
if ~exist('wfn','var'), wfn = []; end
if ~exist('zf','var'),  zf = []; end
if isempty(zf),         zf = 0; end
if ~exist('delay','var'), delay = []; end
if isempty(delay),      delay = 0; end
if length(delay)~=1, error('length(delay)(=%g)~=1',length(delay)); end
if ~exist('lb','var'),  lb = []; end
if isempty(lb),         lb = 0; end
if ~exist('fname','var'), fname = []; end
if ~isempty(fname)
    if ~islogical(fname)
        if ~isempty(regexpi(fname,'\.7$', 'once')), fname = fname(1:end-2); end
        if ~isempty(regexpi(fname,'\.h5$', 'once')), fname = fname(1:end-3); end
        if ~isempty(regexpi(fname,'\.mat$', 'once')), fname = fname(1:end-4); end
    end
end
if ~exist('comb_time','var'), comb_time = []; end
if isempty(comb_time),  comb_time = 0; end
if ~exist('grdwrp','var'),  grdwrp = []; end
if isempty(grdwrp),         grdwrp = true; end
if ~exist('plt','var'),     plt = []; end
if isempty(plt),            plt = true; end
if length(plt)<2
    if isdeployed,          plt(2) = true;
    else,                   plt(2) = false; end
end

timerVal = tic;                        % record execution time


%% reading in data, if pfile name is given
if ~isnumeric(d) || isempty(d)
    if ~exist(d,'file') && ~isempty(d)
        warning('strange input for d (=''%s''); file not existing?',d); 
    end
    [d,h] = read_p(d);
end


%% reading waveform file
if iscell(wfn),  wfn = wfn{h.image.user3}; end
if isempty(wfn), wfn = get_wfname(h); end
if isstruct(wfn)
    wf = wfn;
else
    if ~isempty(regexpi(wfn,'\.wav$', 'once')), wfn = wfn(1:end-4); end
    if isempty(regexpi(wfn,'\.mat$', 'once')),  wfn = [wfn '.mat']; end
    if ~exist(wfn,'file'), error('file not found: wfn=%s',wfn); end
    wf = load(wfn);          % load waveform
end


%% check fields required for reconstruction
fina = {'k','dcf','ind','t','npix','mtx','fov','sampfac','dim'};
for l=1:length(fina)
    if ~isfield(wf,fina{l}), error('wf.%s not existing',fina{l}); end
end
if isempty(wf.dcf),      error('wf.dcf is empty'); end
if any(isnan(wf.dcf)),   error('wf.dcf contains NaN'); end
if any(isinf(wf.dcf)),   error('wf.dcf contains inf'); end
if wf.sampfac<1, warning('wf.sampfac(=%g)<1; ignoring',wf.sampfac); end
if length(wf.mtx)<3, wf.mtx(3) = 1; end
if ~isfield(wf,'inv_grad'), wf.inv_grad = false; end
if wf.inv_grad
    fprintf('Inverted gradient: d=conj(d)')
    d = conj(d);
end


%% matrix size
dim = wf.dim;                          % dimensions encoded by ak_grad
csi_dims = h.rdb_hdr.csi_dims;         % phase encodes via PSD
xcsi = h.rdb_hdr.xcsi;
ycsi = h.rdb_hdr.ycsi;
zcsi = h.rdb_hdr.zcsi;
do_zf = false;                         % zero-filling
zf = zf(:).';

mtx = wf.mtx;
if csi_dims>0
    mtx(1,2) = ycsi;
    mtx(1,3) = zcsi;
end
% fprintf('zf = '); disp(zf)
% fprintf('mtx = '); disp(mtx)
for ll=1:length(zf)
    if zf(1,ll)>mtx(1,ll)
        mtx(1,ll) = zf(1,ll);
        do_zf = true;
    else
        if zf(1,ll)<-0.01
            mtx(1,ll) = round(abs(zf(1,ll))*mtx(1,ll));
            do_zf = true;
        end
    end
end
% fprintf('mtx = '); disp(mtx)
mtx = [mtx , ones(1,3-dim)];           % add 1's for grad pre-allocation
if any(mtx<1), warning('mtx<1');  end


%% misc variables and checks
fov = h.rdb_hdr.fov*1d-3;                   % actual field-of-view [m]
bw = h.rdb_hdr.user0;                       % full acq bandwidth [Hz]
[nexc,~,nphases,nechoes,nslices,ncoils] = size(d); % data size
if nphases~=1, warning('nphases=size(d,3)(=%g)~=1',nphases); end

sumind2 = sum(wf.ind,2);
sumind2 = sumind2(sumind2>0);
if (any(sumind2~=wf.mtx(1)*wf.sampfac(1)) && (dim>1))
    warning('any(sum(wf.ind,2)(1st=%g)~=wf.mtx(1)(=%g)*wf.sampfac(1)(=%g))',...
        sumind2(1,1),wf.mtx(1),wf.sampfac(1));
end
if ((dim>1) && (csi_dims>0))
    error('(dim(=%g)>1) && (csi_dims(=%g)>0)',dim,csi_dims);
end
if ((dim==1) && (csi_dims==0))
    warning('(dim(=%g)==1) && (csi_dims(=%g)==0)',dim,csi_dims);
end
if xcsi>1, warning('xcsi(=%g)>1',xcsi); end
xloc = h.rdb_hdr.user26;
yloc = h.rdb_hdr.user27;
zloc = h.rdb_hdr.user28;
rot90fac = h.data_acq_tab.rotate;
trnsps = h.data_acq_tab.transpose;
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
end


%% off-isocentre shift
% shft = [];
shft = zeros(1,max([dim csi_dims+1]));
switch max([dim csi_dims+1])
    case 1
        n3 = 1;
    case 2
        n3 = nslices;
        if ((abs(xloc)>0.01) || (abs(yloc)>0.01))
            shft = wf.mtx(1,1:2).*[xloc yloc]*1d-3./fov;
        end
    case 3
        n3 = mtx(3);
        if nslices>1, warning('3D and nslices(=%g)>1',nslices); end
        if ((abs(xloc)>0.01) || (abs(yloc)>0.01) || (abs(zloc)>0.01))
            shft = wf.mtx.*[xloc yloc zloc]*1d-3./fov;
        end
    otherwise
        error('dim(=%g)~=1, 2 or 3',dim);
end
if wf.inv_grad, shft = -shft; end
for ld=1:length(shft)
    % if ~isodd(mtx(ld)), shft(ld) = shft(ld) + 0.5; end
    shft(ld) = shft(ld) + 0.5;
end


%% selecting indexed data + reshaping
fprintf('Selecting indexed data + reshaping\n');
% size(dd)= (mtx(2),nfreq,mtx(3),nrep,ncoils)
dtmp = raw2grid(d,ischop(h),wf.k,shft,[],[],wf.ind,delay*1d-6,bw,false);
[~,ns1,nrep,~] = size(dtmp);
if csi_dims==0
    nk_pha1 = wf.mtx(2);
    switch dim
        case 1, error('csi_dims==0 && dim==1: phase encoding missing');
        case 2, nk_pha2 = nslices;
        case 3, nk_pha2 = wf.mtx(3);
    end
    nk_freq = ns1/(wf.mtx(2)*wf.mtx(3));
    if abs(nk_freq-floor(nk_freq))>1d-10
        warning('ns1(=%g) not multiple of nk_phas(=%g)',ns1,nk_freq);
    end
    
    dd = complex(zeros(1,ns1,nslices,nrep,ncoils,class(d)));
    
    for lr=1:nrep
        for lsli=1:nslices
            for lc=1:ncoils
                dd(1,:,lsli,lr,lc) = dtmp(lc,:,lr,lsli);
            end
        end
    end
else     % phase encoding via csicore
    if nslices>1
        warning('(csi_dims(=%g)>0 && (nslices(=%g)>1',csi_dims,nslices);
    end
    nk_freq = ns1;
    nk_pha1 = ycsi;
    nk_pha2 = zcsi;
    nrep = nechoes;
    
    dd = complex(zeros(1,ns1,nk_pha1*nk_pha2,nslices,nrep,ncoils,class(d)));
    for lr=1:nrep
        for lsli=1:nslices
            for lc=1:ncoils
                ii = ((nk_pha1*nk_pha2):-1:1) + (lr-1)*(nk_pha1*nk_pha2);
                dd(1,:,:,lsli,lr,lc) = dtmp(lc,:,ii,lsli);
            end
        end
    end
    dim = csi_dims+1;
    if csi_dims==1, nk_pha2 = nslices; end
end
clear dtmp;
dd = reshape(dd,[nk_freq nk_pha1 nk_pha2 nrep ncoils]);


%% re-sort data for centre-out phase encoding
if isfield(wf,'ind_pe')
    if ~isempty(wf.ind_pe)
        fprintf('Resorting phase encoding\n');
        if any(size(wf.ind_pe)~=[nk_pha1 nk_pha2 2])
            warning('size(wf.ind_pe)(=[%d %d %d])~=[nk_pha1(=%d) nk_pha2(=%d) 2])',...
                size(wf.ind_pe,1),size(wf.ind_pe,2),size(wf.ind_pe,3),...
                nk_pha1,nk_pha2);
        end
        % dd(:,wf.ind_pe(:,1,1),:,:,:,:) = dd;  % 2D
        dtmp = dd;
        for l1=1:nk_pha1
            for l2=1:nk_pha2
                dd(:,wf.ind_pe(l1,l2,1),wf.ind_pe(l1,l2,2),:,:) = ...
                    dtmp(:,l1,l2,:,:,:);
            end
        end
        clear dtmp
    end
end


%% apodisation
if abs(lb)>1d-10
    fprintf('Apodisation with lb=%g\n',lb);
    lbf1 = lb_fun(nk_freq,nk_freq,[0 lb],'h').';
    lbf2 = lb_fun(nk_pha1,nk_pha1,[0 lb],'h');
    dd = bsxfun(@times,dd,lbf1);
    dd = bsxfun(@times,dd,lbf2);
    if nk_pha2>1
        lbf3 = permute(lb_fun(nk_pha2,nk_pha2,[0 lb],'h'),[1 3 2]);
        dd = bsxfun(@times,dd,lbf3);
    end
end


%% zero-filling (Fourier interpolation)
if do_zf
    switch max([dim csi_dims+1])
        case 2, nzf = [mtx(1)*wf.sampfac,mtx(2),n3,nrep,ncoils];
        case 3, nzf = [mtx(1)*wf.sampfac,mtx(2),mtx(3),nrep,ncoils];
    end
    fprintf('Zero-filling data to %s\n',num2str(nzf));
    dd = truma(dd,false,nzf);
end


%% actual reconstruction
fprintf('Reconstructing data: %dD-ifft\n',dim);
switch dim
    case 2
        bb = fftshift(fftshift(ifft(ifft(...
            ifftshift(ifftshift(dd,1),2),[],1),[],2),1),2);
    case 3
        bb = fftshift(fftshift(fftshift(ifft(ifft(ifft(...
            ifftshift(ifftshift(ifftshift(dd,1),2),3),[],1),[],2),[],3),1),2),3);
end
bb = fliplr(bb);


%% truncate oversampling along frequency encoding direction
if wf.sampfac>1
    fprintf('Truncating oversampling along freq-dir\n');
    bb = truma(bb,false,[mtx(1),mtx(2),n3,nrep,ncoils]);
end


%% gradwarp distortion correction
if ischar(grdwrp)
    gradcoil = grdwrp;
    grdwrp = true;
else
    gradcoil = '';
end
if grdwrp
    try
        bb = grad_warp(bb,h,gradcoil,(dim==3));
    catch ME
        warning(sprintf('grad_warp failed with Error\n\t%s',ME.message));
        grdwrp = -1;
    end
end

%% image orientation
bb = mri_orientation(bb,rot90fac,trnsps);


%% RMS coil combination
if ncoils>1
    bbabs = sqrt(mean(bb.*conj(bb),5));
else
    if ~isreal(bb), bbabs = abs(bb);
    else,           bbabs = bb; end
end


%% RMS time-steps combination
if comb_time>0
    bbabs = sqrt(mean(bbabs.^2,4));
    nrep = 1;
end


%% saving result as mat file
if ~isempty(fname)
    fprintf('Saving to file ''%s''\n',fname);
    bb = single(bb);
    bbabs = single(bbabs);
    save([fname '.mat'],'-mat','-v7.3',...
        'bb','bbabs','h','wfn','mtx','delay','fname','zf','lb',...
        'comb_time','grdwrp');
end


%% plotting
if (plt(1)>0) || ~isempty(fname)
    fprintf('Plotting results\n');
    if (plt(1)>0)
        fid = figure(50); clf;
    else
        fid = figure('visible','off');
    end
    set(fid,'DefaultAxesFontSize',14);
    switch dim
        case 2
            bbplt = reshape(bbabs,[size(bbabs,1),size(bbabs,2),n3,nrep]);
            bbplt = permute(bbplt,[1 2 4 3]);
            rshp = false;
            if (nrep==1) && (nslices>8), rshp = true; end
            if (nrep>8) && (nslices==1), rshp = true; end
            imagesc_row(bbplt,[],[],rshp,true);
        case 3
            rot = [];
            if (h.image.plane==2), rot = [0 1 1]; end
            imagesc_ind3d(bbabs,[],[],[],false,true,'gl3',rot);
    end
    colormap(gray);
    figstr = sprintf('P%05d Exam%d Series%d',...
        h.image.rawrunnum,h.exam.ex_no,h.series.se_no);
    set(fid,'name',figstr);
    drawnow
    if ~isempty(fname)
        print([fname '.png'],'-dpng','-r300');
    end
    if ~(plt(1)>0), close(fid); end
end


%% exporting images as dicom into scanner database
if ~isempty(fname) && (plt(2)>0)
    fprintf('Exporting images to dicom\n');
    inp.SeriesDescription  = sprintf('MNSRP;CART;%dD%g;TS%g;%s',...
        dim,mtx(1),nrep,h.series.se_desc);
    inp.MRAcquisitionType = sprintf('%dD',dim);
    % if any(abs(diff(niso))>1d-5)
    %     spacing = round(h.rdb_hdr.fov/mtx_reco(3)/niso(3)*10)/10;
    %     inp.SpacingBetweenSlices = spacing;
    % end
    % matlab-based dicomwrite.m without template
    % dicom header generated from p-file header
    write_dicom(fname,bbabs,h,inp,[],[],export_dcm,3);
end


%% print execution time
ttv = toc(timerVal);
fprintf('recon_cart: runtime = %g [m] %02g [s]\n',...
    floor(ttv/60),round(mod(ttv,60)));

if nargout<1, clear bb; end


end      % recon_cart.m
