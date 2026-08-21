function [bb,bbabs] = recon_grid3d(d,h,wfn,mtx_reco,delay,lb,spks,...
    fname,f0,plt,export_dcm,coco)
%RECON_GRID3D  3D Gridding reconstruction
%[bb,bbabs] = recon_grid3d(d,h,wfn,mtx_reco,delay,lb,spks,fname,f0,...
%                          plt,export_dcm,coco)
%                                                          [unit] (default)
%         d  Raw data (or P-File/ScanArchive fname)
%         h  Header from p-file (or empty)
%       wfn  Location of waveform .mat file (or wf structure)
%       mtx  Reconstruction matrix size (Fourier interp.)              ([])
%     delay  Gradient-acquisition delay                    [us]         (0)
%            1d or 2d=[WHOLE,ZOOM]
%        lb  Apodisation (Gauss time,spatial,exp time) [Hz,mtx,Hz]([0 0 0])
%      spks  Remove spike noise                                     (false)
%     fname  Print <fname>.png and save reco as <fname>.mat            ([])
%            also export dicom if export_dcm
%        f0  Change in centre frequency                    [Hz]         (0)
%       plt   Plotting: 0=off;1=on                isdeployed->(0) else  (1)
%export_dcm  Save images in dicom format:  isdeployed->(true), else (false)
%      coco  Coil combination via mri_coil_combine.m                   ([])
%
%        bb  Reconstructed data    (mtx(1),mtx(2),mtx(3),#timesteps,#coils)
%     bbabs  RMS coil-combined data (mtx(1),mtx(2),mtx(3),#timesteps)
%
%  5/2024  Rolf Schulte
if (nargin<1), help(mfilename); return; end


%% input variables
verb = true;
RFS_MAX_NSCANS = 16382;                % maximum number of scans
if ~exist('wfn','var'),   wfn = ''; end
if ~exist('mtx_reco','var'), mtx_reco = []; end
if isempty(mtx_reco),     mtx_reco = 0; end
if length(mtx_reco)<2, mtx_reco = repmat(mtx_reco,[1 3]); end
mtx_reco = mtx_reco(:).';
if length(mtx_reco)~=3, error('length(mtx_reco)(=%d)~=3',length(mtx_reco)); end
if ~exist('delay','var'), delay = []; end
if isempty(delay),        delay = 0; end
if ~exist('lb','var'),    lb = []; end
if isempty(lb),           lb = 0; end
if length(lb)<2,          lb = [lb 0]; end
if length(lb)<3,          lb = [lb 0]; end
if ~exist('spks','var'),  spks = []; end
if isempty(spks),         spks = false; end
if ~exist('fname','var'), fname = []; end
if ~isempty(fname)
    if ~islogical(fname)
        if ~isempty(regexpi(fname,'\.7$',  'once')), fname = fname(1:end-2); end
        if ~isempty(regexpi(fname,'\.h5$', 'once')), fname = fname(1:end-3); end
        if ~isempty(regexpi(fname,'\.mat$','once')), fname = fname(1:end-4); end
    end
end
if ~exist('f0','var'),    f0 = []; end
if isempty(f0),           f0 = 0; end
if length(f0)~=1, error('length(f0)~=1'); end
if ~exist('plt','var'),   plt = []; end
if isempty(plt)
    if isdeployed,        plt = false;
    else,                 plt = true; end
end
if ~exist('export_dcm','var'),export_dcm = []; end
if isempty(export_dcm)
    if isdeployed,        export_dcm = true;
    else,                 export_dcm = false; end
end
if length(export_dcm)~=1, error('length(export_dcm)~=1'); end
if ~exist('coco','var'),  coco = []; end
if verb>0, timerVal = tic; end         % record execution time
do_loop = false;        % reconstruct time steps at once


%% reading in data, if pfile name is given
if ~isnumeric(d) || isempty(d)
    if ~exist(d,'file') && ~isempty(d)
        warning('strange input for d (=''%s''); file not existing?',d); 
    end
    [d,h] = read_p(d,true);
end
d = single(d);
if size(d,5)>1, error('#slices(=%d)>1',size(d,5)); end


%% reading waveform file
if iscell(wfn), wfn = wfn{h.image.user3}; end
if isempty(wfn), wfn = get_wfname(h); end
wf = load_waveform(wfn,true,verb);


%% different WHOLE and ZOOM delay times (delay->array)
if length(delay)>1
    switch h.rdb_hdr.grad_mode
        case 1, delay = delay(1);
        case 2, delay = delay(2);
        otherwise
            error('length(delay)(=%g)>1 & h.rdb_hdr.grad_mode(=%g) ~= 1 or 2',...
                length(delay),h.rdb_hdr.grad_mode);
    end
end


%% extract fields required for reconstruction
if ~isfield(wf,'mtx'), error('wf.mtx not existing'); end
if ~isfield(wf,'k'),   error('wf.k not existing'); end
if ~isfield(wf,'dcf'), error('wf.dcf not existing'); end
if ~isfield(wf,'ind'), error('wf.ind not existing'); end
if ~isfield(wf,'t'),   error('wf.t not existing'); end
k = wf.k;
dcf = wf.dcf;
ind = wf.ind;
t = wf.t;
ns = sum(ind(end,:),2);
if size(t,1)==1
    if size(t,2)<=size(ind,2)
        t = repmat(t,[size(ind,1) 1]);
    end
end
t = t.'; t = t(:).';

if isempty(dcf),    error('dcf is empty'); end
if any(isnan(dcf)), error('dcf contains NaN'); end
if any(isinf(dcf)), error('dcf contains inf'); end

mtx_acq = wf.mtx(:).';       % nominal matrix resolution of acquisition
if size(mtx_acq,2)<2, mtx_acq = [mtx_acq mtx_acq]; end
if size(mtx_acq,2)<3, mtx_acq = [mtx_acq mtx_acq(1,end)]; end
niso = [1,1,1];              % non-isotropic factor for off-iso shift
if isfield(wf,'fov')
    if length(wf.fov)==3
        if (any(diff(wf.fov)) || any(diff(wf.mtx)))
            niso = wf.fov(1)./wf.fov;
        end
    end
end
clear wf

if abs(mtx_reco(1))<1
    mtx_reco = mtx_acq(:).';
else
    % k = k*mtx_acq/mtx_reco;
    % scale k-space -> Fourier interpolation
    k = bsxfun(@times,k,shiftdim(mtx_acq./mtx_reco,-1));
end
fov = h.rdb_hdr.fov*1d-3;    % actual field-of-view [m]
k = single(k);
dcf = single(dcf);
t = single(t);


%% reshape data for vap_continue_loop
if ((sum(ind(:))>RFS_MAX_NSCANS) && (size(d,4)>1))
    if verb>0, fprintf('Reshaping data for vap_continue_loop\n'); end
    [n1,n2,n3,n4,n5,n6] = size(d);
    tmp = complex(zeros(n1*n4,n2,n3,1,n5,n6,class(d)));
    for l4=1:n4
        ii = (1:n1)+(l4-1)*n1;
        tmp(ii,:,:,1,:,:) = d(:,:,:,l4,:,:);
    end
    d = tmp;
    clear tmp
end


%% misc variables to checkings
bw = h.rdb_hdr.user0;
[nexc,~,n3,~,~,ncoils] = size(d);
xloc = h.rdb_hdr.user26;
yloc = h.rdb_hdr.user27;
zloc = h.rdb_hdr.user28;
rot90fac = h.data_acq_tab.rotate;
trnsps = h.data_acq_tab.transpose;
if ~((trnsps(1)==0)||(trnsps(1)==3))
    warning('h.data_acq_tab.transpose (=%g) not 0 or 3: unknown image orientation',trnsps(1));
end
if n3~=1, warning('size(d,3)(=%g)~=1',n3); end
if nexc<size(ind,1)
    warning('nexc(=%g)<size(ind,1)(=%g)',nexc,size(ind,1));
end
tmp = nexc/size(ind,1);
if abs(tmp-floor(tmp))>1d-10
    warning('nexc(=%g) not multiple of size(ind,1)(=%g); truncating data',...
        nexc,size(ind,1));
    d = d(1:floor(tmp)*size(ind),:,:,:,:,:);
end


%% pre-processing of data
if spks
    if verb>0, fprintf('Removing spike noise\n'); end
    warning('remove_spikes dangerous for 3D radial: removes signal in k centre');
    d = remove_spikes(d); 
end
if ((abs(xloc)>0.01) || (abs(yloc)>0.01) || (abs(zloc)>0.01))
    shft = mtx_reco.*[xloc yloc zloc]*1d-3/fov.*niso;
else
    shft = [];
end

cart_down_fac = [];
dd = raw2grid(d,ischop(h),k,shft,cart_down_fac,[],ind,delay*1d-6,bw,false);
clear d ind
% size of reshaped data dd:
%  dim1=#coils; dim2=#indexed kspace data; dim3=#excitations; dim4=#slices
if any(abs(lb)>0), dd = ak_apodise(dd,k,lb(2),t,[lb(1) lb(3)],false); end
if abs(f0)>0.01
    if verb>0, fprintf('Off-resonance f0(=%g[Hz]) reconstruction\n',f0); end
    phafu = exp(1i*2*pi*t*f0); 
    dd = bsxfun(@times,dd,phafu);
    clear phafu
end
clear t


%% actual reconstruction
ntimesteps = size(dd,3);
if trnsps(1)==3
    if mod(rot90fac(1),2)==1
        mtx0 = mtx_reco(1,[1,2]);
    else
        mtx0 = mtx_reco(1,[2,1]);
    end
else
    if mod(rot90fac(1),2)<1
        mtx0 = mtx_reco(1,[1,2]);
    else
        mtx0 = mtx_reco(1,[2,1]);
    end
end
if trnsps(1)==3, mtx0 = mtx0(1,[2,1]); end
if ~do_loop
    if verb>0, fprintf('Reconstructing/storing data as single (reshape)\n'); end
    dd2 = reshape(permute(dd,[1 3 2]),[ncoils*ntimesteps,size(k,2)]);
    clear dd
    b1 = gridding3dsingle(k,dd2,dcf,mtx_reco);
    clear dd2
    b1 = reshape(b1,[mtx_reco ncoils ntimesteps]);
    bb = complex(zeros(mtx0(1),mtx0(2),mtx_reco(3),...
        ntimesteps,ncoils,'single'));
    for lt=1:ntimesteps
        for lc=1:ncoils
            % image orientation
            for l3=1:mtx_reco(3)
                if trnsps(1)==3
                    bb(:,:,l3,lt,lc) = rot90(fliplr(b1(:,:,l3,lc,lt)).',rot90fac(1));
                else
                    bb(:,:,l3,lt,lc) = rot90(fliplr(b1(:,:,l3,lc,lt)),rot90fac(1));
                end
            end
        end
    end       % for lt=1:ntimesteps
else          % do_loop
    if verb>0, fprintf('Reconstructing/storing data as single (loop)\n'); end
    bb = complex(zeros(mtx0(1),mtx0(2),mtx_reco(3),...
        ntimesteps,ncoils,'single'));
    for lt=1:ntimesteps
        if verb>0, fprintf('gridding %d/%d\n',lt,ntimesteps); end
        b1 = gridding3dsingle(k,dd(:,:,lt),dcf,mtx_reco);
        for lc=1:ncoils
            % image orientation
            for l3=1:mtx_reco(3)
                if trnsps(1)==3
                    bb(:,:,l3,lt,lc) = rot90(fliplr(b1(:,:,l3,lc)).',rot90fac(1));
                else
                    bb(:,:,l3,lt,lc) = rot90(fliplr(b1(:,:,l3,lc)),rot90fac(1));
                end
            end
        end
    end       % for lt=1:ntimesteps
end
clear dd b1 dd1 vv k dcf


%% coil combine
if ncoils>1
    if ~isempty(coco)
        bb = mri_coil_combine(bb,10,'',coco);
        bbabs = abs(bb);
    else
        bbabs = abs(sqrt(mean(bb.*conj(bb),5)));
    end
else
    bbabs = abs(bb);
end


%% saving result as mat file
if ~isempty(fname)
    if verb>0, fprintf('Saving results to ''%s.mat''\n',fname); end
    bb = single(bb);   % single (32bit) precision
    save([fname '.mat'],'-mat','-v7.3',...
        'bb','h','wfn','mtx_reco','delay','lb','spks','fname',...
        'f0','cart_down_fac');
end


%% plotting: 3 orientations
if (plt>0) || ~isempty(fname)
    if verb>0, fprintf('Plotting results\n'); end
    if (plt>0)
        fid = figure;
    else
        fid = figure('visible','off');
    end
    imagesc_ind3d(bbabs,'','','',false,true,'row');
    colormap gray;
    figstr = sprintf('P%05d Exam%d Series%d',...
        h.image.rawrunnum,h.exam.ex_no,h.series.se_no);
    
    set(fid,'name',figstr);
    drawnow
    
    if ~isempty(fname)
        if isdeployed
            print(fid,[fname '.png'],'-dpng','-r300','-painters');
        else
            print(fid,[fname '.png'],'-dpng','-r300');
        end
    end
    if ~(plt>0), close(fid); end
end



%% exporting images as dicom into scanner database
if ~isempty(fname) && (export_dcm>0)
    if verb, fprintf('Exporting images to dicom\n'); end
    inp.SeriesDescription  = sprintf('MNSRP;3D%g;TS%g;%s',...
        mtx_reco(1),ntimesteps,h.series.se_desc);
    inp.MRAcquisitionType = '3D';
    if any(abs(diff(niso))>1d-5)
        spacing = round(h.rdb_hdr.fov/mtx_reco(3)/niso(3)*10)/10;
        inp.SpacingBetweenSlices = spacing;
    end
    
    % matlab-based dicomwrite.m without template
    % dicom header generated from p-file header
    write_dicom(fname,bbabs,h,inp);
end

if verb>0
    ttv = toc(timerVal);
    fprintf('recon_grid3d: runtime = %g [m] %02g [s]\n',...
        floor(ttv/60),round(mod(ttv,60)));
end
if nargout<1, clear bb; end


end      % recon_grid3d.m

