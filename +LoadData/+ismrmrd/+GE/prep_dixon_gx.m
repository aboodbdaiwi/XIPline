
function prep_dixon_gx(d,h,wfnpath,wfn_freq,fname,do_gascor,mtx_reco,delay,lb,spks,...
    f0,nskip)
% Preprocessing 3D Radial GX for offline recon and analysis
%    prep_dixon_gx(d,h,wfn,mtx_reco,delay,lb,spks,...
%                          fname,f0,do_single,do_gascor)
%                                                          [unit] (default)
%         d  Raw data (or P-File/ScanArchive fname)
%         h  Header from p-file (or empty)
%       wfn  Location of waveform .mat file (or wf structure)
%       mtx  Reconstruction matrix resolution (Fourier interp.)   ([])
%     delay  Gradient-acquisition delay                    [us]   (0)
%            1d or 2d=[WHOLE,ZOOM]
%        lb  Gaussian linebroadening (time-3D spatial) [Hz,mtx]   ([0 0])
%      spks  Remove spike noise                                   (false)
%     fname  Print <fname>.png and save reco as <fname>.mat       ([]) 
%            also export dicom if export_dcm
%        f0  Change in centre frequency                    [Hz]   (0)
%     nskip  Number of views to skip                              (0)
% do_gascor  Apply gas contamination correction                   (false)
%
% 1/2022  Andrew Hahn
if (nargin<1), help(mfilename); return; end


%% input variables
verb = true;
RFS_MAX_NSCANS = 16382;                % maximum number of scans
if ~exist('mtx_reco','var'), mtx_reco = []; end
if isempty(mtx_reco),     mtx_reco = 0; end
if ~exist('delay','var'), delay = []; end
if isempty(delay),        delay = 0; end
if ~exist('lb','var'),    lb = []; end
if isempty(lb),           lb = [0 0]; end
lb = lb(:).';
switch length(lb)
    case 1, lb = [lb(1) 0];
    case 2
    otherwise, warning('length(lb)(=%d)~=2',length(lb));
end
if ~exist('spks','var'),  spks = []; end
if isempty(spks),         spks = false; end
% if ~exist('fname','var'), fname = []; end
% if ~isempty(fname)
%     if ~islogical(fname)
%         if ~isempty(regexpi(fname,'\.7$')), fname = fname(1:end-2); end
%         if ~isempty(regexpi(fname,'\.h5$')), fname = fname(1:end-3); end
%         if ~isempty(regexpi(fname,'\.mat$')), fname = fname(1:end-4); end
%     end
% end
if ~exist('f0','var'),    f0 = []; end
if isempty(f0),           f0 = 0; end
if length(f0)~=1, error('length(f0)~=1'); end
if ~exist('nskip','var'), nskip = []; end
if isempty(nskip), nskip = 0; end
if ~exist('do_gascor','var'), do_gascor = []; end
if isempty(do_gascor),    do_gascor = false; end
if verb>0, timerVal = tic; end         % record execution time


%% reading in data, if pfile name is given
if ~isnumeric(d)
    if exist(d,'file')
        [d,h] = LoadData.ismrmrd.GE.Functions.read_p(d,false);
    else
        warning('strange input d/file not existing');
    end
end

d(1:2*nskip,:) = [];

%% reading waveform file: reading waveform file: reading frequencies and flip angles
filePatterns = {'*.mat', '*freq.fdl'};
filePaths = {};
for i = 1:length(filePatterns)
    fileList = dir(fullfile(wfnpath, filePatterns{i}));
    if isempty(fileList)
        fprintf('No file matching %s found in the directory.\n', filePatterns{i});
    else
        % Append each matching file's full path to the list
        for j = 1:length(fileList)
            filePaths{end+1} = fullfile(wfnpath, fileList(j).name);
        end
    end
end

wfn = filePaths{1};
% if isempty(wfn), error('wfn empty'); end
% if isstruct(wfn)
%     wf = wfn;
% else
%     if ~isempty(regexpi(wfn,'\.wav$')), wfn = wfn(1:end-4); end
%     if isempty(regexpi(wfn,'\.mat$')),  wfn = [wfn '.mat']; end
%     if ~exist(wfn,'file'), error('file not found: wfn=%s',wfn); end
%     wf = load(wfn);          % load waveform
% end
wf = LoadData.ismrmrd.GE.Functions.load_waveform(wfn);

% Display the full paths
if isempty(filePaths)
    disp('No matching files found.');
else
    fprintf('Matching files found:\n');
    disp(filePaths');
end
fdl_freq_file = filePaths{2};
freq_array = LoadData.ismrmrd.GE.Functions.read_fdl(fdl_freq_file);
freq_array(1:2*nskip) = [];
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
if ~isfield(wf,'k'),   error('wf.k not existing'); end
if ~isfield(wf,'dcf'), error('wf.dcf not existing'); end
if ~isfield(wf,'ind'), error('wf.ind not existing'); end
if ~isfield(wf,'t'),   error('wf.t not existing'); end

k = wf.k;
dcf = wf.dcf;
ind = wf.ind; ind(1:2*nskip,:) = [];
t = wf.t; 

if size(t,1)==1
    if size(t,2)<=size(ind,2)
        t = repmat(t,[size(ind,1) 1]);
    end
end
t = t.'; t = t(:).';

ppv = size(k,2)/(size(ind,1) + 2*nskip);   % kspace points per view

k = reshape(k, [ppv (size(ind,1) + 2*nskip) size(k,3)]);
k(:,1:2*nskip,:) = [];
k = reshape(k, [1 ppv*size(ind,1) size(k,3)]);

dcf = reshape(dcf, [ppv (size(ind,1) + 2*nskip) size(dcf,3)]);
dcf(:,1:2*nskip,:) = [];
dcf = reshape(dcf, [1 ppv*size(ind,1) size(dcf,3)]);

t = reshape(t, [ppv (size(ind,1) + 2*nskip) size(t,3)]);
t(:,1:2*nskip,:) = [];
t = reshape(t, [1 ppv*size(ind,1) size(t,3)]);

tmp_k = permute(reshape(k, [ppv size(ind,1) size(k,3)]), [3 2 1]);
tmp_k = ipermute(reshape(tmp_k(:,ind(:,1:ppv)), size(tmp_k,1), size(ind,1), sum(ind(1,:))), [3 2 1]);
k = reshape(tmp_k(:,1:2:end,:), [1 size(tmp_k,1)*size(ind,1)/2 size(tmp_k,3)]);

tmp_dcf = permute(reshape(dcf, [ppv size(ind,1) size(dcf,3)]), [3 2 1]);
tmp_dcf = ipermute(reshape(tmp_dcf(:,ind(:,1:ppv)), size(tmp_dcf,1), size(ind,1), sum(ind(1,:))), [3 2 1]);
dcf = reshape(tmp_dcf(:,1:2:end,:), [1 size(tmp_dcf,1)*size(ind,1)/2 size(tmp_dcf,3)]);

tmp_t = permute(reshape(t, [ppv size(ind,1) size(t,3)]), [3 2 1]);
tmp_t = ipermute(reshape(tmp_t(:,ind(:,1:ppv)), size(tmp_t,1), size(ind,1), sum(ind(1,:))), [3 2 1]);
t = reshape(tmp_t(:,1:2:end,:), [1 size(tmp_t,1)*size(ind,1)/2 size(tmp_t,3)]);

clear tmp_k tmp_dcf tmp_t

if isempty(dcf),    error('dcf is empty'); end
if any(isnan(dcf)), error('dcf contains NaN'); end
if any(isinf(dcf)), error('dcf contains inf'); end

mtx_acq = wf.mtx(1);         % nominal matrix resolution of acquisition
clear wf

if mtx_reco <= 0, mtx_reco = mtx_acq; end
k = k*mtx_acq/mtx_reco;      % scale k-space -> Fourier interpolation
fov = h.rdb_hdr.fov*1d-3;    % actual field-of-view [m]


k = double(k);
t = double(t);
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
    d = remove_spikes(d,[],false); 
end
if ((abs(xloc)>0.01) || (abs(yloc)>0.01) || (abs(zloc)>0.01))
    shft = mtx_reco*[xloc yloc zloc]*1d-3/fov;
else
    shft = [];
end

%% Gas contaminant removal
maxk = max(abs(k(:)));
if  maxk > 5
    k2 = (k./maxk).*0.5;
    k = k2;
end
if do_gascor
    dpv = freq_array ~= 0;
    gpv = freq_array == 0;
    
    t_thresh = 0.0015;	% adjust as necessary, this typically performs pretty well though
    demod_freq = -unique(freq_array(dpv));
    % use the data collected after t_thresh ms for correction
    tC = double((1:size(d,2)).*t(1));   % time for fitting 
    tF = tC(tC >= t_thresh);     % time for correction
    gas_views = d(gpv,tC >= t_thresh);
    dissolved_views = d(dpv,tC >= t_thresh);
    
    % nonlinear optimization options
    fitOptions = optimoptions('lsqcurvefit');
    fitOptions.Display = 'iter';
    fitOptions.MaxIter = 10000;
    fitOptions.TolFun=1E-16;
    fitOptions.TolX = 1E-15;
    fitOptions.FinDiffType = 'central';
    fitOptions.Algorithm = 'trust-region-reflective';
    fitOptions.MaxFunEvals = 5000;
    
    % Objective function for data after echo 1, used for estimation of
    % correction parameters
    fun = @(beta,X) X.*exp(2i.*pi.*(beta - repmat(tF, [size(d,1)/2 1]).*(demod_freq)));
    beta0 = 1i.*log(2)/2/pi + 1-2*rand(1);
    beta_fit = lsqcurvefit(fun,beta0,double(gas_views), double(dissolved_views),[],[],fitOptions);

    % Objective funcion for all data, used for applying the correction
    fun = @(beta,X) X.*exp(2i.*pi.*(beta - repmat(tC, [size(d,1)/2 1]).*(demod_freq)));
    d(dpv,:) = d(dpv,:) - fun(beta_fit,d(gpv,:));
    
end

cart_down_fac = []; 

k_gp = LoadData.ismrmrd.GE.Functions.raw2grid(d(1:2:end,:,:,:,:,:),LoadData.ismrmrd.GE.Functions.ischop(h),k,shft,cart_down_fac,[],ind(1:2:end,:),delay*1d-6,bw,false);
k_dp = LoadData.ismrmrd.GE.Functions.raw2grid(d(2:2:end,:,:,:,:,:),LoadData.ismrmrd.GE.Functions.ischop(h),k,shft,cart_down_fac,[],ind(2:2:end,:),delay*1d-6,bw,false);

clear d wf ind
% size of reshaped data dd:
%  dim1=#coils; dim2=#indexed kspace data; dim3=#excitations; dim4=#slices
if any(lb>0), k_gp = ak_apodise(k_gp,k,lb(2),t,lb(1),false); k_dp = ak_apodise(k_dp,k,lb(2),t,lb(1),false); end

% Frequency demodulation
if f0~=0
    phafu = exp(1i*2*pi*t*f0);
    k_gp = k_gp.*phafu;
    k_dp = k_dp.*phafu;
end
clear t


%% Save data to mat file for offline processing
nexc = nexc/2;
nread = numel(k_dp)/nexc;
k_gp = reshape(k_gp,nread,nexc);
k_dp = reshape(k_dp,nread,nexc);

x_gp = reshape(k(:,:,1),nread,nexc);
y_gp = reshape(k(:,:,2),nread,nexc);
z_gp = reshape(k(:,:,3),nread,nexc);

bv = (abs(x_gp) > 0.5) | (abs(y_gp) > 0.5) | (abs(z_gp) > 0.5);

x_gp(any(bv,2),:) = [];
y_gp(any(bv,2),:) = [];
z_gp(any(bv,2),:) = [];
k_dp(any(bv,2),:) = [];
k_gp(any(bv,2),:) = [];

x_dp = x_gp;
y_dp = y_gp;
z_dp = z_gp;


date = ['20' h.rdb_hdr.scan_date(end-1:end) '-' h.rdb_hdr.scan_date(1:2) '-' h.rdb_hdr.scan_date(4:5)];
te_dp = h.rdb_hdr.te;
fov = fov*1e3;
bw = h.rdb_hdr.user0;

% if ~isempty(fname)
    [fp,fn] = fileparts(fname); cd(fp);
    save([fp '/Dixon_' fn '.mat'],'-mat','-v7.3',...
        'k_gp','k_dp','x_gp','y_gp','z_gp','x_dp','y_dp','z_dp','date','te_dp','fov', 'bw');
% end
