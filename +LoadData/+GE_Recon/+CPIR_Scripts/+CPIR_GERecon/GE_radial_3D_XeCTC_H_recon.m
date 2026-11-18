function [Proton] = GE_radial_3D_XeCTC_H_recon(MainInput, GasExchange, Proton)
%GE_RADIAL_3D_XECTC_H_RECON GE data adapter for radial_3D_XeCTC_H_recon.
%
% After GE loading, the rest of this function is copied verbatim from
% radial_3D_XeCTC_H_recon.m (including the original ProtonMontage figure).
% GE data are read using read_p + waveform + raw2grid, NOT recon_grid_interp.
% The original MRD motion, gridding, SOS, bias correction and plotting run.
%
% MainInput.HFullPath             : full path to GE proton ScanArchive .h5
% MainInput.HWaveformFile         : optional waveform .mat override
% MainInput.HDataLocation         : optional, derived from HFullPath
% MainInput.HFileName             : optional, derived from HFullPath
% MainInput.HReconMatrix          : optional matrix override
% MainInput.HScanVersion          : 'XeCTC' or 'CCHMC_V3' override
% MainInput.HGEUseNativeShift     : default true (native GE reader shift)
% MainInput.HGEApplyApodization  : default false (MRD pathway unfiltered)
% MainInput.OutputPath            : optional output folder if GasExchange
%                                   does not already define outputpath
%
% IMPORTANT: This is an unvalidated research conversion. Compare image
% scaling, spatial alignment and motion rejection with native GE images.

if nargin < 1 || ~isstruct(MainInput)
    error('MainInput must be a struct.');
end
if nargin < 2 || isempty(GasExchange), GasExchange = struct(); end
if nargin < 3 || isempty(Proton), Proton = struct(); end

%% GE-ONLY SECTION: Find GE data, read k-space and trajectory
HFullPath = getOption(MainInput,'HFullPath','');
if isempty(HFullPath)
    HFullPath = fullfile(getOption(MainInput,'HDataLocation',''), ...
                         getOption(MainInput,'HFileName',''));
end
assert(isfile(HFullPath),'Cannot find GE proton ScanArchive: %s',HFullPath);
[ProtonDataLocation,stem,ext] = fileparts(HFullPath);
if ~isfield(MainInput,'HFileName') || isempty(MainInput.HFileName)
    MainInput.HFileName = [stem ext];
end
MainInput.Scanner = 'GE';

% The earlier GE proton script records waveform and GE timing settings in
% a matlab_scripter file next to the ScanArchive. No native image is made.
geFiles = findGEProtonSettings(HFullPath,MainInput);
[HKSpace,HTraj,GE] = extractGERadialProton( ...
    HFullPath,geFiles.waveform,geFiles.mtx_reco, ...
    geFiles.interp_factor,geFiles.delay,geFiles.lb,MainInput);

% Match the original MRD variable names and array dimensions.
% HKSpace: [readout x projections x coil]
% HTraj:   [3 x readout x projections]
H_nsamp = size(HKSpace,1);
H_nprof = size(HKSpace,2);
H_coils = size(HKSpace,3);
H_interleaves = 1;
H_AcqMatrix = GE.mtx_acq(1);
H_RecMatrix = getOption(MainInput,'HReconMatrix',GE.mtx_reco(1));
assert(H_RecMatrix==round(H_RecMatrix) && H_RecMatrix>=16, ...
       'HReconMatrix must be an integer >=16.');

% GE k-space trajectory is normalized to GE.mtx_reco. If the user changes
% the MRD matrix, adjust the normalized trajectory consistently.
if H_RecMatrix ~= GE.mtx_reco(1)
    HTraj = HTraj * (GE.mtx_reco(1)/H_RecMatrix);
end
assert(isequal(size(HTraj),[3 H_nsamp H_nprof]), ...
    'GE k-space and trajectory dimensions disagree.');
assert(all(isfinite(real(HKSpace(:)))) && ...
       all(isfinite(imag(HKSpace(:)))) && all(isfinite(HTraj(:))), ...
       'Non-finite GE k-space/trajectory.');
assert(max(abs(HTraj(:)))<=0.5, ...
    'Normalized GE trajectory exceeds +/-0.5; verify waveform scaling.');

% Preserve MRD's scan-version-dependent choice of reconstruction solver.
Institute = getOption(MainInput,'Institute','CCHMC');
if strcmpi(Institute,'XeCTC')
    ScanVersion = 'XeCTC';
else
    ScanVersion = 'CCHMC_V3';
end
ScanVersion = getOption(MainInput,'HScanVersion',ScanVersion);
assert(ismember(ScanVersion,{'XeCTC','CCHMC_V3'}), ...
       'HScanVersion must be XeCTC or CCHMC_V3.');

% Retain the original output location preference: GasExchange.outputpath.
if isfield(GasExchange,'outputpath') && ~isempty(GasExchange.outputpath)
    outputpath = GasExchange.outputpath;
else
    outputpath = getOption(MainInput,'OutputPath', ...
                 fullfile(ProtonDataLocation,'GasExchange_Analysis'));
end
if ~isfolder(outputpath), mkdir(outputpath); end
Proton.outputfolder = outputpath;
Proton.DataFile = HFullPath;
Proton.HeaderFile = geFiles.header;
Proton.WaveformFile = geFiles.waveform;
Proton.MatlabScripterFile = geFiles.scripter;
Proton.GEPreprocessing = rmfield(GE,'h');
Proton.H_nsamp = H_nsamp;
Proton.H_nprof = H_nprof;
Proton.H_coils = H_coils;
Proton.H_interleaves = H_interleaves;

fprintf('GE proton: %d samples x %d views x %d coil(s), matrix %d.\n', ...
    H_nsamp,H_nprof,H_coils,H_RecMatrix);
fprintf('Using unchanged MRD proton reconstruction sections below.\n');

%% Account for Exhale/Motion in Proton Image
%Functions to fit
for coil = 1:size(HKSpace,3)
    HPhase(:,coil) = smooth(angle(HKSpace(1,:,coil)),15,'rloess');
end
HPhaseDynamics = rad2deg(mean(unwrap(HPhase,1),2));
HPulses = (1:length(HPhaseDynamics))';

%Fits to data; use only 1st third to avoid exhale/movement - assuming first 5 seconds are uncorrupted

HDynMean = mean(HPhaseDynamics(1000:round(length(HPulses)/3)));
HDynFit = HDynMean*ones(size(HPhaseDynamics(1000:round(length(HPulses)/3))));
HDynRMSE = sqrt(immse(HDynFit,HPhaseDynamics(1000:round(length(HPulses)/3))));
if isnan(HDynMean)
    HDynMean = mean(HPhaseDynamics(1:round(length(HPulses)/3)));
    HDynFit = HDynMean*ones(size(HPhaseDynamics(1:round(length(HPulses)/3))));
    HDynRMSE = sqrt(immse(HDynFit,HPhaseDynamics(1:round(length(HPulses)/3))));
end

if (3*HDynRMSE<5) 
    HDynRMSE = 5/3; %If small phase varaition, use acceptance window of 5deg
end

%smooth to remove small excursions
HPhaseDynamics_Smoothed = medfilt1(HPhaseDynamics,15,'truncate');

%Determine Proton Data Cutoffs
HDynFit = HDynMean*ones(size(HPulses));%calculate function at all x points
keep = abs(HDynFit - HPhaseDynamics_Smoothed) < 3*HDynRMSE;%label as outlier if >3 rmse away
keep = medfilt1(single(keep),15,'truncate');%smooth to remove small excursions
keep(1:round(length(HPulses)/3)) = 1; %keep all from fit
motion_ind = find(keep==0,1,'first');%find first non keep index

%H Motion
HKSpace_Steady = HKSpace;
HTraj_Steady = HTraj;
if ~isempty(motion_ind)%only remove data if motion
    HKSpace_Steady(:,motion_ind:end,:) = [];
    HTraj_Steady(:,:,motion_ind:end) = [];
else
    motion_ind = length(HPulses);
end

%% Determine FOV Shift
% %not using planned shifts as techs don't seem to do it consistently
% disp('Determining Image offset...')
% PixelShiftinit = [0; 0; 0]; %start with no offset
% converged = 0; %start not converged
% iter = 0;
% while (converged == 0) %run until no change in shift
%     
%     tempIm = GasExchange.Dissolved_HImageRecon(H_RecMatrix,HKSpace_Steady(:,:,1),HTraj_Steady/2,PixelShiftinit); %make image; 2x Resolution
% %     tempIm = GasExchange.Dissolved_LowResImageRecon(Xe_RecMatrix,GasKSpace,XeTraj/2,PixelShiftinit); %make image; 2x Resolution
%     [tempIm, ~] = GasExchange.Vent_BiasCorrection_Coarse(tempIm);
%     tempLab = bwlabeln(imbinarize(abs(tempIm)/(0.5*max(abs(tempIm(:)))))); %make mask label
%     center = regionprops3(tempLab); %get centroids
%     center = sortrows(center,1,'descend'); %order by size
%     if (size(center,1) > 1 && table2array(center(2,1)) > 0.15 * table2array(center(1,1)))        
%         CentVox(1) = (min(center.BoundingBox(:,1))+max(center.BoundingBox(:,1)+center.BoundingBox(:,4)))/2;
%         CentVox(2) = (min(center.BoundingBox(:,2))+max(center.BoundingBox(:,2)+center.BoundingBox(:,5)))/2;
%         CentVox(3) = (min(center.BoundingBox(:,3))+max(center.BoundingBox(:,3)+center.BoundingBox(:,6)))/2;
%     else
%         CentVox(1) = center.BoundingBox(1,1)+center.BoundingBox(1,4)/2;
%         CentVox(2) = center.BoundingBox(1,2)+center.BoundingBox(1,5)/2;
%         CentVox(3) = center.BoundingBox(1,3)+center.BoundingBox(1,6)/2;
%     end
%     PixelShift = 0.5*H_RecMatrix - CentVox'; %get pixel shift
%     PixelShift = PixelShift([3, 2, 1]); %put in same order
%     PixelShift(2:3) = -PixelShift(2:3); %correct for reference
%     PixelShift = PixelShift + PixelShiftinit; %sum previous shifts
%     iter = iter + 1;
%     if (sum(abs(PixelShift - PixelShiftinit) > 1) == 0) %if all offsets changed less than 1 pixel, end
%         converged = 1;
%     elseif (iter == 5) %if reached max iterations
%         converged = 1;
%         if (sum(abs(PixelShift - PixelShiftinit) < 3) ~= 3) %if all offsets changed less than 3 pixels, keep, else set to 0
%             PixelShift = [0; 0; 0]; %move back to no offset
%         end
%         disp('Did not converge on image offset...')
%     else %continue loop
%         PixelShiftinit = PixelShift; %set init for next loop
%     end
% end
% disp('Determining Image Offset Completed.')

% PixelShift = [8.5; 1.5; 2]; %start with no offset
%% Proton Image
%Reconstruct Image
try
    PixelShift = GasExchange.PixelShift; % [1.5;4;-4];
catch
    PixelShift = [0,0,0];
end

disp('Reconstructing UTE Image...')
overgridding_factor = 3;
ProtonImage = zeros(H_RecMatrix*overgridding_factor,H_RecMatrix*overgridding_factor,H_RecMatrix*overgridding_factor,H_coils);%3x overgridding
for i = 1:H_coils
    if strcmp(ScanVersion,'XeCTC') %recon 2x interpolated
        ProtonImage(:,:,:,i) = GasExchangeFunctions.Dissolved_HImageRecon_CTC(H_RecMatrix,HKSpace_Steady(:,:,i),HTraj_Steady/2,PixelShift);
    else
        ProtonImage(:,:,:,i) = GasExchangeFunctions.Dissolved_HImageRecon(H_RecMatrix,HKSpace_Steady(:,:,i),HTraj_Steady,PixelShift);
    end
end
ProtonImage = GasExchangeFunctions.SOS_Coil_Combine(ProtonImage);%Combine Coils
disp('Reconstructing UTE Image Competed.')
if strcmp(MainInput.Scanner,'Siemens')
    ProtonImage = permute(ProtonImage, [2,1,3]);
    ProtonImage = flip(ProtonImage,1);
    ProtonImage = flip(ProtonImage,2);
end
ProtonImage = flip(imrotate(ProtonImage,-90),1);

ProtonImageHR = ProtonImage;
figure; imslice(abs(ProtonImageHR))

%Correct Bias
disp('Correcting Proton Bias...')
ProtonImageLR = medfilt3(ProtonImage, [7 7 7]);%remove undersampling artifacts
ProtonImageLR = imgaussfilt3(ProtonImageLR, 0.5);%reduce noise
[ProtonImageLR, ~] = GasExchangeFunctions.Dissolved_ProtonBiasCorrection(ProtonImageLR);
[ProtonImageHR, ~] = GasExchangeFunctions.Dissolved_ProtonBiasCorrection(ProtonImageHR);

% Determine Levels
ProtonMax = prctile(abs(ProtonImageLR(:)),99.99);

disp('Proton Bias Corrected.')  
cd(outputpath)
%View
ProtonMontage = figure('Name','Proton Image');set(ProtonMontage,'WindowState','minimized');
montage(abs(ProtonImage(H_RecMatrix:overgridding_factor*H_RecMatrix,H_RecMatrix:overgridding_factor*H_RecMatrix,H_RecMatrix:overgridding_factor*H_RecMatrix)),'DisplayRange',[0 ProtonMax])%unregistered for these
savefig('ProtonMontage.fig')
close(gcf)

Proton.Image = double(ProtonImageHR);
Proton.ProtonImageLR = double(ProtonImageLR);
Proton.ProtonImageHR = double(ProtonImageHR);
Proton.filename = MainInput.HFileName;
Proton.folder = ProtonDataLocation;
Proton.H_RecMatrix = H_RecMatrix;
Proton.ProtonMax = ProtonMax;

save(fullfile(outputpath, 'proton_recon_workspace.mat'));

end

%% GE-only helper: get waveform and settings (original GE matlab_scripter)
function p = findGEProtonSettings(scan,MainInput)
[folder,base,ext] = fileparts(scan);
p.header = fullfile(folder,[base '.hdr']);
if ~isfile(p.header)
    found = dir(fullfile(folder,'*.hdr'));
    assert(numel(found)==1,'Missing/ambiguous GE .hdr next to ScanArchive.');
    p.header=fullfile(folder,found(1).name);
end
p.scripter = '';
p.waveform = getOption(MainInput,'HWaveformFile','');
p.mtx_reco = getOption(MainInput,'HGEReconMatrix',[]);
p.interp_factor = getOption(MainInput,'HGEInterpFactor',1);
p.delay = getOption(MainInput,'HGEDelay_us',0);
p.lb = getOption(MainInput,'HGELB',0);

% An explicit waveform path skips the scripter and permits manual settings.
if ~isempty(p.waveform)
    assert(isfile(p.waveform),'HWaveformFile does not exist: %s',p.waveform);
    return;
end

files = dir(fullfile(folder,'*matlab_scripter*'));
assert(~isempty(files),'No GE matlab_scripter found in %s.',folder);
idx=[];
for i=1:numel(files)
    candidate=fullfile(folder,files(i).name);
    try
        if contains(fileread(candidate),[base ext],'IgnoreCase',true)
            idx=i;
            break;
        end
    catch
    end
end
if isempty(idx)
    assert(numel(files)==1, ...
        'Multiple matlab_scripter files exist but none matches GE scan.');
    idx=1;
end
p.scripter=fullfile(folder,files(idx).name);
txt=fileread(p.scripter);
% GE proton anatomy uses recon_grid_interp, unlike GE xenon gas exchange.
token=regexp(txt,'recon_grid_interp\s*\((.*?)\)\s*;', 'tokens','once');
assert(~isempty(token),'No recon_grid_interp(...) command in GE scripter.');
args=splitGEArguments(token{1});
assert(numel(args)>=12,'Expected 12 recon_grid_interp arguments.');
path=strrep(strrep(strtrim(args{3}),'''',''),'\','/');
[~,name,extension]=fileparts(path);
if isempty(extension) || strcmpi(extension,'.wav'), extension='.mat'; end
p.waveform=fullfile(folder,[name extension]);
assert(isfile(p.waveform),'GE proton waveform not found: %s',p.waveform);
p.mtx_reco=parseGEValue(args{4});
p.interp_factor=parseGEValue(args{5});
p.delay=parseGEValue(args{6});
p.lb=parseGEValue(args{7});
end

function args = splitGEArguments(s)
args={}; start=1; depth=0; quoted=false;
for i=1:numel(s)
    c=s(i);
    if c=='''', quoted=~quoted; end
    if ~quoted
        if any(c=='[({'),depth=depth+1;end
        if any(c=='])}'),depth=depth-1;end
        if c==',' && depth==0
            args{end+1}=strtrim(s(start:i-1)); %#ok<AGROW>
            start=i+1;
        end
    end
end
args{end+1}=strtrim(s(start:end));
end

function x = parseGEValue(s)
s=strtrim(s);
if strcmp(s,'[]'),x=[];return;end
if strcmpi(s,'true'),x=true;return;end
if strcmpi(s,'false'),x=false;return;end
if startsWith(s,'[') && endsWith(s,']'),s=s(2:end-1);end
s=strrep(s,',',' ');
parts=strsplit(s,';'); rows=cell(size(parts));
for i=1:numel(parts)
    rows{i}=sscanf(parts{i},'%f').';
    assert(~isempty(rows{i}), 'Unsupported GE numeric setting: %s',s);
end
x=vertcat(rows{:});
end

function value = getOption(s,name,default)
value=default;
if isstruct(s) && isfield(s,name) && ~isempty(s.(name))
    value=s.(name);
end
end

function [HKSpace,HTraj,info] = extractGERadialProton(scanFile,wfn,mtx_reco,interp_factor,delay,lb,MainInput)
% GE read_p + native waveform + raw2grid ONLY. No native GE gridding.
RFS_MAX_NSCANS = 16382;
[d,h] = LoadData.GE_Recon.fidall.read_p(scanFile,true);
d = single(d); % Native recon_grid_interp does NOT conjugate proton data.

if isempty(interp_factor),interp_factor=1;end
if isempty(mtx_reco),mtx_reco=0;end
if isscalar(mtx_reco),mtx_reco=repmat(mtx_reco,1,3);end
if numel(mtx_reco)==2,mtx_reco=[mtx_reco mtx_reco(end)];end
mtx_reco=double(mtx_reco(:).');
if isempty(delay),delay=0;end
if numel(delay)>1,delay=delay(1);end
if isempty(lb),lb=0;end
if isscalar(lb),lb=repmat(lb,1,3);end
if numel(lb)==2,lb=[lb lb(end)];end

if iscell(wfn)
    wfn=wfn{h.image.user3};
end
if isstruct(wfn)
    wf=wfn;
else
    if contains(lower(char(h.image.psdname)),'spiral')
        error('Only 3D radial proton data supported (GE PSD is spiral).');
    end
    if ~endsWith(lower(wfn),'.mat'),wfn=[wfn '.mat'];end
    wf=LoadData.GE_Recon.fidall.load_waveform(wfn,true,true);
end
for fn={'mtx','k','dcf','ind','t'}
    if ~isfield(wf,fn{1}),error('Missing waveform field %s.',fn{1});end
end
mtx_acq=double(wf.mtx(:).');
if isscalar(mtx_acq),mtx_acq=repmat(mtx_acq,1,3);end
if numel(mtx_acq)==2,mtx_acq=[mtx_acq mtx_acq(end)];end
if any(mtx_reco(1:3)<1)
    mtx_reco=interp_factor*mtx_acq(1:3);
    k=wf.k/interp_factor;
else
    k=bsxfun(@times,wf.k,shiftdim(mtx_acq(1:3)./mtx_reco(1:3),-1));
end
ind=logical(wf.ind);
if ndims(k)~=3 || size(k,3)~=3
    error('Expected 3D GE radial waveform k [1 x N x 3], received [%s].',num2str(size(k)));
end
if ~all(isfinite(k(:)))
    error('GE waveform contains NaN/Inf coordinates.');
end
nInterleaves=size(ind,1);
readoutCounts=sum(ind,2);
if any(readoutCounts ~= readoutCounts(1))
    error('Nonuniform GE waveform readout sample counts are unsupported.');
end
nRead=double(readoutCounts(1));
if nRead<8,error('Too few GE radial readout samples (%d).',nRead);end
if size(k,2)~=nRead*nInterleaves
    error(['Waveform trajectory length %d does not equal %d x %d. ' ...
        'Review GE ind / radial projection organization.'], ...
        size(k,2),nRead,nInterleaves);
end
if sum(ind(:)) > RFS_MAX_NSCANS && size(d,4)>1
    [n1,n2,n3,n4,n5,n6]=size(d);
    tmp=complex(zeros(n1*n4,n2,n3,1,n5,n6,'like',d));
    for ix=1:n4
        sel=(1:n1)+(ix-1)*n1;
        tmp(sel,:,:,:,:,:)=d(:,:,:,ix,:,:);
    end
    d=tmp;
end
[ne,nread,nphases,nEcho,nslices,ncoils]=size(d);
if nslices~=1 || nphases~=1
    error('Only single 3D proton slab supported. GE data dims [%s].',num2str(size(d)));
end
if ne < nInterleaves
    error('GE acquisitions (%d) fewer than waveform spokes (%d).',ne,nInterleaves);
end

bw=h.rdb_hdr.user0;
if abs(bw)<500
    bw=h.rdb_hdr.bw*2e3;
end
% Native GE phase cycling, fractional sampling delay, header-derived shift.
shft = zeros(1,3);
fov=double(h.rdb_hdr.fov)*1e-3;
niso=[1 1 1];
if isfield(wf,'fov') && numel(wf.fov)==3 && ...
        (any(diff(wf.fov)) || any(diff(wf.mtx)))
    niso=double(wf.fov(1))./double(wf.fov(:).');
end
if abs(h.rdb_hdr.user26)>0.01 || abs(h.rdb_hdr.user27)>0.01 || abs(h.rdb_hdr.user28)>0.01
    shft=mtx_reco.*double([h.rdb_hdr.user26 h.rdb_hdr.user27 h.rdb_hdr.user28])*1e-3/fov.*niso;
end
% The GE reader may perform the native half-voxel shift before MRD gridding.
if getOption(MainInput,'HGEUseNativeShift',true)
    shft=shft + 0.5; % native GE proton half-voxel shift
else
    shft=[0 0 0];   % leave image-centering to MRD PixelShift
end

k=double(k);
dd=LoadData.GE_Recon.fidall.raw2grid(d, ...
    LoadData.GE_Recon.fidall.ischop(h),k,shft,[],[],ind,delay*1e-6,bw,false);
clear d

% Native optional apodization (retained because user asked to preserve GE
% preprocessing); MRD kernel/DCF are subsequently used for reconstruction.
if getOption(MainInput,'HGEApplyApodization',false) && any(abs(lb)>0)
    t=wf.t;
    if size(t,1)==1 && size(t,2)<=size(ind,2)
        t=repmat(t,[size(ind,1),1]);
    end
    t=t.';t=t(:).';
    if numel(t)~=size(k,2)
        error('GE waveform time length inconsistent with trajectory.');
    end
    dd=LoadData.GE_Recon.fidall.ak_apodise( ...
        dd,k,lb(2),t,[lb(1) lb(3)],false);
end

% GE native gridding treats dd dimension 2 as a full volume's indexed
% k-space samples and dd dimension 3 as temporal frames / echoes. The MRD
% approach uses a single [readout x projection x coil] dataset.
if size(dd,2)~=size(k,2) || size(dd,3)~=1 || size(dd,4)~=1
    error(['Unexpected GE raw2grid dimensions [%s]. ' ...
        'Expected [%d coils, %d k-points, 1 frame, 1 slab]. ' ...
        'Do not silently mix time frames / echoes.'], ...
        num2str(size(dd)),ncoils,size(k,2));
end
HKSpace=complex(zeros(nRead,nInterleaves,ncoils,'single'));
for ic=1:ncoils
    HKSpace(:,:,ic)=reshape(dd(ic,:,1,1),nRead,nInterleaves);
end
HTraj=permute(reshape(k,1,nRead,nInterleaves,3),[4,2,3,1]);
HTraj=reshape(HTraj,3,nRead,nInterleaves);

info=struct();
info.h=h;
info.mtx_reco=mtx_reco;
info.mtx_acq=mtx_acq;
info.niso=niso;
info.shft=shft;
info.fov_m=fov;
info.nread=nRead;
info.nproj=nInterleaves;
info.ncoils=ncoils;
info.ninterleaves=1;
info.nviewsPerWaveform=nInterleaves;
info.bw=bw;
info.nativeDCFAvailable=isfield(wf,'dcf') && ~isempty(wf.dcf);
info.nativeApodizationApplied=logical(getOption(MainInput,'HGEApplyApodization',false));
info.nativeShiftApplied=logical(getOption(MainInput,'HGEUseNativeShift',true));

fprintf('GE proton extracted: [%d samples x %d projections x %d coils].\n', ...
    nRead,nInterleaves,ncoils);
end
