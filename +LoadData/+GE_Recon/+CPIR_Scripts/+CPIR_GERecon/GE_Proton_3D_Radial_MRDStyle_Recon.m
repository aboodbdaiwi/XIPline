function Proton = GE_Proton_3D_Radial_MRDStyle_Recon(MainInput)

% GE_Proton_3D_Radial_MRDStyle_Recon
% GE ScanArchive/waveform/raw2grid preprocessing, with MRD-style proton
% gridding (GasExchangeFunctions.Dissolved_HImageRecon[_CTC]).
%
% Original GE anatomy output contract is retained. The new reconstruction
% uses normalized GE waveforms, iterative DCF, per-coil reconstruction,
% SOS coil combination, optional k0-phase motion rejection, central FOV
% cropping, and the original GE orientation/bias-correction pathways.
%
% NO native GE gridding (recon_grid_interp) is called.
%
% IMPORTANT: Experimental research pipeline, not clinically validated.
% Validate orientation, spatial scaling, image sharpness and residual motion
% against the native GE reconstruction for representative acquisitions.
%
% Optional MainInput fields:
%   HReconMatrix      Positive scalar, overrides GE reconstruction matrix.
%   HMRDMethod        ''CTC'' or ''standard''. Default CTC for XeCTC,
%                     standard otherwise (as in the MRD function).
%   HMotionFilter     Logical, default true. Conservative phase-based cutoff.
%   HExtraTrajScale   Positive scalar, default 1. NEVER multiply by matrix.
%   HPixelShift       [3x1] MRD shift in voxels, default [0;0;0].
%   HApplyGradwarp    Logical, default from native GE scripter ''grdwrp''.
%   HUseGEOrientation Logical, default true.
%   HKeepKspace       Logical, default false (avoid large output struct).
%   HKeepCoils        Logical, default false (avoid large 4D data).
%   HPreReconOnly     Logical, default false (pre-recon data/QC mode).
%   HSaveQC           Logical, default false (save QC PNG and FIG).
%
% Requires XIPline GE fidall and GasExchangeFunctions packages + MEX. ''CTC''
% outputs a 3x oversampled reconstruction internally; this wrapper crops
% the central requested N^3 FOV before bias correction and XIPline output.
%
% Input:
%   MainInput.HFullPath      - Full path to GE proton ScanArchive .h5 file
%   MainInput.HDataLocation  - Data folder (optional)
%   MainInput.HFileName      - Data filename (optional)
%
% Output:
%   Proton - structure containing reconstructed proton images,
%                 reconstruction settings, source files, and metadata.

if nargin < 1 || isempty(MainInput)
    error('MainInput is required.');
end

if ~isfield(MainInput,'HFullPath') || isempty(MainInput.HFullPath)
    error('MainInput.HFullPath is missing or empty.');
end

d = MainInput.HFullPath;

if ~isfile(d)
    error('GE proton data file was not found:\n%s',d);
end

[dataFolder,dataName,dataExt] = fileparts(d);

fprintf('\n============================================================\n');
fprintf('GE Proton 3D Radial MRD-style Reconstruction\n');
fprintf('============================================================\n');
fprintf('Data file:\n%s\n',d);

%% Find corresponding GE header file

expectedHeader = fullfile(dataFolder,[dataName '.hdr']);

if isfile(expectedHeader)

    h = expectedHeader;

else

    headerFiles = dir(fullfile(dataFolder,'*.hdr'));

    if isempty(headerFiles)
        error('No GE .hdr file was found in:\n%s',dataFolder);
    end

    matchIndex = [];

    for ii = 1:numel(headerFiles)

        [~,headerName] = fileparts(headerFiles(ii).name);

        if strcmpi(headerName,dataName)
            matchIndex = ii;
            break;
        end

    end

    if isempty(matchIndex)

        if numel(headerFiles) == 1
            matchIndex = 1;
        else
            error(['Could not uniquely identify the GE header file.\n' ...
                   'Expected:\n%s'],expectedHeader);
        end

    end

    h = fullfile( ...
        headerFiles(matchIndex).folder, ...
        headerFiles(matchIndex).name);

end

fprintf('Header file:\n%s\n',h);

%% Find matlab_scripter file

scripterFiles = dir(fullfile(dataFolder,'*matlab_scripter*'));

if isempty(scripterFiles)
    error(['Could not find a file containing "matlab_scripter" in:\n%s'], ...
        dataFolder);
end

scripterIndex = [];

for ii = 1:numel(scripterFiles)

    candidatePath = fullfile( ...
        scripterFiles(ii).folder, ...
        scripterFiles(ii).name);

    try

        candidateText = fileread(candidatePath);

        if contains(candidateText,[dataName dataExt], ...
                'IgnoreCase',true)

            scripterIndex = ii;
            break;

        end

    catch
    end

end

if isempty(scripterIndex)

    if numel(scripterFiles) == 1
        scripterIndex = 1;
    else
        error(['Multiple matlab_scripter files were found, but none could ' ...
               'be uniquely matched to:\n%s'], ...
               [dataName dataExt]);
    end

end

scripterPath = fullfile( ...
    scripterFiles(scripterIndex).folder, ...
    scripterFiles(scripterIndex).name);

fprintf('MATLAB scripter:\n%s\n',scripterPath);

scripterText = fileread(scripterPath);

%% Extract native GE reconstruction settings (NOT native gridding)

commandToken = regexp( ...
    scripterText, ...
    'recon_grid_interp\s*\((.*?)\)\s*;', ...
    'tokens', ...
    'once');

if isempty(commandToken)
    error(['Could not find a valid recon_grid_interp(...) command in:\n' ...
           '%s'],scripterPath);
end

commandArguments = strtrim(commandToken{1});

%% Split MATLAB function arguments safely

args = {};
currentArg = '';
insideQuote = false;
bracketDepth = 0;

for ii = 1:length(commandArguments)

    c = commandArguments(ii);

    if c == ''''

        insideQuote = ~insideQuote;
        currentArg(end+1) = c;

    elseif ~insideQuote && (c == '[' || c == '(' || c == '{')

        bracketDepth = bracketDepth + 1;
        currentArg(end+1) = c;

    elseif ~insideQuote && (c == ']' || c == ')' || c == '}')

        bracketDepth = bracketDepth - 1;
        currentArg(end+1) = c;

    elseif c == ',' && ~insideQuote && bracketDepth == 0

        args{end+1} = strtrim(currentArg);
        currentArg = '';

    else

        currentArg(end+1) = c;

    end

end

if ~isempty(strtrim(currentArg))
    args{end+1} = strtrim(currentArg);
end

if numel(args) < 12
    error(['The recon_grid_interp command does not contain the expected ' ...
           '12 arguments.\n\nParsed command:\n%s'], ...
           commandArguments);
end

% Expected recon_grid_interp arguments:
%
% 1  = d
% 2  = h
% 3  = wfn
% 4  = mtx_reco
% 5  = interp_factor
% 6  = delay
% 7  = lb
% 8  = fname
% 9  = mask_fac
% 10 = coco
% 11 = plt
% 12 = grdwrp

%% Extract waveform filename

waveformArg = strtrim(args{3});

if startsWith(waveformArg,'''') && endsWith(waveformArg,'''')
    waveformOriginalPath = waveformArg(2:end-1);
else
    waveformOriginalPath = waveformArg;
end

waveformNormalized = strrep(waveformOriginalPath,'\','/');

slashIndex = find(waveformNormalized == '/',1,'last');

if isempty(slashIndex)
    waveformFileName = waveformNormalized;
else
    waveformFileName = waveformNormalized(slashIndex+1:end);
end

wfn = fullfile(dataFolder,waveformFileName);

if ~isfile(wfn)

    waveformFiles = dir(fullfile(dataFolder,'*.mat'));

    waveformIndex = [];

    for ii = 1:numel(waveformFiles)

        if strcmpi(waveformFiles(ii).name,waveformFileName)
            waveformIndex = ii;
            break;
        end

    end

    if isempty(waveformIndex)
        error(['Waveform file specified by matlab_scripter:\n%s\n\n' ...
               'was not found in:\n%s'], ...
               waveformFileName,dataFolder);
    end

    wfn = fullfile( ...
        waveformFiles(waveformIndex).folder, ...
        waveformFiles(waveformIndex).name);

end

fprintf('Waveform file:\n%s\n',wfn);

%% Extract reconstruction parameters

mtx_reco = parseMatlabValue(args{4});
interp_factor = parseMatlabValue(args{5});
delay = parseMatlabValue(args{6});
lb = parseMatlabValue(args{7});
mask_fac = parseMatlabValue(args{9});
coco = parseMatlabValue(args{10});
plt = parseMatlabValue(args{11});
grdwrp = parseMatlabValue(args{12});

% Do not use the original Linux fname path.
% XIPline will manage the returned reconstruction. Native GE gridding is not used.
fname = [];

fprintf('\nReconstruction settings:\n');
fprintf('  mtx_reco      = %s\n',mat2str(mtx_reco));
fprintf('  interp_factor = %s\n',mat2str(interp_factor));
fprintf('  delay         = %s\n',mat2str(delay));
fprintf('  lb            = %s\n',mat2str(lb));
fprintf('  mask_fac      = %s\n',mat2str(mask_fac));
fprintf('  coco          = %s\n',mat2str(coco));
fprintf('  plt           = %s\n',mat2str(plt));
fprintf('  grdwrp        = %s\n',mat2str(grdwrp));


%% GE raw preparation (NO native GE gridding)

fprintf('\nReading GE proton k-space and waveform trajectories...\n');
[HKSpace, HTraj, geInfo] = extractGERadialProton( ...
    d, wfn, mtx_reco, interp_factor, delay, lb);

[nReadout,nProjections,nCoils] = size(HKSpace);
maxAxis = max(abs(HTraj(:)));
maxRadius = max(sqrt(sum(HTraj.^2,1)),[], 'all');

fprintf('GE proton data: %d readout samples x %d projections x %d coils\n', ...
    nReadout,nProjections,nCoils);
fprintf('Normalized |k| max radius: %.6f, max axis: %.6f\n', ...
    maxRadius,maxAxis);

% Recompute trajectory safely only if an explicit scaling factor is supplied.
trajScale = getOption(MainInput,'HExtraTrajScale',1);
assert(isscalar(trajScale) && isfinite(trajScale) && trajScale>0, ...
    'MainInput.HExtraTrajScale must be a finite positive scalar.');
HTraj = HTraj * trajScale;

% The MRD CTC reconstruction uses a 3x oversampled Cartesian grid internally.
% The legacy sparse_gridding_distance_mex has unsigned bounds and can HANG
% MATLAB for out-of-grid coordinates. Explicitly prevent unsafe calls.
H_RecMatrix = getOption(MainInput,'HReconMatrix',geInfo.mtx_reco(1));
assert(isnumeric(H_RecMatrix) && isscalar(H_RecMatrix) && ...
    isfinite(H_RecMatrix) && H_RecMatrix>=16 && ...
    H_RecMatrix==round(H_RecMatrix), ...
    'HReconMatrix must be a positive integer >= 16.');
assert(all(geInfo.mtx_reco(1:3)==geInfo.mtx_reco(1)), ...
    'The MRD proton reconstruction requires an isotropic native GE matrix.');

% GE waveform k was normalized for the native GE mtx_reco. If a different
% MRD reconstruction matrix is requested, change normalized coordinates
% consistently (same acquired physical k-space; different voxel size).
if H_RecMatrix ~= geInfo.mtx_reco(1)
    fprintf('Rescaling normalized GE trajectory for matrix %d -> %d.\n', ...
        geInfo.mtx_reco(1),H_RecMatrix);
    HTraj = HTraj .* (geInfo.mtx_reco(1) / H_RecMatrix);
end

maxAxis = max(abs(HTraj(:)));
maxRadius = max(sqrt(sum(HTraj.^2,1)),[], 'all');
safetyMargin = 3 / (3 * H_RecMatrix); % avoid unsigned out-of-range MEX bounds
if maxAxis > (0.5 - safetyMargin)
    error(['UNSAFE GE trajectory: max absolute normalized axis = %.5f. ' ...
           'MRD gridding MEX needs fractional k-space coordinates within +/-0.5. ' ...
           'Do not scale k by the reconstruction matrix.'], maxAxis);
end
if maxRadius < 0.01
    warning('Very small k-space radius (%.5f). Verify GE trajectory scaling.',maxRadius);
end
if ~all(isfinite(real(HKSpace(:)))) || ...
   ~all(isfinite(imag(HKSpace(:)))) || ~all(isfinite(HTraj(:)))
    error('GE raw k-space or trajectory contains NaN/Inf.');
end

% A conservative version of the MRD k0 phase change-point detector.
% Uses the actual per-projection trajectory-center sample; does not assume
% k-space center must be sample 1. No profiles are removed without evidence
% for a persistent late-phase shift.
motionEnabled = logical(getOption(MainInput,'HMotionFilter',true));
[motionKeep,phaseDynamics,motionInfo] = ...
    findProtonMotionCutoff(HKSpace,HTraj,motionEnabled);
HKSpace_Steady = HKSpace(:,motionKeep,:);
HTraj_Steady = HTraj(:,:,motionKeep);
if nnz(motionKeep) < 0.25*nProjections
    error('Motion selection kept fewer than 25%% of proton projections.');
end
fprintf('Motion selection: retained %d/%d projections (reason: %s)\n', ...
    nnz(motionKeep),nProjections,motionInfo.reason);

% Intermediate diagnostic mode avoids entering the expensive compiled MEX.
preReconOnly = logical(getOption(MainInput,'HPreReconOnly',false));
if preReconOnly
    Proton = struct();
    Proton.PreReconOnly = true;
    Proton.HKSpace = HKSpace_Steady;
    Proton.HTraj = HTraj_Steady;
    Proton.GEInfo = geInfo;
    Proton.Motion = motionInfo;
    Proton.HPhaseDynamics = phaseDynamics;
    Proton.H_RecMatrix = H_RecMatrix;
    Proton.H_coils = nCoils;
    Proton.DataFile = d;
    fprintf('HPreReconOnly=true: returning raw k-space; gridding not run.\n');
    return;
end

%% Reconstruct proton using custom MRD gridding

if isfield(MainInput,'HMRDMethod') && ~isempty(MainInput.HMRDMethod)
    reconMode = lower(char(MainInput.HMRDMethod));
elseif isfield(MainInput,'Institute') && strcmpi(MainInput.Institute,'XeCTC')
    reconMode = 'ctc';
else
    reconMode = 'standard';
end
if ~ismember(reconMode,{'ctc','standard'})
    error('HMRDMethod must be ''CTC'' or ''standard''.');
end

PixelShift = double(getOption(MainInput,'HPixelShift',[0;0;0]));
PixelShift = PixelShift(:);
assert(numel(PixelShift)==3 && all(isfinite(PixelShift)), ...
    'MainInput.HPixelShift must be a finite [3x1] vector.');
if any(PixelShift ~= 0)
    warning('GE raw2grid already applied native GE spatial shifts. Check that HPixelShift does not duplicate them.');
end

if H_RecMatrix > 128
    warning(['Requested reconstruction matrix N=%d produces an internal ' ...
             '(3N)^3 complex volume per coil. Consider N=96 or 128 ' ...
             'for an initial test to avoid excessive memory use.'],H_RecMatrix);
end
fprintf('\nMRD-style proton: mode=%s, requested matrix=%d, coils=%d\n', ...
    reconMode,H_RecMatrix,nCoils);
if strcmp(reconMode,'ctc')
    fprintf('Kernel: sigma=0.30, extent=0.60, overgrid=3, iterative DCF=15\n');
else
    fprintf('Kernel: sigma=0.175, extent=1.575, overgrid=3, iterative DCF=15\n');
end

% Each custom MRD function returns a 3N^3 *complex* volume (crop=false).
% To preserve the original GE Proton output dimensions, crop to center N^3
% immediately after each coil. This also limits peak memory use.
% Avoid the MRD script's explicit 384x384x384xNcoils allocation.
sosSum = zeros([H_RecMatrix H_RecMatrix H_RecMatrix],'single');
keepCoils = logical(getOption(MainInput,'HKeepCoils',false));
if keepCoils
    coilsComplex = complex(zeros(H_RecMatrix,H_RecMatrix,H_RecMatrix,nCoils,'single'));
else
    coilsComplex = [];
end
singleCoilComplex = [];

for ic = 1:nCoils
    fprintf('MRD-style 3D gridding: coil %d/%d\n',ic,nCoils);
    ticCoil = tic;
    if strcmp(reconMode,'ctc')
        imgOvergridded = GasExchangeFunctions.Dissolved_HImageRecon_CTC( ...
            H_RecMatrix, double(HKSpace_Steady(:,:,ic)), HTraj_Steady, PixelShift);
    else
        imgOvergridded = GasExchangeFunctions.Dissolved_HImageRecon( ...
            H_RecMatrix, double(HKSpace_Steady(:,:,ic)), HTraj_Steady, PixelShift);
    end
    imageDims = size(imgOvergridded);
    imageDims(end+1:3) = 1;
    if ~isequal(imageDims(1:3),3*H_RecMatrix*ones(1,3))
        error('Unexpected MRD output size [%s], expected [%d %d %d].', ...
            num2str(size(imgOvergridded)),3*H_RecMatrix*ones(1,3));
    end
    croppedCoil = single(centerCrop3D(imgOvergridded,H_RecMatrix));
    clear imgOvergridded
    sosSum = sosSum + abs(croppedCoil).^2;
    if nCoils==1
        singleCoilComplex = croppedCoil;
    end
    if keepCoils
        coilsComplex(:,:,:,ic) = croppedCoil;
    end
    clear croppedCoil
    fprintf('Coil %d complete (%.1f sec)\n',ic,toc(ticCoil));
end

% SOS combines magnitude without imposing unmeasured inter-coil phase.
protonSOS = sqrt(sosSum);
clear sosSum HKSpace_Steady

%% GE post-reconstruction orientation, distortion correction, bias correction

% The custom MRD helper applies rot90(flip(...)) internally. Undo that
% before reproducing GE's native post-gridding fliplr/mri_orientation and
% original GE_Proton wrapper permute/rotate/flip orientation.

useGEOrientation = 0; % logical(getOption(MainInput,'HUseGEOrientation',true));
if useGEOrientation
    protonSOS = undoMRDRotation(protonSOS);
    if ~isempty(singleCoilComplex)
        singleCoilComplex = undoMRDRotation(singleCoilComplex);
    end
    if keepCoils
        for ic=1:nCoils
            coilsComplex(:,:,:,ic) = undoMRDRotation(coilsComplex(:,:,:,ic));
        end
    end
end

% GE native gridding performs mask -> gradient warp -> orientation.
if ~isempty(mask_fac) && any(mask_fac>0)
    try
        M = LoadData.GE_Recon.fidall.calc_mask(protonSOS,2,0.1*mask_fac,[],[],true);
        protonSOS = protonSOS .* single(M);
    catch ME
        warning('Optional native GE mask was not applied: %s',ME.message);
    end
end

applyGradwarp = logical(getOption(MainInput,'HApplyGradwarp', ...
    ~isempty(grdwrp) && islogicalOrNumericTrue(grdwrp)));
gradwarpApplied = false;
if applyGradwarp
    try
        if ischar(grdwrp) || isstring(grdwrp)
            gradcoil = char(grdwrp);
        else
            gradcoil = '';
        end
        if isfield(geInfo,'niso') && numel(geInfo.niso)>=3
            zfov = double(geInfo.h.rdb_hdr.fov)/double(geInfo.niso(3));
        else
            zfov = double(geInfo.h.rdb_hdr.fov);
        end
        protonSOS = LoadData.GE_Recon.fidall.grad_warp( ...
            protonSOS,geInfo.h,gradcoil,true,'',zfov);
        gradwarpApplied = true;
    catch ME
        warning('GE gradwarp failed; keeping uncorrected proton image: %s',ME.message);
    end
end

if useGEOrientation
    protonSOS = applyNativeGEProtonOrientation(protonSOS,geInfo.h);
    if ~isempty(singleCoilComplex)
        singleCoilComplex = applyNativeGEProtonOrientation(singleCoilComplex,geInfo.h);
    end
    if keepCoils
        for ic=1:nCoils
            coilsComplex(:,:,:,ic) = applyNativeGEProtonOrientation( ...
                coilsComplex(:,:,:,ic),geInfo.h);
        end
    end
end

ProtonImageRaw = double(abs(protonSOS));
ProtonImageRaw=flip(rot90(ProtonImageRaw,-1),1);
% ProtonImageRaw=flip(ProtonImageRaw,3);
ProtonImageHR = ProtonImageRaw;
ProtonImageLR = medfilt3(ProtonImageRaw,[7 7 7]);
ProtonImageLR = imgaussfilt3(ProtonImageLR,0.5);
% figure; imslice(ProtonImageRaw)
% Preserve original GE proton N4 correction workflow, with an explicit
% success flag; failure must not silently imply correction was performed.
protonBiasApplied = false;
try
    [ProtonImageLR,~] = ...
        GasExchangeFunctions.Dissolved_ProtonBiasCorrection(ProtonImageLR);
    [ProtonImageHR,~] = ...
        GasExchangeFunctions.Dissolved_ProtonBiasCorrection(ProtonImageHR);
    protonBiasApplied = true;
catch ME
    warning('GE proton bias correction failed: %s',ME.message);
end

ProtonMax = prctile(abs(ProtonImageLR(:)),99.99);

%% Save Proton Montage

ProtonMontage = figure('Name','Proton Image');
set(ProtonMontage,'WindowState','minimized');

% Display reconstructed proton image before bias correction
montage(abs(ProtonImageRaw),'DisplayRange',[0 ProtonMax]);

% Save in the reconstruction output folder
if ~isfolder(MainInput.OutputPath)
    mkdir(MainInput.OutputPath);
end

savefig(ProtonMontage, ...
    fullfile(MainInput.OutputPath,'ProtonMontage.fig'));

close(ProtonMontage);

%% Store results: original GE Proton output contract

Proton = struct();
Proton.Image = double(ProtonImageHR);
Proton.ProtonImageLR = double(ProtonImageLR);
Proton.ProtonImageHR = double(ProtonImageHR);
Proton.ProtonMax = double(ProtonMax);
Proton.H_RecMatrix = H_RecMatrix;

% Original GE convention: this is the reconstructed coil-combined magnitude
% *before* LR smoothing/N4; it is not a phase-resolved combined image.
Proton.ImageCombined = double(ProtonImageRaw);
if nCoils==1
    Proton.ImageComplex = singleCoilComplex;
else
    Proton.ImageComplex = [];
end
if keepCoils
    Proton.CoilImagesComplex = coilsComplex;
end
Proton.H_AcqMatrix = geInfo.mtx_acq(1);
Proton.H_nsamp = nReadout;
Proton.H_nprof = nProjections;
Proton.H_coils = nCoils;
Proton.H_interleaves = geInfo.ninterleaves;
Proton.HTraj = HTraj;
Proton.PixelShift = PixelShift;
Proton.GENativeShift = geInfo.shft;
Proton.HPhaseDynamics = phaseDynamics;
Proton.Motion = motionInfo;
Proton.MotionKeep = motionKeep;
Proton.MotionRejectedProjections = nnz(~motionKeep);
Proton.ReconstructionMethod = 'GE raw2grid + MRD 3D LSQGridded iterative DCF + SOS';
Proton.MRDReconMode = reconMode;
Proton.TrajectoryUnits = 'normalized cycles-per-voxel (fractional k-space)';
Proton.GradwarpApplied = gradwarpApplied;
Proton.GEOrientationApplied = useGEOrientation;
Proton.ProtonBiasCorrectionApplied = protonBiasApplied;
Proton.CroppedFromOvergrid = true;

if logical(getOption(MainInput,'HKeepKspace',false))
    Proton.HKSpace = HKSpace;
end

% Preserve original source file and parameter metadata
Proton.DataFile = d;
Proton.HeaderFile = h;
Proton.WaveformFile = wfn;
Proton.WaveformFileName = waveformFileName;
Proton.MatlabScripterFile = scripterPath;
Proton.OriginalWaveformPath = waveformOriginalPath;
Proton.ReconCommand = ['recon_grid_interp(' commandArguments ');'];
Proton.mtx_reco = mtx_reco;
Proton.interp_factor = interp_factor;
Proton.delay = delay;
Proton.lb = lb;
Proton.mask_fac = mask_fac;
Proton.coco = coco;
Proton.plt = plt;
Proton.grdwrp = grdwrp;
Proton.AcquisitionType = '3D_Radial';

if isfield(MainInput,'HFileName') && ~isempty(MainInput.HFileName)
    Proton.filename = MainInput.HFileName;
else
    Proton.filename = [dataName dataExt];
end
if isfield(MainInput,'HDataLocation') && ~isempty(MainInput.HDataLocation)
    Proton.folder = MainInput.HDataLocation;
else
    Proton.folder = dataFolder;
end
if isfield(MainInput,'OutputPath') && ~isempty(MainInput.OutputPath)
    Proton.outputfolder = MainInput.OutputPath;
else
    Proton.outputfolder = dataFolder;
end
Proton.FileName = Proton.filename;
Proton.DataLocation = Proton.folder;
if isfield(MainInput,'Scanner')
    Proton.Scanner = MainInput.Scanner;
else
    Proton.Scanner = 'GE';
end
if isfield(MainInput,'Institute')
    Proton.Institute = MainInput.Institute;
else
    Proton.Institute = '';
end
Proton.ImageSize = size(Proton.Image);
Proton.CombinedImageSize = size(Proton.ImageCombined);

fprintf('\nGE MRD-style proton reconstruction completed: [%s], %d coils.\n', ...
    num2str(Proton.ImageSize),nCoils);

if logical(getOption(MainInput,'HSaveQC',false))
    try
        if ~isfolder(Proton.outputfolder),mkdir(Proton.outputfolder);end
        qcFig = figure('Color','w','Name','GE Proton MRD-style QC','Visible','off');
        montage(abs(Proton.Image),'DisplayRange',[]);
        title('GE Proton | MRD-style');
        savefig(qcFig,fullfile(Proton.outputfolder,'GE_Proton_MRDStyle.fig'));
        exportgraphics(qcFig,fullfile(Proton.outputfolder,'GE_Proton_MRDStyle.png'), ...
            'Resolution',300);
        close(qcFig);
    catch ME
        warning('Unable to save proton QC figures: %s',ME.message);
    end
end
fprintf('============================================================\n\n');
end

%% Local helper functions

function [HKSpace,HTraj,info] = extractGERadialProton(scanFile,wfn,mtx_reco,interp_factor,delay,lb)
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
shft=shft + 0.5; % native GE proton half-voxel shift

k=double(k);
dd=LoadData.GE_Recon.fidall.raw2grid(d, ...
    LoadData.GE_Recon.fidall.ischop(h),k,shft,[],[],ind,delay*1e-6,bw,false);
clear d

% Native optional apodization (retained because user asked to preserve GE
% preprocessing); MRD kernel/DCF are subsequently used for reconstruction.
if any(abs(lb)>0)
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

fprintf('GE proton extracted: [%d samples x %d projections x %d coils].\n', ...
    nRead,nInterleaves,ncoils);
end

function [keep,phaseDynamics,info] = findProtonMotionCutoff(K,T,enabled)
[nRead,nProj,nCoils]=size(K);
rad=squeeze(sqrt(sum(T.^2,1)));
[~,idx0]=min(rad,[],1);
centerSignals=complex(zeros(nProj,nCoils));
for ic=1:nCoils
    v=K(:,:,ic);
    centerIdx=sub2ind([nRead nProj],idx0,1:nProj);
    centerSignals(:,ic)=reshape(v(centerIdx),[],1);
end
phaseByCoil=unwrap(angle(centerSignals),[],1);
phaseRef=median(phaseByCoil(1:max(1,floor(nProj/3)),:),1);
phaseDynamics=rad2deg(median(phaseByCoil-phaseRef,2));
phaseSmoothed=movmedian(phaseDynamics,15,'Endpoints','shrink');
keep=true(nProj,1);
info=struct('enabled',enabled,'cutoff',nProj+1,'threshold_deg',NaN, ...
    'reason','disabled','referencePhase_deg',NaN);
if ~enabled || nProj < 90
    if enabled && nProj<90,info.reason='too few projections for cutoff';end
    return;
end
firstThird=max(10,floor(nProj/3));
refRange=max(1,round(nProj*0.05)):firstThird;
refPhase=median(phaseSmoothed(refRange));
rmse=sqrt(mean((phaseSmoothed(refRange)-refPhase).^2));
threshold=max(5,3*rmse);
info.threshold_deg=threshold;
info.referencePhase_deg=refPhase;

bad=abs(phaseSmoothed-refPhase)>threshold;
minConsecutive=max(15,round(0.015*nProj));
cutoff=[];
for ii=firstThird+1:nProj-minConsecutive+1
    if all(bad(ii:ii+minConsecutive-1))
        cutoff=ii;
        break;
    end
end
if ~isempty(cutoff)
    keep(cutoff:end)=false;
    info.cutoff=cutoff;
    info.reason='persistent phase change after reference window';
else
    info.reason='no persistent late phase change';
end
end

function v = centerCrop3D(v,N)
sz=size(v);
assert(numel(sz)>=3 && all(sz(1:3)>=N),'3D crop size exceeds reconstructed volume.');
a=floor((sz(1:3)-N)/2)+1;
v=v(a(1):a(1)+N-1,a(2):a(2)+N-1,a(3):a(3)+N-1);
end

function v = undoMRDRotation(v)
% MRD helper ends with rot90(flip(Image_Out)): invert exactly.
% v=flip(rot90(v,-1),1);
% v=flip(rot90(v,-1),1);
end

function v = applyNativeGEProtonOrientation(v,h)
% Replicate GE recon_grid_interp post-gridding orientation plus original
% GE_Proton_3D_Radial_Recon wrapper's orientation.
v=fliplr(v);
v=LoadData.GE_Recon.fidall.mri_orientation(v, ...
    h.data_acq_tab.rotate,h.data_acq_tab.transpose);
v=permute(v,[2 1 3]);
v=rot90(v,1); % equivalent to imrotate(...,90) without interpolation
v=flip(v,1);
end

function y = getOption(s,field,default)
if isfield(s,field) && ~isempty(s.(field))
    y=s.(field);
else
    y=default;
end
end

function tf = islogicalOrNumericTrue(v)
tf=(isnumeric(v)||islogical(v)) && ~isempty(v) && any(v(:)~=0);
end


function value = parseMatlabValue(arg)

% Parse numeric/scalar/vector/logical values extracted from the
% matlab_scripter command without using eval.

arg = strtrim(arg);

if isempty(arg) || strcmp(arg,'[]')
    value = [];
    return;
end

if strcmpi(arg,'true')
    value = true;
    return;
elseif strcmpi(arg,'false')
    value = false;
    return;
end

% Scalar numeric
scalarValue = str2double(arg);

if ~isnan(scalarValue)
    value = scalarValue;
    return;
end

% Numeric vector/matrix such as:
% [0 0 0]
% [1 0 1 16]
% [0 1]

if startsWith(arg,'[') && endsWith(arg,']')

    content = strtrim(arg(2:end-1));

    if isempty(content)
        value = [];
        return;
    end

    rows = strsplit(content,';');
    parsedRows = cell(size(rows));

    for rr = 1:numel(rows)

        rowText = strtrim(rows{rr});

        % Allow commas or spaces
        rowText = strrep(rowText,',',' ');

        nums = sscanf(rowText,'%f').';

        if isempty(nums)
            error('Unable to parse reconstruction parameter: %s',arg);
        end

        parsedRows{rr} = nums;

    end

    try
        value = vertcat(parsedRows{:});
    catch
        error('Unable to parse reconstruction matrix/vector: %s',arg);
    end

    return;
end

error('Unsupported reconstruction parameter in matlab_scripter: %s',arg);

end