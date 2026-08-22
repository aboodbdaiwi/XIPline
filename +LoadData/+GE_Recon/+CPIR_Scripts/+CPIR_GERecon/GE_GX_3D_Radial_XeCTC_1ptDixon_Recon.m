function GasExchange = GE_GX_3D_Radial_XeCTC_1ptDixon_Recon(MainInput)

% GE_GX_3D_Radial_XeCTC_1ptDixon_Recon
% Reconstruct GE 129Xe gas-exchange 3D radial XeCTC 1-point Dixon data for XIPline.
%
% Input:
%   MainInput.XeFullPath     - Full path to GE ScanArchive .h5 file
%   MainInput.XeDataLocation - Data folder (optional)
%   MainInput.XeFileName     - Data filename (optional)
%
% Output:
%   GasExchange - structure containing reconstructed gas-exchange images,
%                 reconstruction settings, source files, and metadata.

if nargin < 1 || isempty(MainInput)
    error('MainInput is required.');
end

if ~isfield(MainInput,'XeFullPath') || isempty(MainInput.XeFullPath)
    error('MainInput.XeFullPath is missing or empty.');
end

d = MainInput.XeFullPath;

if ~isfile(d)
    error('GE gas-exchange data file was not found:\n%s',d);
end

[dataFolder,dataName,dataExt] = fileparts(d);

fprintf('\n============================================================\n');
fprintf('GE 129Xe Gas Exchange 3D Radial XeCTC 1pt Dixon Reconstruction\n');
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

%% Extract recon_grid_interp command

commandToken = regexp( ...
    scripterText, ...
    'recon_grid3d_gx\s*\((.*?)\)\s*;', ...
    'tokens', ...
    'once');

if isempty(commandToken)
    error(['Could not find a valid recon_grid3d_gx(...) command in:\n' ...
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

% Expected recon_grid3d_gx arguments:
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
% XIPline will manage the returned reconstruction.
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

%% Perform reconstruction

fprintf('\nReconstructing GE gas-exchange 3D radial XeCTC 1pt Dixon data...\n');

[bb,bbabs,gx] = ...
    LoadData.GE_Recon.CPIR_Scripts.CPIR_GERecon.recon_grid3d_gx_XIPline( ...
    d, ...
    h, ...
    wfn, ...
    mtx_reco, ...
    interp_factor, ...
    delay, ...
    lb, ...
    fname, ...
    mask_fac, ...
    coco, ...
    plt, ...
    grdwrp);

fprintf('GE gas-exchange 3D radial XeCTC 1pt Dixon reconstruction completed.\n');

%% Store results

GasExchange = struct();

%% Store reconstructed gas and dissolved data separately

% recon_grid3d_gx returns the two XeCTC 1-point-Dixon image streams in
% dimension 4.  Determine which stream is gas from its much larger total
% signal, rather than hard-coding odd/even ordering.
bbSize = size(bbabs);

if numel(bbSize) < 4 || bbSize(4) ~= 2
    error(['Expected the GE CPIR gas-exchange reconstruction to contain ' ...
           'exactly two image streams in dimension 4, but size(bbabs) = [%s].'], ...
           num2str(bbSize));
end

streamSignal = zeros(1,2);
for tt = 1:2
    tmp = double(abs(bbabs(:,:,:,tt)));
    streamSignal(tt) = sum(tmp(:));
end

[~,gasIdx] = max(streamSignal);
dissolvedIdx = 3-gasIdx;

GasExchange.GasStreamIndex = gasIdx;
GasExchange.DissolvedStreamIndex = dissolvedIdx;
GasExchange.StreamSignal = streamSignal;

% Extract complex images before taking magnitude.  The downstream
% gas-exchange analysis uses phase information from the dissolved image.
GasImage = squeeze(bb(:,:,:,gasIdx,:));
DissolvedImage = squeeze(bb(:,:,:,dissolvedIdx,:));

GasImageCombined = squeeze(bbabs(:,:,:,gasIdx));
DissolvedImageCombined = squeeze(bbabs(:,:,:,dissolvedIdx));

% Apply the same XIPline orientation to both components.
GasImage = orientGXVolume(GasImage);
DissolvedImage = orientGXVolume(DissolvedImage);
GasImageCombined = orientGXVolume(GasImageCombined);
DissolvedImageCombined = orientGXVolume(DissolvedImageCombined);

% Store the directly reconstructed images.
GasExchange.GasImage = GasImage;
GasExchange.DissolvedImage = DissolvedImage;
GasExchange.GasImageCombined = double(abs(GasImageCombined));
GasExchange.DissolvedImageCombined = double(abs(DissolvedImageCombined));

% The uncorrected ventilation image is the reconstructed gas-phase image.
GasExchange.UncorrectedVentImage = GasImage;

% Apply the same ventilation bias-correction approach used by the Philips
% gas-exchange pipeline. If it cannot be applied, retain the original gas
% image rather than silently failing the reconstruction.
try
    VentMask = abs(GasImage) > (prctile(abs(GasImage(:)),95)/3);
    VentMask = imerode(VentMask,strel('sphere',5));

    Xeregions = regionprops3(logical(VentMask),'all');
    Xeregions = sortrows(Xeregions,1,'descend');

    VentMask2 = false(size(VentMask));

    if ~isempty(Xeregions)
        VentMask2(Xeregions.VoxelIdxList{1,1}) = true;

        if size(Xeregions,1) > 1 && ...
                Xeregions.Volume(2) > 0.33*Xeregions.Volume(1)
            VentMask2(Xeregions.VoxelIdxList{2,1}) = true;
        end

        VentMask2 = logical(imdilate(VentMask2,strel('sphere',6)));
        [VentImage,~] = GasExchangeFunctions.Vent_BiasCorrection( ...
            GasImage,VentMask2);
        GasExchange.VentMask = VentMask2;
    else
        VentImage = GasImage;
        GasExchange.VentMask = true(size(GasImage));
    end
catch ME
    warning('GE ventilation bias correction failed: %s',ME.message);
    VentImage = GasImage;
    GasExchange.VentMask = true(size(GasImage));
end

GasExchange.VentImage = VentImage;

% Generic XIPline image field remains the bias-corrected ventilation
% magnitude image.
GasExchange.Image = double(abs(VentImage));
GasExchange.ImageCombined = double(abs(GasImageCombined));

% Preserve the complete reconstructed complex two-stream dataset.
GasExchange.ImageComplex = bb;

%% Store GE raw trajectory and k-space needed by gas-exchange analysis

GasExchange.XeTraj = gx.XeTraj;

if ~isempty(gx.Interleave1KSpace) && ~isempty(gx.Interleave2KSpace)

    kSignal = [ ...
        mean(abs(gx.Interleave1KSpace(1,:))), ...
        mean(abs(gx.Interleave2KSpace(1,:)))];

    [~,gasKIdx] = max(kSignal);
    dissKIdx = 3-gasKIdx;

    if gasKIdx == 1
        GasKSpace = gx.Interleave1KSpace;
        DissolvedKSpace = gx.Interleave2KSpace;
    else
        GasKSpace = gx.Interleave2KSpace;
        DissolvedKSpace = gx.Interleave1KSpace;
    end

    GasExchange.GasKSpace = GasKSpace;
    GasExchange.DissolvedKSpace = DissolvedKSpace;

    % Remove the approach to steady state, matching the Philips pipeline.
    SS_ind = min(60,max(0,size(GasKSpace,2)-1));

    GasKSpace_SS = GasKSpace(:,SS_ind+1:end);
    DissolvedKSpace_SS = DissolvedKSpace(:,SS_ind+1:end);
    XeTraj_SS = gx.XeTraj(:,:,SS_ind+1:end);

    GasExchange.SS_ind = SS_ind;
    GasExchange.GasKSpace_SS = GasKSpace_SS;
    GasExchange.DissolvedKSpace_SS = DissolvedKSpace_SS;
    GasExchange.XeTraj_SS = XeTraj_SS;

    GasExchange.XeImg_nsamp = size(GasKSpace,1);
    GasExchange.Xe_nprof = size(GasKSpace,2);
    GasExchange.Xe_interleaves = 2;
    GasExchange.XeOrder = 1:size(GasKSpace,2);

else
    warning(['GE raw interleaves could not be reduced to [readout x projection]. ' ...
             'Raw arrays are still available in GasExchange.GERaw.']);

    GasExchange.GasKSpace = [];
    GasExchange.DissolvedKSpace = [];
    GasExchange.GasKSpace_SS = [];
    GasExchange.DissolvedKSpace_SS = [];
    GasExchange.XeTraj_SS = [];
    GasExchange.SS_ind = [];
    GasExchange.XeImg_nsamp = gx.nacq;
    GasExchange.Xe_nprof = gx.nproj;
    GasExchange.Xe_interleaves = 2;
    GasExchange.XeOrder = [];
end

GasExchange.GERaw = gx;

%% Match the Philips gas-exchange metadata/output contract where supported

% GE image positioning is already handled in recon_grid3d_gx through the
% header-based spatial shift and mri_orientation. No additional Philips
% pixel recentering should be applied.
GasExchange.PixelShift = [0;0;0];

% GE header TE is stored in microseconds. Store TE90/ActTE90 in ms to match
% the Philips GasExchange structure.
try
    TE90 = double(gx.header.rdb_hdr.te)*1e-3;
catch
    TE90 = NaN;
end
GasExchange.TE90 = TE90;
GasExchange.ActTE90 = TE90;

% Use the prescribed GE flip angle when available.
try
    DisFlipAngle = double(gx.header.image.mr_flip);
catch
    DisFlipAngle = NaN;
end
GasExchange.DisFlipAngle = DisFlipAngle;

% XeCTC 1-point Dixon frequency jump used by the CCHMC implementation.
% This matches the CCHMC value used in the Philips pipeline.
GasExchange.freq_jump = 7143;

% Dwell time from the GE spectral width.
try
    GasExchange.dwell_s = 1/double(gx.header.rdb_hdr.spectral_width);
catch
    GasExchange.dwell_s = NaN;
end

% The GE reconstruction applies the physical trajectory/header FOV shift
% internally; expose the effective reconstruction matrix.
GasExchange.Xe_RecMatrix = gx.mtx_reco;

%% Load GE calibration spectroscopy for gas-exchange processing

% The GE XeCTC gas-exchange image reconstruction itself does not provide the
% dissolved-phase spectral model needed for gas-contamination correction,
% RBC/barrier reference information, and later gas-exchange analysis.
%
% Use the calibration dataset supplied through MainInput:
%
%   MainInput.CalFullPath
%   MainInput.CalDataLocation
%   MainInput.CalFileName
%   MainInput.Cal_name
%   MainInput.CalDataext

GasExResults = [];
CalResults = [];
AppendedDissolvedNMRFit = [];

if isfield(MainInput,'CalFullPath') && ...
        ~isempty(MainInput.CalFullPath) && ...
        isfile(MainInput.CalFullPath)

    fprintf('\nAnalyzing GE calibration spectroscopy...\n');

    try
        [GasExResults, CalResults] = ...
            LoadData.GE_Recon.CPIR_Scripts.CPIR_GERecon.GE_Calibration_Spect_CCHMC( ...
            MainInput);

        fprintf('GE calibration spectroscopy completed.\n');

    catch ME
        warning(['GE calibration spectroscopy failed. Gas-exchange images ' ...
                 'will still be returned, but calibration-dependent ' ...
                 'corrections will be skipped.\n%s'],ME.message);
    end

else
    warning(['MainInput.CalFullPath was not supplied or the calibration ' ...
             'file does not exist. Calibration-dependent GE gas-exchange ' ...
             'processing will be skipped.']);
end

% Store the complete calibration outputs for downstream XIPline analysis.
GasExchange.GasExResults = GasExResults;
GasExchange.CalResults = CalResults;

% Use the calibration dissolved spectral fit as the GE equivalent of the
% Philips AppendedDissolvedNMRFit.
if ~isempty(GasExResults) && ...
        isstruct(GasExResults) && ...
        isfield(GasExResults,'DisFit') && ...
        ~isempty(GasExResults.DisFit)

    AppendedDissolvedNMRFit = GasExResults.DisFit;

    % Expose calibration quantities directly in GasExchange as well.
    if isfield(GasExResults,'RbcBarRatio')
        GasExchange.RbcBarRatio = GasExResults.RbcBarRatio;
    else
        GasExchange.RbcBarRatio = NaN;
    end

    if isfield(GasExResults,'GasDisRatio')
        GasExchange.GasDisRatio = GasExResults.GasDisRatio;
    else
        GasExchange.GasDisRatio = NaN;
    end

else
    GasExchange.RbcBarRatio = NaN;
    GasExchange.GasDisRatio = NaN;
end

GasExchange.AppendedDissolvedNMRFit = AppendedDissolvedNMRFit;

% Store useful calibration values without replacing acquisition-specific
% gas-exchange parameters.
if ~isempty(CalResults) && isstruct(CalResults)

    if isfield(CalResults,'te90')
        GasExchange.CalibrationTE90 = CalResults.te90;
    else
        GasExchange.CalibrationTE90 = NaN;
    end

    if isfield(CalResults,'flip_angle')
        GasExchange.CalibrationFlipAngle = CalResults.flip_angle;
    else
        GasExchange.CalibrationFlipAngle = NaN;
    end

    if isfield(CalResults,'FlipScaleFactor')
        GasExchange.FlipScaleFactor = CalResults.FlipScaleFactor;
    else
        GasExchange.FlipScaleFactor = NaN;
    end

    if isfield(CalResults,'freq_target')
        GasExchange.CalibrationFrequencyTarget = CalResults.freq_target;
    else
        GasExchange.CalibrationFrequencyTarget = NaN;
    end

    if isfield(CalResults,'targetTG')
        GasExchange.targetTG = CalResults.targetTG;
    else
        GasExchange.targetTG = NaN;
    end

else
    GasExchange.CalibrationTE90 = NaN;
    GasExchange.CalibrationFlipAngle = NaN;
    GasExchange.FlipScaleFactor = NaN;
    GasExchange.CalibrationFrequencyTarget = NaN;
    GasExchange.targetTG = NaN;
end

% Gas contamination correction requires a dissolved spectral fit. Apply it
% only when that information is available; otherwise preserve the
% uncorrected dissolved k-space/image explicitly.
GasExchange.CorrectedDissKSpace_SS = GasExchange.DissolvedKSpace_SS;
GasExchange.CorrDissolvedImage = DissolvedImage;

if ~isempty(AppendedDissolvedNMRFit) && ...
        ~isempty(GasExchange.DissolvedKSpace_SS) && ...
        ~isempty(GasExchange.GasKSpace_SS)

    try
        gasFlipAngle = 0.5;

        CorrectedDissKSpace_SS = ...
            GasExchangeFunctions.GasPhaseContaminationRemoval( ...
            GasExchange.DissolvedKSpace_SS, ...
            GasExchange.GasKSpace_SS, ...
            GasExchange.dwell_s, ...
            -GasExchange.freq_jump, ...
            AppendedDissolvedNMRFit.phase(3), ...
            AppendedDissolvedNMRFit.area(3), ...
            gasFlipAngle);

        GasExchange.CorrectedDissKSpace_SS = CorrectedDissKSpace_SS;

    catch ME
        warning('GE gas-phase contamination correction failed: %s',ME.message);
    end
end

% These Philips fields require dedicated spectroscopy and/or the Philips
% keyhole reconstruction pathway. Keep the fields present for XIPline
% compatibility, but do not fabricate measurements that are not available
% from the GE image reconstruction.
GasExchange.RBC2Bar_struct = [];
GasExchange.RBCOsc_High_Image = [];
GasExchange.RBCOsc_Low_Image = [];
GasExchange.RBCOsc_Normalization = [];
GasExchange.DissolvedNMR = [];
GasExchange.SigDynamics = [];
GasExchange.XeSin = [];

% Reconstruction settings
GasExchange.mtx_reco = mtx_reco;
GasExchange.interp_factor = interp_factor;
GasExchange.delay = delay;
GasExchange.lb = lb;
GasExchange.mask_fac = mask_fac;
GasExchange.coco = coco;
GasExchange.plt = plt;
GasExchange.grdwrp = grdwrp;

% Source files
GasExchange.DataFile = d;
GasExchange.HeaderFile = h;
GasExchange.WaveformFile = wfn;
GasExchange.WaveformFileName = waveformFileName;
GasExchange.MatlabScripterFile = scripterPath;

% Preserve original GE paths/settings for traceability
GasExchange.OriginalWaveformPath = waveformOriginalPath;
GasExchange.ReconCommand = ...
    ['recon_grid3d_gx(' commandArguments ');'];

% Dimensions
GasExchange.ImageSize = size(GasExchange.Image);
GasExchange.CombinedImageSize = size(GasExchange.ImageCombined);

if isfield(MainInput,'OutputPath')
    GasExchange.outputpath = MainInput.OutputPath;
end

% Calibration source metadata
if isfield(MainInput,'CalFullPath')
    GasExchange.CalFullPath = MainInput.CalFullPath;
else
    GasExchange.CalFullPath = '';
end

if isfield(MainInput,'CalDataLocation')
    GasExchange.CalDataLocation = MainInput.CalDataLocation;
else
    GasExchange.CalDataLocation = '';
end

if isfield(MainInput,'CalFileName')
    GasExchange.CalFileName = MainInput.CalFileName;
else
    GasExchange.CalFileName = '';
end

if isfield(MainInput,'Cal_name')
    GasExchange.Cal_name = MainInput.Cal_name;
else
    GasExchange.Cal_name = '';
end

if isfield(MainInput,'CalDataext')
    GasExchange.CalDataext = MainInput.CalDataext;
else
    GasExchange.CalDataext = '';
end

% MainInput metadata
if isfield(MainInput,'Scanner')
    GasExchange.Scanner = MainInput.Scanner;
else
    GasExchange.Scanner = 'GE';
end

if isfield(MainInput,'Institute')
    GasExchange.Institute = MainInput.Institute;
else
    GasExchange.Institute = '';
end

if isfield(MainInput,'XeFileName')
    GasExchange.FileName = MainInput.XeFileName;
else
    GasExchange.FileName = [dataName dataExt];
end

if isfield(MainInput,'XeDataLocation')
    GasExchange.DataLocation = MainInput.XeDataLocation;
else
    GasExchange.DataLocation = dataFolder;
end

fprintf('\nReconstructed gas-exchange image size:\n');
disp(size(GasExchange.Image));

fprintf('============================================================\n\n');

end


function vol = orientGXVolume(vol)

% Apply the XIPline GE orientation to the first three spatial dimensions.
% Input is expected to be a 3-D complex/magnitude volume.

vol = squeeze(vol);

if ndims(vol) ~= 3
    error('Expected a 3-D gas-exchange image volume; received size [%s].', ...
        num2str(size(vol)));
end

vol = permute(vol,[2,1,3]);
vol = imrotate(vol,-90);
vol = flip(vol,1);

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