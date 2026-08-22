function Ventilation = GE_Vent_2D_Spiral_Recon(MainInput)

% GE_Vent_2D_Spiral_Recon
% Reconstruct GE 2D Spiral hyperpolarized 129Xe ventilation data for XIPline.
%
% Input:
%   MainInput.XeFullPath     - Full path to GE ScanArchive .h5 file
%   MainInput.XeDataLocation - Data folder (optional)
%   MainInput.XeFileName     - Data filename (optional)
%
% Output:
%   Ventilation - structure containing reconstructed ventilation images,
%                 reconstruction settings, source files, and metadata.

if nargin < 1 || isempty(MainInput)
    error('MainInput is required.');
end

if ~isfield(MainInput,'XeFullPath') || isempty(MainInput.XeFullPath)
    error('MainInput.XeFullPath is missing or empty.');
end

d = MainInput.XeFullPath;

if ~isfile(d)
    error('GE ventilation data file was not found:\n%s',d);
end

[dataFolder,dataName,dataExt] = fileparts(d);

fprintf('\n============================================================\n');
fprintf('GE 129Xe 2D Spiral Ventilation Reconstruction\n');
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

fprintf('\nReconstructing GE 2D Spiral ventilation data...\n');

[bb,bbabs] = ...
    LoadData.GE_Recon.CPIR_Scripts.CPIR_GERecon.recon_grid_interp( ...
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

fprintf('GE 2D Spiral ventilation reconstruction completed.\n');

%% Store results

Ventilation = struct();

%% Store reconstructed data

% -------------------------------------------------------------------------
% Full reconstructed data
% -------------------------------------------------------------------------
Image = double(abs(bb));

Image = permute(Image,[2,1,3]);
Image = imrotate(Image,90);
Image = flip(Image,1);

Ventilation.Image = Image;
% imslice(Ventilation.Image)
% -------------------------------------------------------------------------
% Coil-combined reconstructed data
% -------------------------------------------------------------------------
ImageCombined = double(abs(bbabs));

ImageCombined = permute(ImageCombined,[2,1,3]);
ImageCombined = imrotate(ImageCombined,90);
ImageCombined = flip(ImageCombined,1);

Ventilation.ImageCombined = ImageCombined;

% -------------------------------------------------------------------------
% Keep original complex reconstruction if needed for future processing
% -------------------------------------------------------------------------
Ventilation.ImageComplex = bb;

% Reconstruction settings
Ventilation.mtx_reco = mtx_reco;
Ventilation.interp_factor = interp_factor;
Ventilation.delay = delay;
Ventilation.lb = lb;
Ventilation.mask_fac = mask_fac;
Ventilation.coco = coco;
Ventilation.plt = plt;
Ventilation.grdwrp = grdwrp;

% Source files
Ventilation.DataFile = d;
Ventilation.HeaderFile = h;
Ventilation.WaveformFile = wfn;
Ventilation.WaveformFileName = waveformFileName;
Ventilation.MatlabScripterFile = scripterPath;

% Preserve original GE paths/settings for traceability
Ventilation.OriginalWaveformPath = waveformOriginalPath;
Ventilation.ReconCommand = ...
    ['recon_grid_interp(' commandArguments ');'];

% Dimensions
Ventilation.ImageSize = size(Ventilation.Image);
Ventilation.CombinedImageSize = size(Ventilation.ImageCombined);

% MainInput metadata
if isfield(MainInput,'Scanner')
    Ventilation.Scanner = MainInput.Scanner;
else
    Ventilation.Scanner = 'GE';
end

if isfield(MainInput,'Institute')
    Ventilation.Institute = MainInput.Institute;
else
    Ventilation.Institute = '';
end

if isfield(MainInput,'XeFileName')
    Ventilation.FileName = MainInput.XeFileName;
else
    Ventilation.FileName = [dataName dataExt];
end

if isfield(MainInput,'XeDataLocation')
    Ventilation.DataLocation = MainInput.XeDataLocation;
else
    Ventilation.DataLocation = dataFolder;
end

fprintf('\nReconstructed ventilation image size:\n');
disp(size(Ventilation.Image));

fprintf('============================================================\n\n');

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