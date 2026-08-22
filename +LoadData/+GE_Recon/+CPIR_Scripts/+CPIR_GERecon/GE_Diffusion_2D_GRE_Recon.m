function Diffusion = GE_Diffusion_2D_GRE_Recon(MainInput)

% GE_Diffusion_Recon
% Reconstruct GE hyperpolarized 129Xe diffusion MRI data for XIPline.
%
% Input:
%   MainInput.XeFullPath     - Full path to GE ScanArchive .h5 file
%   MainInput.XeDataLocation - Data folder (optional)
%   MainInput.XeFileName     - Data filename (optional)
%
% Output:
%   Diffusion - structure containing reconstructed diffusion images,
%               reconstruction settings, source files, and metadata.

if nargin < 1 || isempty(MainInput)
    error('MainInput is required.');
end

if ~isfield(MainInput,'XeFullPath') || isempty(MainInput.XeFullPath)
    error('MainInput.XeFullPath is missing or empty.');
end

d = MainInput.XeFullPath;

if ~isfile(d)
    error('GE diffusion data file was not found:\n%s',d);
end

[dataFolder,dataName,dataExt] = fileparts(d);

fprintf('\n============================================================\n');
fprintf('GE 129Xe Diffusion Reconstruction\n');
fprintf('============================================================\n');
fprintf('Data file:\n%s\n',d);

% -------------------------------------------------------------------------
% Find corresponding GE header file
% -------------------------------------------------------------------------

expectedHeader = fullfile(dataFolder,[dataName '.hdr']);

if isfile(expectedHeader)
    h = expectedHeader;
else

    headerFiles = dir(fullfile(dataFolder,'*.hdr'));

    if isempty(headerFiles)
        error('No GE .hdr file was found in:\n%s',dataFolder);
    end

    % First try to find a header containing the same ScanArchive basename
    matchIndex = [];

    for ii = 1:numel(headerFiles)

        [~,headerName] = fileparts(headerFiles(ii).name);

        if strcmpi(headerName,dataName)
            matchIndex = ii;
            break;
        end

    end

    if isempty(matchIndex)

        % If there is exactly one header in the folder, use it
        if numel(headerFiles) == 1
            matchIndex = 1;
        else
            error(['Could not uniquely identify the GE header file.\n' ...
                   'Expected:\n%s'],expectedHeader);
        end

    end

    h = fullfile(headerFiles(matchIndex).folder, ...
                 headerFiles(matchIndex).name);

end

fprintf('Header file:\n%s\n',h);

% -------------------------------------------------------------------------
% Find matlab_scripter file
% -------------------------------------------------------------------------

scripterFiles = dir(fullfile(dataFolder,'*matlab_scripter*'));

if isempty(scripterFiles)
    error(['Could not find a file containing "matlab_scripter" in:\n%s'], ...
        dataFolder);
end

% If several matlab_scripter files exist, select the one containing
% the current ScanArchive filename.
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

% If no file explicitly references this ScanArchive, only accept the
% result automatically when there is exactly one scripter file.
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

% -------------------------------------------------------------------------
% Extract recon_cart_diffusion command
% -------------------------------------------------------------------------

commandToken = regexp( ...
    scripterText, ...
    'recon_cart_diffusion\s*\((.*?)\)\s*;', ...
    'tokens', ...
    'once');

if isempty(commandToken)
    error(['Could not find a valid recon_cart_diffusion(...) command in:\n' ...
           '%s'],scripterPath);
end

commandArguments = strtrim(commandToken{1});

% -------------------------------------------------------------------------
% Split MATLAB function arguments while respecting quoted paths
% -------------------------------------------------------------------------

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

if numel(args) < 8
    error(['The recon_cart_diffusion command does not contain the expected ' ...
           '8 arguments.\n\nParsed command:\n%s'], ...
           commandArguments);
end

% Expected:
%
% 1 = d
% 2 = h
% 3 = waveform
% 4 = zf
% 5 = delay
% 6 = lb
% 7 = fname
% 8 = comb_time

% -------------------------------------------------------------------------
% Extract waveform filename from original GE Linux path
% -------------------------------------------------------------------------

waveformArg = strtrim(args{3});

if startsWith(waveformArg,'''') && endsWith(waveformArg,'''')
    waveformOriginalPath = waveformArg(2:end-1);
else
    waveformOriginalPath = waveformArg;
end

% fileparts on Windows will not reliably interpret "/" from a Linux path,
% so normalize it before extracting the filename.
waveformNormalized = strrep(waveformOriginalPath,'\','/');

slashIndex = find(waveformNormalized == '/',1,'last');

if isempty(slashIndex)
    waveformFileName = waveformNormalized;
else
    waveformFileName = waveformNormalized(slashIndex+1:end);
end

% The exported waveform is expected to be in the SAME directory as the
% ScanArchive data.
wfn = fullfile(dataFolder,waveformFileName);

if ~isfile(wfn)

    % Case-insensitive fallback search
    waveformFiles = dir(fullfile(dataFolder,'*.mat'));

    waveformIndex = [];

    for ii = 1:numel(waveformFiles)

        if strcmpi(waveformFiles(ii).name,waveformFileName)
            waveformIndex = ii;
            break;
        end

    end

    if isempty(waveformIndex)
        error(['Waveform file was specified by matlab_scripter as:\n%s\n\n' ...
               'but could not be found in the local data folder:\n%s'], ...
               waveformFileName,dataFolder);
    end

    wfn = fullfile( ...
        waveformFiles(waveformIndex).folder, ...
        waveformFiles(waveformIndex).name);

end

fprintf('Waveform file:\n%s\n',wfn);

% -------------------------------------------------------------------------
% Extract reconstruction parameters
% -------------------------------------------------------------------------

zf = str2double(strtrim(args{4}));
delay = str2double(strtrim(args{5}));
lb = str2double(strtrim(args{6}));
comb_time = str2double(strtrim(args{8}));

if isnan(zf)
    error('Unable to extract zf from matlab_scripter.');
end

if isnan(delay)
    error('Unable to extract delay from matlab_scripter.');
end

if isnan(lb)
    error('Unable to extract lb from matlab_scripter.');
end

if isnan(comb_time)
    error('Unable to extract comb_time from matlab_scripter.');
end

% -------------------------------------------------------------------------
% fname
%
% The fname contained in matlab_scripter points to the original Linux
% reconstruction computer and therefore should NOT be passed directly.
%
% For XIPline we do not need recon_cart_diffusion to independently save
% PNG/MAT/DICOM files, so use [] here. The reconstruction results are
% returned directly and stored in Diffusion.
% -------------------------------------------------------------------------

fname = [];

fprintf('\nReconstruction settings:\n');
fprintf('  zf        = %g\n',zf);
fprintf('  delay     = %g\n',delay);
fprintf('  lb        = %g\n',lb);
fprintf('  comb_time = %g\n',comb_time);

% -------------------------------------------------------------------------
% Perform GE diffusion reconstruction
% -------------------------------------------------------------------------

fprintf('\nReconstructing GE diffusion data...\n');

[bb,bbabs] = ...
    LoadData.GE_Recon.CPIR_Scripts.CPIR_GERecon.recon_cart_diffusion( ...
    d, ...
    h, ...
    wfn, ...
    zf, ...
    delay, ...
    lb, ...
    fname, ...
    comb_time);

fprintf('GE diffusion reconstruction completed.\n');

% -------------------------------------------------------------------------
% Store results for XIPline
% -------------------------------------------------------------------------

Diffusion = struct();

% Main reconstructed image requested by XIPline
Diffusion.Image = double(abs(bb));

% Coil-combined reconstruction from recon_cart_diffusion
Diffusion.ImageCombined = double(abs(bbabs));

% Preserve complex reconstruction if it is needed later for advanced
% diffusion processing. Remove this field if memory becomes an issue.
Diffusion.ImageComplex = bb;

% Reconstruction parameters
Diffusion.zf = zf;
Diffusion.delay = delay;
Diffusion.lb = lb;
Diffusion.comb_time = comb_time;

% Source files
Diffusion.DataFile = d;
Diffusion.HeaderFile = h;
Diffusion.WaveformFile = wfn;
Diffusion.WaveformFileName = waveformFileName;
Diffusion.MatlabScripterFile = scripterPath;

% Preserve the original waveform path recorded by the GE scanner
Diffusion.OriginalWaveformPath = waveformOriginalPath;

% Preserve the complete original reconstruction command
Diffusion.ReconCommand = ...
    ['recon_cart_diffusion(' commandArguments ');'];

% Basic reconstruction dimensions
Diffusion.ImageSize = size(Diffusion.Image);
Diffusion.CombinedImageSize = size(Diffusion.ImageCombined);

% Scanner / dataset metadata available from MainInput
if isfield(MainInput,'Scanner')
    Diffusion.Scanner = MainInput.Scanner;
else
    Diffusion.Scanner = 'GE';
end

if isfield(MainInput,'Institute')
    Diffusion.Institute = MainInput.Institute;
else
    Diffusion.Institute = '';
end

if isfield(MainInput,'XeFileName')
    Diffusion.FileName = MainInput.XeFileName;
else
    Diffusion.FileName = [dataName dataExt];
end

if isfield(MainInput,'XeDataLocation')
    Diffusion.DataLocation = MainInput.XeDataLocation;
else
    Diffusion.DataLocation = dataFolder;
end

fprintf('\nReconstructed image size:\n');
disp(size(Diffusion.Image));

fprintf('============================================================\n\n');

end