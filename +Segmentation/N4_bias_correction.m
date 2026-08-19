function [N4, Bias] = N4_bias_correction(MR, parentPath)

% Potential Inputs
% Set location of exe and temp images
%location of N4BiasFieldCorrection

N4Path = mfilename('fullpath');
idcs = strfind(N4Path,filesep); %determine location of file separators
N4Path = [N4Path(1:idcs(end)-1),filesep]; %remove file

if isempty(N4Path)
    N4Path = 'C:\XIPline\';
end
maskarray = ones(size(MR));
% RH: was hardcoded to '\' which broke on macOS/Linux paths; use filesep so this works cross-platform
if parentPath(end) ~= filesep
    parentPath = [parentPath, filesep];
end

% Export Image and Mask temporarily for use in N4
niftiwrite(abs(MR+0.001),[parentPath,'Image.nii']); % Image
niftiwrite(abs(maskarray),[parentPath,'Weight.nii']); % Weight

% Run Bias Correction
% general settings (Abdullah)

% RH: N4BiasFieldCorrection.exe is a Windows-only binary and can't run on macOS/Linux.
% On those platforms, run N4 through ANTsPy (same xe_gx conda env used for Registration.ANTs)
% instead of shelling out to a missing native CLI.
if ispc
    cmd = ['"',N4Path,'N4BiasFieldCorrection.exe"',...%run bias correction
        ' -d 3 -i "',parentPath,'Image.nii"',... % set to 3 dimensions and input image of Image.nii
        ' -s 2',... % shrink by factor of 1
        ' -w "',parentPath,'Weight.nii" ',... % import mask called Weight.nii
        ' -o ["',parentPath,'CorrectedImage.nii","',parentPath,'Bias.nii"]']; % output corrected image and bias field
else
    XIPlineRoot = fullfile(getenv('HOME'), 'XIPline');
    registration_path = fullfile(XIPlineRoot, 'Registration');
    pythonExe = Registration.ANTs.getAntsPythonExe(registration_path);
    n4Script = fullfile(N4Path, 'n4_bias_correction.py');
    cmd = sprintf('"%s" "%s" "%s" "%s" %d "%s" "%s"', ...
        pythonExe, n4Script, [parentPath,'Image.nii'], [parentPath,'Weight.nii'], ...
        2, [parentPath,'CorrectedImage.nii'], [parentPath,'Bias.nii']);
end
system(cmd);

%% Import Results
N4 = double(niftiread(fullfile(parentPath,'CorrectedImage.nii')));
Bias = double(niftiread(fullfile(parentPath,'Bias.nii')));

%% Delete Temp Files
delete([parentPath,'Image.nii']);
delete([parentPath,'Weight.nii']);
delete([parentPath,'CorrectedImage.nii']);
delete([parentPath,'Bias.nii']);

end