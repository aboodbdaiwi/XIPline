% Call AssetRecon script with the calibration and non-calibration Pfiles.
% Data is relative to this script, so setup filepaths so the script can be
% called from anywhere.

directory = fileparts(mfilename('fullpath'));
pfile = fullfile(directory, 'Data/Asset/P94208.7');
calibration = fullfile(directory, 'Data/Asset/Calibration/P89088.7');

AssetRecon(pfile, calibration);
