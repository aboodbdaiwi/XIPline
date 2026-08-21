% Run a calibration - will generate files that can be used for
% calibration based processing (ASSET, PURE, MCSI).

directory = fileparts(mfilename('fullpath'));
pfile = fullfile(directory, 'Data/Calibration/P39424.7');

GERecon('Calibration.Process', pfile);