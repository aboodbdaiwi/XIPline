% MCSI - Using McsiCalibration.h5 in Pfile directory for calibration data

directory = fileparts(mfilename('fullpath'));
pfile = fullfile(directory, 'Data/Spectro-MCSI/P08704.7');
calibration = fullfile(directory, 'Data/Spectro-MCSI/McsiCalibration.h5');

SpectroMultiChannelMultiVoxel(pfile, calibration);