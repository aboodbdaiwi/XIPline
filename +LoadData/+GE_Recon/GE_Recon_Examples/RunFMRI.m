% fMRI With Dynamic Phase Correction (phase correction coefficients are
% updated throughout the long scan to minimize ghost levels which may
% elevate with temperature changes as the scan progresses)
directory = fileparts(mfilename('fullpath'));
pfile = fullfile(directory, 'Data/fMRI-Asset/P41472.7');
EpiMultiPhaseRecon(pfile);