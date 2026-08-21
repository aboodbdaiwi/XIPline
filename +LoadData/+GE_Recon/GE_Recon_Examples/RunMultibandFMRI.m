% fMRI With Dynamic Phase Correction and Multiband acquisition (phase correction coefficients are
% updated throughout the long scan to minimize ghost levels which may
% elevate with temperature changes as the scan progresses)
directory = fileparts(mfilename('fullpath'));
pfile = fullfile(directory, 'Data/fMRI-Multiband/P16896.7');
EpiMultiPhaseRecon(pfile);
