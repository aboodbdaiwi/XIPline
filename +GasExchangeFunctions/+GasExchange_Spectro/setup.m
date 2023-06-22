% Get current path
rootDistribDir = pwd();

% Add GE package to path
disp('Adding spectroscopy package to MATLAB path...');

path(genpath([rootDistribDir filesep() 'SpectroscopyPackage']) ,path);