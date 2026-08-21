% Call CartesianRecon script with the specified Pfile
directory = fileparts(mfilename('fullpath'));
pfile = fullfile(directory, 'Data/Cartesian/P02048.7');

CartesianRecon(pfile);