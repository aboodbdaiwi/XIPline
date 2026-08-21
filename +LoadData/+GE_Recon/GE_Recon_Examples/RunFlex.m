% Call FlexRecon script with the Pfile and kacq file

directory = fileparts(mfilename('fullpath'));
pfile = fullfile(directory, 'Data/Flex/P11776.7');
kacq = fullfile(directory, 'Data/Flex/kacq_yz.txt.906151242');

FlexRecon(pfile, kacq);