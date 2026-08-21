% 3DASL
directory = fileparts(mfilename('fullpath'));
scanArchiveFile = fullfile(directory, 'Data/3DASL/ScanArchive.h5');
ASL3dScanArchiveRecon(scanArchiveFile);