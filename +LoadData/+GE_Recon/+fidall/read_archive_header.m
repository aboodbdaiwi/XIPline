function [h,archive]=read_archive_header(fname)
%READ_ARCHIVE Read hdf5 MR raw data file header using Orchestra tools
%[h,archive] = read_archive_header(fname)
%      fname   filename              (default=last h5 one in dir)
%          h   header structure; compatible to p-files
%    archive   original archive structure
%
% 11/2022 Rolf Schulte
%See also GERecon, archive2header.
if (nargout<1), help(mfilename); return; end


%% input parameters
verb = true;
if ~exist('fname','var'), fname = []; end
if isempty(fname)
    xx = dir('*.h5');
    if length(xx)>1
        warning('multiple files found; chosing last');
    end
    if isempty(xx)
        error('file ''%s'' not found',fname);
    end
    if isunix
        dsep = '/';
    else 
        dsep = '\';
    end
    fname = [xx(end).folder dsep xx(end).name];
    if ~exist(fname,'file'), error('file %s not found',fname); end
    if verb>0, fprintf('Loading header fname=\n\t%s\n',fname); end
end
if ~exist(fname,'file'),  error('file ''%s'' not found',fname); end


%% open archive and load header
archive = GERecon('Archive.Load', fname);
h = [];


%% checks
if ~strcmpi(archive.ScanType,'scan')
    fprintf('ScanType = ''%s'' -> skipping\n',archive.ScanType);
    return
end
if ~isfield(archive.DownloadData,'rdb_hdr_rec')
    GERecon('Archive.Close', archive);
    fprintf('ScanType = ''%s''\n',archive.ScanType);
    warning('No rdb_hdr_rec field');
    return
end


%% Close the archive to release it from memory
try
    GERecon('Archive.Close', archive);
catch ME
    fprintf('Error using GERecon\n');
    fprintf(ME.message);
    return
end


%% convert archive structure to p-file header structure
h = archive2header(archive);


end      % read_archive_header.m
