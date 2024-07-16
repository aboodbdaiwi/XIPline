function [h,archive]=read_archive_header(fname)
%READ_ARCHIVE Read hdf5 MR raw data file header using Orchestra tools
%[h,archive] = read_archive_header(fname)
%      fname   filename              (default=last h5 one in dir)
%          h   header structure; compatible to p-files
%    archive   original archive structure
%
% 1/2021  Rolf Schulte
%See also GERecon, archive2header.
if (nargout<1), help(mfilename); return; end


%% input parameters
if ~exist('fname','var'), fname = []; end
if isempty(fname)
    tmp=ls('*.h5');
    fname = tmp(end,:);
    if verb>0, fprintf('fname = ''%s''\n',fname); end
end
if ~exist(fname,'file'),  error('file ''%s'' not found',fname); end


%% open archive and load header
archive = GERecon('Archive.Load', fname);

if ~isfield(archive.DownloadData,'rdb_hdr_rec')
    GERecon('Archive.Close', archive);
    fprintf('ScanType = ''%s''\n',archive.ScanType);
    error('No rdb_hdr_rec field');
end


%% Close the archive to release it from memory
GERecon('Archive.Close', archive);


%% convert archive structure to p-file header structure
h = archive2header(archive);


end      % read_archive_header.m
