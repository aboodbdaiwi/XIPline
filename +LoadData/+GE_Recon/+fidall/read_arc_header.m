function header = read_arc_header(fname,verb,do_hdr_file,max_ntries)
%READ_ARC_HEADER  Read ScanArchive Header (w/o Ox)
%     header = read_arc_header(fname,verb,do_hdr_file,max_ntries)
%      fname   ScanArchive file name
%       verb   Verbose mode                                          (true)
%do_hdr_file   Always create a hdr from h5 file                     (false)
% max_ntries   Maximum #tries for h5read (each with 5s pause)          (30)
%     header   GE rawdata header structure
%
%  5/2023 Matteo Cencini
% 10/2024 Rolf Schulte
if (nargout<1), help(mfilename); return; end


%% input
if ~exist('verb','var'),        verb = []; end
if isempty(verb),               verb = true; end
if ~exist('do_hdr_file','var'), do_hdr_file = []; end
if isempty(do_hdr_file),        do_hdr_file = false; end
if ~exist('max_ntries','var'),  max_ntries = []; end
if isempty(max_ntries),         max_ntries = 30; end

if isempty(regexpi(fname,'\.h5$', 'once'))
    warning('file lacks suffix .h5: ''%s''',fname);
end
if ~exist(fname,'file'), error('file not found: ''%s''',fname); end
mver = get_mver;
if mver<9.1
    warning('strings/base64decode requires Matlab R2016b/V9.1; mver=%g',mver);
end


%% extracting header and writing to hdr file
fname_hdr = [fname(1:end-3) '.hdr'];
if ~exist(fname_hdr,'file') || do_hdr_file
    % load structure
    if verb, fprintf('Loading ScanArchive header from ''%s''\n',fname); end
    % assuming ieee-le byte ordering
    for ll=1:max_ntries
        % workaround for matlab bug
        try
            fdict = h5read(fname, '/DownloadData.xml');
            if verb, fprintf('h5read successful\n'); end
            break;
        catch ME
            fprintf('Error using h5read; trying again in 5sec (%d/%d)\n',...
                ll,max_ntries);
            fprintf('Error message: %s\n',ME.message);
        end
        if ll==max_ntries, error('%s',ME.message); end
        pause(5);
    end
    xmlstr = string(fdict{1}(:).');  % strings introduced R2016b (9.1)
    hstr = extractBetween(xmlstr, '<PoolHeaderBuffer', '</PoolHeaderBuffer>');
    if isempty(hstr), error('isempty(hstr); merge Matteo''s code back in'); end

    % extract header
    if verb, disp('Extracting header'); end
    % extract header bytes
    hstr = extractBetween(hstr, '<Data>', '</Data>');
    hstr = char(hstr);
    packed_data = matlab.net.base64decode(hstr);
    

    % write header to binary file
    if verb
        fprintf('Writing header to file ''%s''\n',fname_hdr);
    end
    [fid,errmsg] = fopen(fname_hdr,'w');
    if fid<0
        fprintf('Warning: fopen failed: ''%s''; trying tempdir\n',errmsg);
        [~,fn] = fileparts(fname_hdr);
        if isunix
            fname_hdr = [tempdir '/' fn '.hdr'];
        else
            fname_hdr = [tempdir '\' fn '.hdr'];
        end
        if verb, fprintf('Writing header to file ''%s''\n',fname_hdr); end
        [fid,errmsg] = fopen(fname_hdr,'w');
        if fid<0, error('Cannot open file ''%s'': ''%s''',fname_hdr,errmsg); end
    end
    fwrite(fid,packed_data);
    fclose(fid);
end


%% read in header with Matlab GE p-file reading routines
if verb
    fprintf('read_MR_headers(''%s'')\n',fname_hdr);
end
if ~exist(fname_hdr,'file')
    error('file not found: ''%s''',fname_hdr); 
end
header = read_MR_headers(fname_hdr);


end      % read_arc_header.m
