function [data,header] = read_arc(fname,verb,sort,get_mrcfg)
%READ_ARC Read ScanArchive data (w/o Ox)
%[data,header] = read_arc(fname,verb,sort,get_mrconfig)
%        fname   ScanArchive filename          (default=last h5 one in dir)
%         verb   Verbose mode (0-3)                                     (1)
%         sort   Sort data                                           (true)
%    get_mrcfg   Extract MRconfig and attach to header               (true)
%
%         data   Raw data (same dimensionality as p-files)
%                [nframes,frame_size,nphases,nechoes,nslices,nreceivers]
%       header   Header structure (read_MR_headers)
%
%  5/2023 Matteo Cencini
%  7/2024 Rolf Schulte
if (nargout<1), help(mfilename); return; end


%% input
if ~exist('verb','var'), verb = []; end
if isempty(verb),        verb = 1; end
if ~exist('sort','var'), sort = []; end
if isempty(sort),        sort = true; end
if ~exist('get_mrcfg','var'), get_mrcfg = []; end
if isempty(get_mrcfg),   get_mrcfg = true; end

if ~exist('fname','var'), fname = []; end
if isempty(fname)
    tmp = ls('*.h5');
    fname = tmp(end,:);
    if verb>0, fprintf('fname = ''%s''\n',fname); end
end
if isempty(regexpi(fname,'\.h5$', 'once'))
    warning('file lacks suffix .h5: ''%s''; adding',fname);
end
if ~exist(fname,'file'), error('file not found: ''%s''',fname); end


%% read ScanArchive header
if verb>0, fprintf('Reading ScanArchive header from ''%s''\n',fname); end
header = read_arc_header(fname,verb>1);
nn = [header.rdb_hdr.nframes,...
    header.rdb_hdr.frame_size,...
    header.rdb_hdr.nphases,...
    header.rdb_hdr.nechoes,...
    header.rdb_hdr.nslices,...
    (header.rdb_hdr.dab(2)-header.rdb_hdr.dab(1))+1];

if verb>0
    fprintf('psd_iname = ''%s''\n',header.image.psd_iname);
    fprintf(['nframes=%g  frame_size=%g  nphases=%g  nechoes=%g  '...
        'nslices=%g  nreceivers=%g\n'],nn);
    fprintf('data size = %g [MB] (single precision)\n',8*prod(nn)*1d-6);
end


%% read ScanArchive data
if verb>0, fprintf('Loading ScanArchive data from ''%s''\n',fname); end
h5data = loadh5(fname);


%% extract MRconfig from ScanArchive
if get_mrcfg
    if verb>0, fprintf('Extract MRconfig and attach to header\n'); end
    try
        header.mrconfig = sub_get_mrconfig(h5data);
    catch ME
        warning('Cannot extract mrconfig from data\n');
        fprintf('Error message: %s\n',ME.message);
    end
end


%% load controls
if verb>0, fprintf('Loading controls\n'); end
controls = sub_read_controls(h5data);
view_acq = controls.view_idx;
header.data_acq_tab.view_acq = view_acq(view_acq>0);


%% load frames
if verb>0, fprintf('Loading frames\n'); end
data = sub_read_frames(h5data, header, controls);
clear h5data


%% sort frames
if sort>0
	if verb>0, fprintf('Sorting data\n'); end
	[view_idx, slice_idx, echo_idx, op_idx, opcode] = ...
		sub_get_frame_ordering(controls);
	data = sub_sort_frames(data, view_idx, slice_idx, echo_idx, op_idx, ...
		opcode, header, verb);
end


%% checks
if isempty(data), warning('data is empty'); end
for ll=1:6
    if nn(ll) ~=size(data,ll)
        warning('nn(%d)(=%d) ~= size(data,%d)(=%d)',ll,nn(ll),ll,size(data,ll));
    end
end
if any(any(any(any(any(any(isnan(data)))))))
    warning('data contains NaN; replacing by 0');
    data(isnan(data)) = 0;
end
if any(any(any(any(any(any(isinf(data)))))))
    warning('data contains inf; replacing by 0');
    data(isinf(data)) = 0;
end
if isreal(data), warning('data is real'); end


end      % read_arc.m


%% ************************************************************************
% sub-functions
%
function controls = sub_read_controls(input)

% get scan name
archive_name = sub_get_archive_name(input);

% get controls
cbytes = sub_get_segments(input.Acquisition, strcat(archive_name, '_Block0_Segment'));

% get RawControl size
rawsize = double(cbytes(1));

% remove last 16 bytes and check remaining is a multiple of RawControl size
dummy = rem(length(cbytes), rawsize);
cbytes = cbytes(1:end-dummy);

% reshape to [rawcontrolsize, ncontrols]
cbytes = reshape(cbytes, [rawsize, length(cbytes) / rawsize]);
cbytes = cbytes(9:end, :);
cbytes = flip(cbytes, 1);

% parse controls
tmp = cbytes(end, :);
controls.opcode = typecast(tmp, 'uint8');
tmp = cbytes(end-8:end-7, :);
controls.slice_idx = typecast(tmp(:), 'int16');
tmp = cbytes(end-9, :);
controls.echo_idx = typecast(tmp, 'uint8');
tmp = cbytes(end-10, :);
controls.op_idx = typecast(tmp, 'uint8');
tmp = cbytes(end-12:end-11, :);
controls.view_idx = typecast(tmp(:), 'int16');

% find scan stop
tmp = flip(controls.opcode);
stop = find(tmp ~= 0, 1, 'first') - 1;
stop = length(tmp) - stop + 1;

% select
controls.opcode   = controls.opcode(1:stop);
controls.slice_idx = controls.slice_idx(1:stop);
controls.echo_idx = controls.echo_idx(1:stop);
controls.op_idx   = controls.op_idx(1:stop);
controls.view_idx = controls.view_idx(1:stop);

end      % sub_read_controls


%% ************************************************************************
function [view_idx, slice_idx, echo_idx, op_idx, opcode] = ...
    sub_get_frame_ordering(controls)
opcode = controls.opcode;
view_idx = controls.view_idx;
slice_idx = controls.slice_idx + 1;
echo_idx = controls.echo_idx + 1;
op_idx = controls.op_idx;

end      % sub_get_frame_ordering


%% ************************************************************************
% extract data frames
function frames = sub_read_frames(archive, h, ctrl)
%archive  archive structure
%      h  header structure
%   ctrl  control structure
% frames  data

% get scan name
archive_name = sub_get_archive_name(archive);

% get controls
cbytes = sub_get_segments(archive.Acquisition, strcat(archive_name, '_Block1_Segment'));
clear archive

% get info
ncoils = double(h.rdb_hdr.dab(2) - h.rdb_hdr.dab(1) + 1);
ndat = double(h.rdb_hdr.frame_size);
ptsize = double(h.rdb_hdr.point_size); % Either 2 (data stored in short int format, int16) or 4 (extended precision)

% select precision
switch ptsize
    case 2, precision = 'int16';
    case 4, precision = 'int32';
    otherwise, error('ptsize(=%g) ~= 2 or 4',ptsize);
end

% get header size
header_size = 48;

% get frame header starts
hstart = strfind(cbytes, cbytes(13:24)) - 12;    % slower
% hstart = regexp(char(cbytes), char(cbytes(13:24)))-12;  % bug
% test data: C:\BoxSync\data\202307xx_lise_qti_satrec\Exam2296anon\Series6

% get number of frames
nframes = length(hstart);

% get frame headers
fheaders = zeros(header_size, nframes, 'uint8');
for n = 1:nframes
    fheaders(:, n) = cbytes(hstart(n):hstart(n)+header_size-1);
end

% get sequence number
tmp = fheaders(end-7:end-4, :);
seq_idx = typecast(tmp(:), 'int32') + 1;

% sort
[~, order] = sort(seq_idx);
fstart = hstart(order) + header_size;

% select programmable packets
idx = logical((ctrl.opcode == 1) + (ctrl.opcode == 16));
idx = ctrl.opcode(idx).' == 1;

% actual selection
fstart = fstart(idx);
nframes = length(fstart);

% initialize frames
frame_size = 2 * ndat * ncoils * ptsize;
frames = zeros(nframes * frame_size, 1, 'uint8');

% get frames
for n = 1:nframes
    start = (n-1) * frame_size + 1;
    stop = n * frame_size;
    frames(start:stop) = cbytes(fstart(n):fstart(n)+frame_size-1);
end
clear cbytes

% cast
frames = typecast(frames, precision);

% get real and imaginary
% rframes = single(frames(1:2:end));
% iframes = single(frames(2:2:end));

% build complex
% frames = rframes + 1i * iframes;
frames = single(frames(1:2:end)) + 1i * single(frames(2:2:end));



% reshape
frames = reshape(frames, [ndat, ncoils, nframes]);

end      % sub_read_frames


%% ************************************************************************
function output = sub_sort_frames(input, view_idx, slice_idx, echo_idx, ...
    op_idx, opcode, h, verb)

% get psd name
psdname = deblank(h.image.psd_iname);
if verb, fprintf('psd_iname = ''%s''\n', psdname); end

% select frames with opcode = 1
keep = opcode == 1;

% actual selection
view_idx = view_idx(keep);
slice_idx = slice_idx(keep);
echo_idx = echo_idx(keep);
op_idx = op_idx(keep);
opcode = opcode(keep);

% get sizes from input
[ndat, ncoils, nframes] = size(input);
nviews = h.rdb_hdr.nframes;
nslices = h.rdb_hdr.nslices;
nechoes = h.rdb_hdr.nechoes;
nphases = h.rdb_hdr.nphases;
if nphases > 1
    error('nphases(=%g) > 1; not implemented', nphases);
end


% reshape if data consecutive
rshp = 0;     % default=full (slow) reordering in loop
if all(opcode==1) && all(op_idx==0)
    if all(diff(view_idx)==0)
        rshp = 1;
        if nslices > 1
            warning('view_idx consecutive but nslices(=%g) > 1',nslices);
            rshp = 0;
        end
        if nechoes > 1
            warning('view_idx consecutive but nechoes(=%g) > 1',nechos);
            rshp = 0;
        end
    else
        if (nechoes > 1) && (nslices == 1)
            rshp = 2;
            if ~all(slice_idx==1)
                warning('~all(slice_idx==1)');
                rshp = 0;
            end
            if any(diff(echo_idx)>1) || any(diff(echo_idx)<0)
                warning('echo_idx');
                rshp = 0;
            end
        end
    end
end
% reshape if data consecutive
switch rshp
    case 1
        if verb>0
            fprintf('fast reshaping\n');
        end
        output = permute(input,[3 1 4 5 6 2]);
    case 2         % nechoes>1
        if verb>0
            fprintf('fast reshaping with nechoes=%g\n',nechoes);
        end
        output = permute(reshape(input,[ndat ncoils nviews nechoes]),...
            [3 1 5 4 6 2]);
    otherwise
        if verb>0
            fprintf('full (slow) reordering in loop\n');
        end
        % initialize data
        output = complex(zeros(nviews, ndat, nphases, nechoes, nslices, ...
            ncoils, 'single'));

        % prepare datasets
        if verb > 1
            fprintf('Opcode SliceNum EchoNum Operation ViewNum\n');
        end
        % input = permute(input,[3 1 4 5 6 2]);
        for ll = 1:nframes
            if opcode(ll) == 1
                % get operation, echo, slice and view indexes
                op = op_idx(ll);
                iecho = echo_idx(ll);
                islice = slice_idx(ll);
                iview = view_idx(ll);

                % print
                if verb > 1
                    fprintf('%g %g %g %g %g %g\n', ...
                        opcode(ll), islice, iecho, op, iview, ll);
                end

                % apply operation
                if iview > 0
                    dtmp = reshape(input(:, :, ll), [1, ndat, 1, 1, 1, ncoils]);
                    switch op
                        case 0       % DABSTORE
                            output(iview, :, 1, iecho, islice, :) = dtmp;
                        case 1       % DABADD
                            output(iview, :, 1, iecho, islice, :) = ...
                                output(iview, :, 1, iecho, islice, :) + dtmp;
                        case 2       % DABSUBCNTS
                            output(iview, :, 1, iecho, islice, :) = ...
                                -output(iview, :, 1, iecho, islice, :) + dtmp;
                        case 3       % DABSUBXCVR
                            output(iview, :, 1, iecho, islice, :) = ...
                                output(iview, :, 1, iecho, islice, :) - dtmp;
                        otherwise
                            warning('ll=%g: unknown control.operation(=%g)',ll,op);
                    end

                end
            end
            if verb > 1
                fprintf('\n');
            end
        end
end
clear input

% adjust slice ordering
if nslices > 1
    do_slice_order = true;
else
    do_slice_order = false;
end
if ~isempty(regexpi(psdname,'3drad', 'once')) || ...
        ~isempty(regexpi(psdname,'burzte', 'once'))
    do_slice_order = false;
end

if do_slice_order
    slice_order = zeros(nslices, 1);
    for ll = 1:nslices
        slice_order(ll) = h.data_acq_tab.slice_in_pass(ll);
    end
    if verb > 0
        fprintf('slice_order = %s\n',num2str(slice_order(:).'));
    end
    if any(slice_order==0) || any(diff(sort(slice_order))~=1)
        warning('corrupt slice_order -> skipping');
    else
        output = output(:, :, :, :, slice_order, :);
    end
end

end      % sub_sort_frames


%% ************************************************************************
function aname = sub_get_archive_name(input)
mstr = input.Acquisition.StorageMetaData_0x2E_xml;

% get name string
if size(mstr, 1) > size(mstr, 2)
    mstr = mstr.';
end
mstr = string(mstr);

% get actual name
aname = extractBetween(mstr, '<Name>', '</Name>');

% find actual name
nameidx = [];
for n = 1:length(aname)
    nameidx = [nameidx, length(aname{n})];
end
[~, nameidx] = min(nameidx);
aname = char(aname{min(nameidx)});

end      % sub_get_archive_name


%% ************************************************************************
function seg = sub_get_segments(input, str)

mstr = input.StorageMetaData_0x2E_xml;

% get name string
if size(mstr, 1) > size(mstr, 2)
    mstr = mstr.';
end
mstr = string(mstr);

% get segments name
segnames = extractBetween(mstr, '<Name>', '</Name>');

% get desired indices
segidx = [];
for ll = 1:length(segnames)
    segidx = [segidx, contains(segnames{ll}, str)];
end
segidx = logical(segidx);

segnames = sort(segnames(segidx));

% get segment bytes
seg = [];
for ll = 1:length(segnames)
    seg = [seg, input.(segnames{ll})];
end

end      % sub_get_segments


%% extract MRconfig from data
function mrconfig = sub_get_mrconfig(data)

str = data.w.config.MRconfig_0x2E_cfg.';

iinl = [0 regexpi(str,'\n')];          % newline index list

for ll=1:length(iinl)-1
    ii = (iinl(ll)+1):(iinl(ll+1)-1);
    if ~isempty(ii)
        line = str(ii);
        iieq = regexpi(line,'=','once');
        if isempty(iieq)
            warning('isempty(iieq): line=%s',line);
        else
            name = strtrim(line(1:(iieq-1)));
            name = strrep(name,'.','_');
            val = strtrim(line((iieq+1):end));
            val = strrep(val,'"','');
            if ~isnan(str2double(val))
                val = str2double(val);
            end
            mrconfig.(name) = val;
        end
    end
end

end      % sub_get_mrconfig
