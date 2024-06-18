function [d,h,archive]=read_archive(fname,verb,do_single,no_add,do_slice_order)
%READ_ARCHIVE Read hdf5 MR raw data file using Orchestra tools
%[d,h,archive] = read_archive(fname,verb,do_single,no_add,do_slice_order)
%     fname  filename              (default=last h5 one in dir)
%      verb  verbose mode          (0=off;1=simple,2=more,3=all)
% do_single  d in single precision (default=false; true for 3drad/burzte)
%    no_add  force control.operation=0 (default=false)
%do_slice_order  perform slice ordering 
%                (default=true (false for 3drad/burzte))
%
%         d  raw data matrix (same dimensionality as p-files)
%            [nframes,frame_size,nphases,nechoes,nslices,nreceivers]
%         h  header structure; compatible to p-files
% archive  original archive structure
%
% 10/2020  Rolf Schulte
%See also GERecon, archive2header.
if (nargout<1), help(mfilename); return; end

if ~exist('verb','var'),  verb = []; end
if isempty(verb),         verb = 1; end
if ~exist('fname','var'), fname = []; end
if isempty(fname)
    tmp=ls('*.h5');
    fname = tmp(end,:);
    if verb>0, fprintf('fname = ''%s''\n',fname); end
end
if ~exist(fname,'file'),  error('file ''%s'' not found',fname); end
if ~exist('do_single','var'), do_single = []; end
if ~exist('no_add','var'), no_add = []; end
if isempty(no_add),        no_add = false; end
if ~exist('do_slice_order','var'), do_slice_order = []; end


%% open archive and load header
archive = LoadData.ismrmrd.GE.GERecon('Archive.Load', fname);

if ~isfield(archive.DownloadData,'rdb_hdr_rec')
    GERecon('Archive.Close', archive);
    fprintf('ScanType = ''%s''\n',archive.ScanType);
    error('No rdb_hdr_rec field');
end

psd_iname = archive.DownloadData.rdb_hdr_image.psd_iname;
if verb>0, fprintf('psd_iname = ''%s''\n',psd_iname); end

if isempty(do_single)
    if ~isempty(regexpi(psd_iname,'3drad')) || ...
            ~isempty(regexpi(psd_iname,'burzte'))
        do_single = true;
    else
        do_single = false;
    end
end


%% Scan parameters
stop = archive.DownloadData.rdb_hdr_rec.rdb_hdr_dab.stop_rcv;
start = archive.DownloadData.rdb_hdr_rec.rdb_hdr_dab.start_rcv;
nreceivers = stop - start + 1;
nframes = archive.DownloadData.rdb_hdr_rec.rdb_hdr_nframes;
frame_size = archive.DownloadData.rdb_hdr_rec.rdb_hdr_frame_size;
nphases = archive.DownloadData.rdb_hdr_rec.rdb_hdr_nphases;
nechoes = archive.DownloadData.rdb_hdr_rec.rdb_hdr_nechoes;
nslices = archive.DownloadData.rdb_hdr_rec.rdb_hdr_nslices;
if nphases>1
    GERecon('Archive.Close', archive);
    error('nphases(=%g)>1; not implemented',nphases);
end
view_offset = [];

%% pre-allocate memory for data
if do_single
    d = complex(zeros(nframes,frame_size,nphases,nechoes,nslices,nreceivers,'single'));
else
    d = complex(zeros(nframes,frame_size,nphases,nechoes,nslices,nreceivers));
end
arr = zeros(1,7);


%% Loop through each control, sorting frames if applicable
if verb>1
    fprintf('nframes=%g  frame_size=%g  nphases=%g  nechoes=%g  nslices=%g  nreceivers=%g\n',...
        nframes,frame_size,nphases,nechoes,nslices,nreceivers);
    fprintf('l opcode echoTrainIndex sliceNum echoNum operation viewNum FrameCount\n');
end
for l = 1:archive.ControlCount
    control = GERecon('Archive.Next', archive);
    if verb>1, fprintf('%g  %g',l,control.opcode); end
    
    % Sort only programmable packets in range
    if (control.opcode == 1)
        if verb>1
            fprintf('  %g %g %g %g %g %g',...
                control.echoTrainIndex,control.sliceNum,...
                control.echoNum,control.operation,control.viewNum,control.FrameCount);
        end
        if verb>0
            tmp = [control.opcode control.echoTrainIndex ...
                control.sliceNum control.echoNum control.operation ...
                control.viewNum control.FrameCount];
            arr = max([arr ; tmp],[],1);
        end
        % if ~exist('view_offset','var')
        if isempty(view_offset)
            switch control.viewNum
                case 0, view_offset = 1;
                case {1,3}, view_offset = 0;
                otherwise
                    % GERecon('Archive.Close', archive);
                    warning('control.viewNum(=%g) ~= 0 or 1',control.viewNum);
                    view_offset = 0;
            end
            if verb>0, fprintf('\nview_offset=%g\n',view_offset); end
        end
        
        % d(control.viewNum+view_offset,:,1,control.echoNum+1,control.sliceNum+1,:) = ...
        %   control.Frame.Data;
        if do_single
            % data = single(reshape(control.Frame.Data,[1 frame_size 1 1 1 nreceivers]));
            data = single(reshape(control.Data,[1 frame_size 1 1 1 nreceivers]));
        else
            % data = double(reshape(control.Frame.Data,[1 frame_size 1 1 1 nreceivers]));
            data = double(reshape(control.Data,[1 frame_size 1 1 1 nreceivers]));
        end
        l1 = control.viewNum+view_offset;
        l4 = control.echoNum+1;
        l5 = control.sliceNum+1;
        % loaddab -> dabop -> control.operation; according to epic.h:
        % #define DABSTORE 0   /* DAB operation to store xcvr data */
        % #define DABADD 1     /* DAB operation to accumulate data */
        % #define DABSUBCNTS 2 /* DAB operation to subtract contents from xcvr data */
        % #define DABSUBXCVR 3 /* DAB operation to subtract xcvr data from contents */
        dabop = control.operation;
        if no_add, dabop = 0; end
        switch dabop
            case 0, d(l1,:,1,l4,l5,:) = data;                     % DABSTORE
            case 1, d(l1,:,1,l4,l5,:) = d(l1,:,1,l4,l5,:) + data; % DABADD
            case 2, d(l1,:,1,l4,l5,:) = data - d(l1,:,1,l4,l5,:); % DABSUBCNTS
            case 3, d(l1,:,1,l4,l5,:) = d(l1,:,1,l4,l5,:) - data; % DABSUBXCVR
            otherwise
                warning('l=%g: unknown control.operation(=%g)',...
                    l,control.operation);
        end
    end
    if verb>1, fprintf('\n'); end
    
    if (control.opcode==0) % end of pass and/or scan
        if l~=archive.ControlCount
            warning('control.opcode==0: l(=%g)~=archive.ControlCount(=%g)',...
                l,archive.ControlCount);
        end
    end
end


if verb>0
    fprintf('l opcode echoTrainIndex sliceNum echoNum operation viewNum FrameCount\n');
    fprintf('max= %g  %g  %g  %g  %g  %g  %g  %g\n',l,arr);
    fprintf('nframes=%g  frame_size=%g  nechoes=%g  nslices=%g  nreceivers=%g\n',...
        nframes,frame_size,nechoes,nslices,nreceivers);
end

%% adjust slice ordering
if isempty(do_slice_order)
    if nslices>1, do_slice_order = true;
    else,         do_slice_order = false; end
    if ~isempty(regexpi(psd_iname,'3drad')) || ...
            ~isempty(regexpi(psd_iname,'burzte'))
        do_slice_order = false;
    end
end

if do_slice_order
    slice_order = zeros(1,nslices);
    for l=1:nslices, slice_order(1,l) = archive.DownloadData.rdb_hdr_data_acq_tab(l).slice_in_pass; end
    if verb>0, fprintf('Slice_order = '); disp(slice_order); end
    if any(slice_order==0) || any(diff(sort(slice_order))~=1)
        warning('corrupt slice_order -> skipping');
    else
        d = d(:,:,:,:,slice_order,:);
    end
end


%% convert archive structure to p-file header structure
if nargout>1
    h = LoadData.ismrmrd.GE.Functions.archive2header(archive,verb>2);
    for lslice=1:nslices
        arcInfo = GERecon('Archive.Info', archive, lslice);
        h.orientation{lslice} = arcInfo.Orientation;
        h.corners{lslice} = arcInfo.Corners;
    end
end


%% Close the archive to release it from memory
% GERecon('Archive.Close', archive);
% if closed: error in write_dicom_ox.m

end   % main function read_archive.m
