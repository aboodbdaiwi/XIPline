function [raw, header, ec] = read_MR_rawdata(in_name, save_baseline, phaselist, echolist, slicelist, rcvrlist)
%read_MR_rawdata - Read MR rawdata files (P-files)
%
% The tool suppports rawdata files with multiple temporal phases, echos,
% slices, and receivers.
%
%   [A, HEADER] = read_MR_rawdata(FILENAME);    returns the
%   MR rawdata in FILENAME into A, and returns the headers
%   into HEADER. Array A may be up to 6-D, with dimensions of
%   (ky, kx, phase, echo, slice, rcvr).
%
%   A = read_MR_rawdata(FILENAME);    returns the rawdata into A;
%   since the headers are not returned, execution is much faster.
%
%   A = read_MR_rawdata(FILENAME, BASELINE_FLAG, PHASE, ECHO, SLICE, RCVR);
%   returns the rawdata into A, using the chosen echo, slice and rcvr.
%   PHASE, ECHO, SLICE, and RCVR can be a single value or a range of values.
%   The BASELINE_FLAG must be either 'sb' for "save baseline" or 'db' for
%   "discard baseline".
%
% If the rawfile contains only one image worth of rawdata, or if only one
% phase, echo, slice, and receiver is requested, the returned matrix A will be 2-D.
% The index order of the returned matrix will be:
% [YRES XRES PHASES ECHOES SLICES RECEIVERS]
%

% Copyright (c) 2019 by General Electric Company. All rights reserved.

ec=0;
non_standard_pfile=0;
% Determine if the input parameter IN_NAME points to a directory of a file
fprintf('\nREAD_MR_RAWDATA:\n')
switch exist(in_name)
    case 7        % it's a DIRECTORY -- open gui to get filename
        [ filename, dirname ] = uigetfile( fullfile(in_name, '*.7'), 'Select a P-file to open' );
        if isequal(filename,0) || isequal(dirname,0)
            errordlg('No file selected. No data returned.','ERROR!','modal');
            raw = [];
            ec=1;
            return;
        end
    case 2        % it's a FILE
        [ dirname, filename, fileext] = fileparts(in_name);
        filename = [filename fileext];
        [~, num_char]=size(filename);
        % replace strmatch with startsWith when do not have to support < R2016b
        P_check = strmatch('P',filename);
        ext_check = strcmp(fileext,'.7');
        if (num_char ~= 8 || P_check ~= 1 || ext_check ~= true)
            non_standard_pfile=1;
        end
    otherwise
        errordlg('Parameter is not a file or directory. No data returned.','ERROR!','modal');
        raw = [];
        ec=1;
        return;
end

fullfilename = fullfile( dirname, filename);
header = read_MR_headers( fullfilename, 'all' );

% Get endian from header structure
endianID = header.endian;
% Pull useful info from rdb_hdr
nfphases   = header.image.fphase;
if nfphases>6d4   % RFS fix for broken header
    warning('header.image.fphase(=%d)>6d4: using header.rdb_hdr.nphases(=%d) instead',...
        header.image.fphase,header.rdb_hdr.nphases);
    nfphases = header.rdb_hdr.nphases;
end
psdname    = deblank(header.image.psd_iname);
npasses    = header.rdb_hdr.npasses;
da_xres    = header.rdb_hdr.da_xres;
da_yres    = header.rdb_hdr.da_yres;
nslices    = header.rdb_hdr.nslices;
nechoes    = header.rdb_hdr.nechoes;
point_size = header.rdb_hdr.point_size;
if isfield(header.rdb_hdr, 'sspsave')    % RFS bug fix
    sspsave    = header.rdb_hdr.sspsave;
else
    sspsave = 0;
end
if isfield(header.rdb_hdr, 'dacq_data_type')
    dacq_data_type = header.rdb_hdr.dacq_data_type;
else
    dacq_data_type = 0;
end
receivers = (header.rdb_hdr.dab(2)-header.rdb_hdr.dab(1))+1;

  % RFS bug fix:
  if ~isfield(header.rdb_hdr,'raw_pass_size')
      warning('RFS bug fix: nslices_per_pass set to %g',...
          sum( header.data_acq_tab.pass_number == 0 ));
      nslices_per_pass = sum( header.data_acq_tab.pass_number == 0 );  % Number of slices in each pass
  else
      % original
      nslices_per_pass = header.rdb_hdr.raw_pass_size / (receivers * nechoes * da_yres * da_xres * point_size * 2);
  end

switch dacq_data_type
    case 0
        switch point_size
            case 2
                raw_data_type = 'int16';   % 16 bit data
            case 4
                raw_data_type = 'int32';   % 32 bit data
            otherwise
                fprintf('ERROR: Invalid point size (%d)\n', point_size);
                ec = 2;
                return;
        end
    case 1
        raw_data_type = 'int16';    % 16 bit data
    case 2
        raw_data_type = 'int32';    % 32 bit data
    case 3
        raw_data_type = 'float';    % floating data
    case 4
        raw_data_type = 'double';   % double data
    otherwise
        fprintf('ERROR: Invalid DACQ data type (%d)\n', dacq_data_type);
        ec = 2;
        return;
end


   % +++ RFS +++ bug fix
  if header.rdb_hdr.zcsi>1
      warning('RFS bug fix for 3D CSI: zcsi=%g',header.rdb_hdr.zcsi);
      fprintf('  Setting nslices(=%g) & nslices_per_pass(=%g) to 1',...
          nslices,nslices_per_pass);
      nslices = 1;
      nslices_per_pass = 1;
  end
  if ~isempty(regexpi(header.image.psd_iname,'3drad')) || ...
          ~isempty(regexpi(header.image.psd_iname,'silen'))|| ...
          ~isempty(regexpi(header.image.psd_iname,'burzte'))
      fprintf('psd_iname=%s\n',header.image.psd_iname);
      warning('RFS bug fix for 3dradial: reading data manually');
      if exist('save_baseline','var')
          if ~strcmp(save_baseline,'db'), error('save_baseline not ''db'''); end
      end
      fid = fopen(fullfile( dirname, filename),'r',header.endian);
      status = fseek(fid, header.total_length, 'bof');
      if status==0, warning('fseek'); end
      raw = fread(fid,inf,raw_data_type); % rawdata=raw;
      fclose(fid);
      % RAWDATA sorting
      raw = reshape(raw(1:2:end)+1i*raw(2:2:end),length(raw)/receivers/2,receivers);
      % raw = reshape(raw(1:da_xres*da_yres*nslices,:),da_xres,da_yres,nslices,receivers);
      % raw = reshape(permute(raw(:,2:da_yres,:,:),[2 1 3 4]),...
      %     (da_yres-1),da_xres,1,1,nslices,receivers);
      raw = reshape(raw(1:da_xres*da_yres*nslices*nechoes,:),da_xres,da_yres,nslices,receivers,nechoes);
      
      raw = reshape(permute(raw(:,2:da_yres,:,:,:),[2 1 3 4 5]),...
          (da_yres-1),da_xres,1,nechoes,nslices,receivers);
      return
  end
  % --- RFS ---


fprintf('nfphases = %d psdname = %s npasses = %d da_xres = %d da_yres = %d nslices = %d nechoes = %d point_size = %d receivers= %d nslices_per_pass = %d\n', ...
    nfphases, psdname, npasses, da_xres, da_yres, nslices, nechoes, point_size, receivers, nslices_per_pass);

is_grafimage_psd = 0;
if (strcmp(psdname, 'GRAFIMAGE'))
    is_grafimage_psd = 1;
end

if (non_standard_pfile == 1 && npasses ~= 1 && is_grafimage_psd == 0)
    errordlg('A non-standard P file name must contain only 1 pass.','ERROR!','modal');
    raw = [];
    ec=2;
    return;
end

if (nfphases == 0)
    nphases = 1;
else
    nphases = nfphases;
end

% Unless passed at least 6 input parameters, read ALL the data in the rawfile
if (nargin == 1)
    phaselist = 1:nphases;
    echolist  = 1:nechoes;
    slicelist = 1:nslices;
    rcvrlist  = 1:receivers;
    save_baseline = 'db';
elseif (nargin == 2)
    phaselist = 1:nphases;
    echolist  = 1:nechoes;
    slicelist = 1:nslices;
    rcvrlist  = 1:receivers;
elseif (nargin ~= 6)
    errordlg('You must call read_MR_radata with either 1 or 6 arguments.','ERROR!','modal');
    raw = [];
    ec=2;
    return;
end

% Compute (byte) size of frames, echoes and slices
% keep definition, but do not execute : elements  = da_xres*2;
frame_sz  = 2 * da_xres * point_size;     % Size of one frame of rawdata (one Kx line)
echo_sz   = frame_sz * da_yres;           % Size of one image (one Kx * Ky matrix)
slice_sz  = echo_sz  * nechoes;           % Size of nechoes images (one raw matrix * number of echoes)
mslice_sz = slice_sz * receivers;	      % (one slice * number of receivers)

fprintf('length slicelist = %d nphases = %d nslices_per_pass = %d npasses = %d \n', ...
    length(slicelist), nphases, nslices_per_pass, npasses);
% Determine which of the 3 types of temporal phase raw data it is
temporal_type='na';
if (nphases > 1)                         % It's a multi temporal phase scan
    if (nslices_per_pass == nphases)    % It's a Sequential Multi-Phase acq. Each raw file contains 1 slice all phases.
        temporal_type='seq';
        if (nargin == 1)
            slicelist=1:npasses;
        end
    elseif (npasses == nphases)         % It's an Interleaved Multi-Phase acq. Each raw file contains all slices at one temporal phase.
        temporal_type='int';
        if (nargin == 1)
            slicelist=1:nslices_per_pass;
        end
    elseif strcmp(psdname,'EPI')      % It's a fMRI scan. Only one raw file.
        temporal_type='fMRI';
    end
    
    if (nechoes > 1)
        errordlg('read_MR_rawdata will not support multi-echo multi-phase data.','ERROR!','modal');
        raw=[];
        ec=3;
        return;
    end
end

if (is_grafimage_psd == 1)
    temporal_type = 'grafimage';
end

fprintf('length slicelist = %d nphases = %d nslices_per_pass = %d npasses = %d temporal_type = %s\n', ...
    length(slicelist), nphases, nslices_per_pass, npasses, temporal_type);

% Compute offset in bytes to start of raw_data, and the act_yres to use.
if strcmp(save_baseline,'sb')
    baseline_offset=0;
    act_yres=da_yres;
elseif strcmp(save_baseline,'db')
    baseline_offset   = da_xres * 2 * point_size;
    act_yres=da_yres-1;
else
    errordlg('You must specify sb or db for baseline argument.','ERROR!','modal');
    raw=[];
    ec=4;
    return;
end

% use within for-loops instead of re-calculating each time
NumberDataElementsToFread = 2 * da_xres * act_yres ;

if (nphases > 1)
    if strcmp(temporal_type,'seq')
        total_bytes = header.total_length + echo_sz * nphases * receivers;
        raw_file_skip=nphases;
    elseif strcmp(temporal_type,'int')
        total_bytes = header.total_length + echo_sz * nslices_per_pass * receivers;
        raw_file_skip=nslices_per_pass;
    elseif strcmp(temporal_type,'fMRI')
        total_bytes = header.total_length + echo_sz * nphases * nslices_per_pass * receivers;
        raw_file_skip=1;
    else
        errordlg('The raw data set is multiphase but is not a recognized temporal type.','ERROR!','modal');
        raw=[];
        ec=5;
        return;
    end
else
    total_bytes = header.total_length + nslices_per_pass * mslice_sz + sspsave;
    raw_file_skip=1;
end

if (strcmp(temporal_type, 'grafimage'))
    total_bytes = header.total_length + npasses * nslices_per_pass * mslice_sz;
    fprintf('total_length = %d npasses = %d nslices_per_pass = %d mslice_sz = %d\n', ...
        header.total_length, npasses, nslices_per_pass, mslice_sz);
    npasses = 1;  % all P files are concatenated into one in grafimage data
    raw_file_skip=1;
end

% Check for complete raw data set: all files present and of the expected length
fprintf('\nRawfile Dataset information:\n');
fprintf('  File format:            %10s\n', header.format );
fprintf('  Endian ID:              %10s\n', header.endian );
fprintf('  Header size, bytes:     %10d\n', header.total_length );
fprintf('  Slice size, bytes:      %10d\n', mslice_sz );
fprintf('  Number of Time Phases:  %10d\n', nphases);
fprintf('  Expected Pfile, bytes:  %10d\n', total_bytes );
fprintf('Check rawfiles Status:\n');
file_inc=raw_file_skip;
for filei = 1:npasses
    if (non_standard_pfile == 0)
        filename = sprintf('P%05d.7', header.rdb_hdr.run_int+file_inc-raw_file_skip);
        file_inc=file_inc+raw_file_skip;
        fidi = fopen( fullfile( dirname, filename), 'r', endianID );
    else
        fidi = fopen( fullfilename, 'r', endianID);
    end
    if fidi < 3
        status_str = 'Could not be opened';
    else
        fseek( fidi, 0, 'eof' );
        filesize = ftell( fidi );
        fclose( fidi );
        if filesize~=total_bytes
            status_str = sprintf('Wrong rawfile size: %d bytes (actual) \n', filesize);
        else
            status_str = 'OK';
        end
    end
    fprintf( '  Pass %2d rawfile:%12s   %s\n', filei, filename, status_str);
end

if strcmp(temporal_type,'na')
    fprintf('  Echo Slice Rcvr Pass Filename     Offset   Status\n');
    fprintf('  -------------------------------------------------\n');
    raw = zeros( act_yres, da_xres, max(phaselist),max(echolist), max(slicelist), max(rcvrlist) );
    for echonum = min(echolist):max(echolist)
        echo = echolist(echonum - min(echolist) + 1);
        echo_offset = echo_sz * (echo-1);
        for slicenum = min(slicelist):max(slicelist)
            slice_requested = slicelist(slicenum - min(slicelist) + 1);
            pass = header.data_acq_tab.pass_number(slice_requested);
            slice_offset = slice_sz * (header.data_acq_tab.slice_in_pass(slice_requested)-1);
            for rcvrnum = min(rcvrlist):max(rcvrlist)
                rcvr = rcvrlist(rcvrnum - min(rcvrlist) + 1);
                receiver_offset = nslices_per_pass * slice_sz * (rcvr-1);
                file_offset = header.total_length + receiver_offset + slice_offset + echo_offset + baseline_offset;
                if (non_standard_pfile == 0)
                    filename = sprintf('P%05d.7', header.rdb_hdr.run_int+header.data_acq_tab.pass_number(slice_requested));
                    fid = fopen( fullfile( dirname, filename), 'r', endianID);
                else
                    fid = fopen( fullfilename, 'r', endianID);
                end
                % execute below even if file does not exist, to have
                % workflow of fprintf information to user, like in other types, especially status
                if fid<3
                    status_str = 'File not found';
                else
                    status = fseek(fid, file_offset, 'bof');
                    if status==-1
                        status_str = 'Offset outside of file';
                    else
                        [raw_data, count] = fread( fid, NumberDataElementsToFread, raw_data_type);
                        if count~=NumberDataElementsToFread
                            status_str = 'Could not read sufficient data';
                        else
                            status_str = 'OK';
                            raw_data = complex( raw_data(1:2:count), raw_data(2:2:count) );
                            raw_data = reshape( raw_data, da_xres, act_yres );
                            raw(:,:,1,echonum, slicenum, rcvrnum) = transpose( raw_data );
                        end
                    end
                end % if fid<3
                
                if fid>=3
                    fclose(fid);
                end
                
% RFS                fprintf(' %4d %4d %4d %4d %10s %10.0f   %s\n', ...
% RFS                    echo, slice_requested, rcvr, pass+1, filename, file_offset, status_str);
                 
            end
        end
    end
    
elseif strcmp(temporal_type,'grafimage')
    
    % filename does not change within for loops, so fopen file once outside of for loops;
    % fseek originally and always done from bof, so is not iteration-dependent
    
    if (non_standard_pfile == 0)
        USEfilename = sprintf('P%05d.7', header.rdb_hdr.run_int);
        USEfilename = fullfile(dirname, USEfilename);
    else
        USEfilename = fullfilename;
    end
    
    if exist( USEfilename , 'file' )~=2
        errordlg('Specified file not found. No data returned.','ERROR!','modal');
        raw = [];
        ec=1;
        return;
    end
    
    fid = fopen( USEfilename, 'r', endianID);
    
    if fid<3
        errordlg('File found but cannot be opened. No data returned.','ERROR!','modal');
        raw = [];
        ec=1;
        return;
    end
    
    % only fprintf if file exists and continue execution :
    fprintf('act_yres = %d da_xres = %d phaselist = %d echolist = %d slicelist = %d rcvrlist = %d\n', ...
        act_yres, da_xres, length(phaselist),length(echolist), length(slicelist), length(rcvrlist) );
    
    fprintf('  Echo Slice Rcvr Pass Filename     Offset   Status\n');
    fprintf('  -------------------------------------------------\n');
    
    raw = zeros( act_yres, da_xres, length(phaselist),length(echolist), length(slicelist), length(rcvrlist) );
    for echonum = 1:length(echolist)
        echo = echolist(echonum);
        echo_offset = echo_sz * (echo-1);
        for slicenum = 1:length(slicelist)
            slice_requested = slicelist(slicenum);
            pass = 0;
            slice_offset = slice_sz * (header.data_acq_tab.slice_in_pass(slice_requested)-1+floor((slice_requested-1)/nslices_per_pass)*nslices_per_pass);
            
            for rcvrnum = 1:length(rcvrlist)
                rcvr = rcvrlist(rcvrnum);
                receiver_offset = nslices_per_pass * slice_sz * (rcvr-1);
                file_offset = header.total_length + receiver_offset + slice_offset + echo_offset + baseline_offset;
                
                status = fseek(fid, file_offset, 'bof');
                if status==-1
                    status_str = 'Offset outside of file';
                else
                    [raw_data, count] = fread(fid, NumberDataElementsToFread, raw_data_type);
                    if count~=NumberDataElementsToFread
                        status_str = 'Could not read sufficient data';
                    else
                        status_str = 'OK';
                        raw_data = complex( raw_data(1:2:count), raw_data(2:2:count) );
                        raw_data = reshape( raw_data, da_xres, act_yres );
                        raw(:,:,1,echonum, slicenum, rcvrnum) = transpose( raw_data );
                    end
                end
                
                fprintf(' %4d %4d %4d %4d %10s %10.0f   %s\n', ...
                    echo, slice_requested, rcvr, pass+1, filename, file_offset, status_str);
            end
        end
    end
    fclose(fid);
    
elseif strcmp(temporal_type,'seq')
    fprintf('  Slice Phase Rcvr Filename     Offset   Status\n');
    fprintf('  -----------------------------------------------------------\n');
    raw = zeros( act_yres, da_xres, max(phaselist), max(echolist), max(slicelist), max(rcvrlist) );
    rcvr=min(rcvrlist);
    for slicenum = min(slicelist):max(slicelist)
        slice_requested = slicelist(slicenum - min(slicelist) + 1);
        for phasenum = min(phaselist):max(phaselist)
            phase_requested = phaselist(phasenum - min(phaselist) + 1);
            phase_offset = echo_sz * rcvr * (phasenum -1);
            for rcvrnum = min(rcvrlist):max(rcvrlist)
                rcvr = rcvrlist(rcvrnum - min(rcvrlist) + 1);
                receiver_offset = echo_sz * (rcvr-1);
                filename = sprintf('P%05d.7', header.rdb_hdr.run_int + ((slicenum-1) * nslices_per_pass));
                file_offset = header.total_length + phase_offset + receiver_offset + baseline_offset;
                fid = fopen( fullfile( dirname, filename), 'r', endianID);
                % execute below even if file does not exist, to keep same
                % workflow of fprintf information to user, especially status
                if fid<3
                    status_str = 'File not found';
                else
                    status = fseek(fid, file_offset, 'bof');
                    if status==-1
                        status_str = 'Offset outside of file';
                    else
                        [raw_data, count] = fread( fid, NumberDataElementsToFread, raw_data_type);
                        if count~=NumberDataElementsToFread
                            status_str = 'Could not read sufficient data';
                        else
                            status_str = 'OK';
                            raw_data = complex( raw_data(1:2:count), raw_data(2:2:count) );
                            raw_data = reshape( raw_data, da_xres, act_yres );
                            raw(:,:,phasenum,1,slicenum,rcvrnum) = transpose( raw_data );
                        end
                    end
                end
                
                fprintf(' %4d %4d %4d %10s %10d   %s\n', ...
                    slice_requested, phase_requested, rcvr, filename, file_offset, status_str);
                
                if fid>=3
                    fclose(fid);
                end
                
            end
        end
    end
    
elseif strcmp(temporal_type,'int')
    fprintf('  Slice Phase Rcvr Filename     Offset   Status\n');
    fprintf('  -----------------------------------------------------------\n');
    raw = zeros( act_yres, da_xres, max(phaselist), max(echolist), max(slicelist), max(rcvrlist) );
    for phasenum = min(phaselist):max(phaselist)
        phase_requested = phaselist(phasenum - min(phaselist) + 1);
        for rcvrnum = min(rcvrlist):max(rcvrlist)
            rcvr = rcvrlist(rcvrnum - min(rcvrlist) + 1);
            receiver_offset = echo_sz * (rcvr-1) * nslices_per_pass;
            for slicenum = min(slicelist):max(slicelist)
                slice_requested = slicelist(slicenum - min(slicelist) + 1);
                slice_offset = echo_sz * (header.data_acq_tab.slice_in_pass(slice_requested)-1);
                filename = sprintf('P%05d.7', header.rdb_hdr.run_int+header.data_acq_tab.pass_number(slice_requested));
                file_offset = header.total_length + slice_offset + receiver_offset + baseline_offset;
                fid = fopen( fullfile( dirname, filename), 'r', endianID);
                % execute below even if file does not exist, to keep same
                % workflow of fprintf information to user, especially status
                if fid<3
                    status_str = 'File not found';
                else
                    status = fseek(fid, file_offset, 'bof');
                    if status==-1
                        status_str = 'Offset outside of file';
                    else
                        [raw_data, count] = fread( fid, NumberDataElementsToFread, raw_data_type);
                        if count~=NumberDataElementsToFread
                            status_str = 'Could not read sufficient data';
                        else
                            status_str = 'OK';
                            raw_data = complex( raw_data(1:2:count), raw_data(2:2:count) );
                            raw_data = reshape( raw_data, da_xres, act_yres );
                            raw(:,:,phasenum,1,slicenum,rcvrnum) = transpose( raw_data );
                        end
                    end
                end
                
                fprintf(' %4d %4d %4d %10s %10d   %s\n', ...
                    slice_requested, phase_requested, rcvr, filename, file_offset, status_str);
                
                if fid>=3
                    fclose(fid);
                end
                
            end
        end
    end
    
elseif strcmp(temporal_type,'fMRI')
    fprintf('  Slice Phase Rcvr Filename     Offset   Status\n');
    fprintf('  -----------------------------------------------------------\n');
    raw = zeros( act_yres, da_xres, max(phaselist), max(echolist), max(slicelist), max(rcvrlist) );
    
    % filename does not change within for loops, so fopen file once outside of for loops;
    % fseek originally and always done from bof, so is not iteration-dependent
    
    filename = sprintf('P%05d.7', header.rdb_hdr.run_int);
    
    if exist( fullfile( dirname, filename) , 'file' )~=2
        errordlg('File not found. No data returned.','ERROR!','modal');
        raw = [];
        ec=1;
        return;
    end
    
    fid = fopen( fullfile( dirname, filename), 'r', endianID);
    
    if fid<3
        errordlg('File found but cannot be opened. No data returned.','ERROR!','modal');
        raw = [];
        ec=1;
        return;
    end
    
    for phasenum = min(phaselist):max(phaselist)
        phase_requested = phaselist(phasenum - min(phaselist) + 1);
        phase_offset = (phase_requested - 1) * nslices * max(rcvrlist) * echo_sz;
        for rcvrnum = min(rcvrlist):max(rcvrlist)
            rcvr = rcvrlist(rcvrnum - min(rcvrlist) + 1);
            receiver_offset = echo_sz * (rcvr-1) * nslices_per_pass;
            for slicenum = min(slicelist):max(slicelist)
                slice_requested = slicelist(slicenum - min(slicelist) + 1);
                slice_offset = echo_sz * (header.data_acq_tab.slice_in_pass(slice_requested) - 1);
                file_offset = header.total_length + phase_offset + slice_offset + receiver_offset + baseline_offset;
                
                status = fseek(fid, file_offset, 'bof');
                if status==-1
                    status_str = 'Offset outside of file';
                else
                    [raw_data, count] = fread( fid, NumberDataElementsToFread, raw_data_type);
                    if count~=NumberDataElementsToFread
                        status_str = 'Could not read sufficient data';
                    else
                        status_str = 'OK';
                        raw_data = complex( raw_data(1:2:count), raw_data(2:2:count) );
                        raw_data = reshape( raw_data, da_xres, act_yres );
                        raw(:,:,phasenum,1,slicenum,rcvrnum) = transpose( raw_data );
                    end
                end
                fprintf(' %4d %4d %4d %10s %10d   %s\n', ...
                    slice_requested, phase_requested, rcvr, filename, file_offset, status_str);
            end
        end
    end
    
    fclose(fid);
    
    if (header.rdb_hdr.ref == 3) %New fMRI with extra ref views per raw data set
        if strcmp(save_baseline,'db')
            if (header.rdb_hdr.extra_frames_top ~= 0)
                last_extra_view=header.rdb_hdr.extra_frames_top;
                raw(1:last_extra_view,:,:,:,:,:)=[];
            elseif (header.rdb_hdr.extra_frames_bot ~= 0)
                last_extra_view=act_yres-header.rdb_hdr.extra_frames_bot+1;
                raw(last_extra_view:act_yres,:,:,:,:,:)=[];
            end
        end
    end
end
return


%{
Revision History
Date (YYYY mmm DD) Author
Comment

2003-MAY-29  Matthew Eash.
Rev 2.3  "Initial version"

2003-DEC-17  Steve Huff.
Rev 2.4  "Added Multi Temporal Phase Capablity."

2004-FEB-25  Steve Huff.
Rev 2.5  "Numerous updates for specifying ranges."

2004-MAR=01  Steve Huff.
Rev 2.6  "Added non-standard P file reading and
         spead up stripping of extra views for fMRI".

2017-JUN-06  M.Miyazaki
PRICE    "Supported floating point raw data type"

2019-FEB-04  Andreas Ebel
VEGA     "HCSDM00540564: Added check of non_standard_pfile flag
          when checking for complete raw data set.
          Included rhsspsave in calculation of expected file size.
          Fixed indentation and deleted commented/unused code."

2019-APR-17 Barry Fetics
VEGA     "HCSDM00551525 improve fopen/fclose locations,
          prevent crash : Too many open files. Close files to prevent MATLAB instability."
          Collect golden data for, and test all 5 P file types
            - na, grafimage, int, seq, fMRI
          General clean-up.
%}
