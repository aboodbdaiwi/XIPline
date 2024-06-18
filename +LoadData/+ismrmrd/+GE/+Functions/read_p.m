function [data,header] = read_p(fname,do_single,use_ox,correction,verb)
% READ_P Read GE rawdata: P-files (*.7 or *.7.gz) or hdf5 files (*.h5)
%
%[data,header] = read_p(fname,do_single,use_ox,correction,verb)
%                                                             (default)
%        fname   filename (full name with/without directory,'*.h5','*.7')
%    do_single   read data in single precision                (false)
%       use_ox   Use Orchestra for reading p-files            (false)
%   correction   apply correction: []->determine from header  ([])
%            0   none
%            1   ptx: if ((header.rdb_hdr.user2==13) && (f0>35d6) && 
%                         (header.rdb_hdr.rdbm_rev<15))
%                correct for updown convert: f0=160d6-f0 & data=conj(data)
%            2   mr901: if
%                f0=f0+1931067306
%            3   reverse odd burzte acquisitions
%         verb   Verbose mode                                 (true)
%
% 7/2020  Rolf Schulte
% See also read_MR_rawdata, read_MR_headers, read_archive, archive2header.
if (nargout<1), help(mfilename); return; end

%% default input
if ~exist('fname','var'),       fname = []; end
if ~exist('do_single','var'),   do_single = []; end
if isempty(do_single),          do_single = false; end
if ~exist('use_ox','var'),      use_ox = []; end
if isempty(use_ox),             use_ox = false; end
if ~exist('corrections','var'), correction = []; end
if ~exist('verb','var'),        verb = []; end
if isempty(verb),               verb = true; end

f0corrmr901 = 193106730.6;        % correction freq for MR901
if verb>0, timerVal = tic; end    % record execution time


%% look for suitable file to load
if isempty(fname), fname = '*.7'; end            % default=look for p-file
if exist(fname)==7, fname = [fname '/*.7']; end  % if directory -> look for p-file 
if ~isempty(strfind(fname,'*'))
    if isunix
        fname = ls(fname);
        if isempty(fname), error('file %s not found',fname); end
        if ~isempty(regexp(fname,'\n','once'))
            fname = fname(1:(regexp(fname,'\n','once')-1));
            warning('multiple files found; chosing fname=%s',fname);
        end
        fname = regexprep(fname,'\*','');
    else
        xx = dir(fname);
        if length(xx)>1
            warning('multiple files found; chosing first');
        end
        if isempty(xx)
            error('file ''%s'' not found',fname);
        end
        fname = [xx(1).folder '\' xx(1).name];
    end    
    if ~exist(fname,'file'), error('Bug: file %s not found',fname); end
    if verb>0, fprintf('Loading fname=\n\t%s\n',fname); end
end


%% check if file found
if ~exist(fname,'file'), error('file %s not found',fname); end


%% if suffix=='.h5' read as hdf5 file from image archive
if ~isempty(regexpi(fname,'\.h5$'))
    [data,header] = LoadData.ismrmrd.GE.Functions.read_archive(fname,verb,do_single);
    if verb>0
        fprintf('read_p.m finished: ');
        fprintf('total reading time = %.2f [s]\n',toc(timerVal));
    end
    return
end


%% enable gunzip for Pxxx.7.gz
removefile = false;
if ~isempty(regexpi(fname,'\.gz$'))
    if verb>0, fprintf('Unzipping first\n'); end
    if exist(fname(1:end-3),'file') 
        error('gunzip-ped file %s already existing',fname); 
    end
    gunzip(fname);
    removefile = true;
    oldfname = fname;
    fname = fname(1:end-3);
end


%% check header revision
fid = fopen(fname,'r','l');
rdbm_rev = fread(fid, 1, 'float32');
rdbm_rev = round(rdbm_rev*1000)/1000;
fclose(fid);
if verb, fprintf('rdbm_rev=%g\n',rdbm_rev); end


%% reading p-file
use_ox = 1;
if (exist('GERecon','file') && use_ox)
    if verb, fprintf('Reading p-file using Orchstra (GERecon)\n'); end
    [data,header] = LoadData.ismrmrd.GE.Functions.read_pfile(fname,do_single,verb);
else
    if verb
        fprintf('Reading p-file using Matlab (read_MR_rawdata.m) routines\n');
    end
    [data,header] = LoadData.ismrmrd.GE.Functions.read_MR_rawdata(fname);
    if do_single
        if verb
            fprintf('Warning: read_MR_rawdata uses double\n');
            fprintf('\tconverting to single afterwards\n');
        end
        data = single(data);
    end
end


%% cleaning up
if removefile
    if exist(oldfname,'file')
        delete(fname); 
    else
        error('old file (%s) not existing anymore',oldfname);
    end
end


%% check headers if correction necessary
if isempty(correction)
    f0 = header.rdb_hdr.ps_mps_freq/10;
    correction = 0;
    % ptx
    if ((header.rdb_hdr.user2==13) && (f0>35d6) && (header.rdb_hdr.rdbm_rev<15))
        correction = 1;
    end
    % MR901
    if ((f0<0) && (round(header.rdb_hdr.rdbm_rev*10)==144))
        correction = 2;
    end
end


%% apply miscelaneous correcions
switch correction
    case 0
    case 1         % PTX updown converter for 13C
        fprintf('\nAttention: up-down-converter: Nuc=13,f0=%d\n',round(f0));
        fprintf('\t-> f0=160d6-f0  and  data=conj(data)\n\n');
        header.rdb_hdr.ps_mps_freq = (160d6 - f0)*10;
        data = conj(data);
    case 2         % MR901
        fprintf('\nAttention: MR901; correcting f0=f0(=%g)+f0cor(=%g)=%g\n\n',...
            round(f0),round(f0corrmr901),round(f0+f0corrmr901));        
        header.rdb_hdr.ps_mps_freq = (f0+f0corrmr901)*10;
    case 3
        fprintf('\nAttention: reversing odd burzte acquisitions\n');
    otherwise
        fprintf('Warning: correction(=%g) not found\n',correction);
end


%% fill missing fields
[n1,n2,n3,n4,n5,n6] = size(data);
if ~isfield(header.image,'specnuc')
    fprintf('\nAttention: h.image.specnuc field missing; copying from h.rdb_hdr.user2\n');
    header.image.specnuc = header.rdb_hdr.user2;
end
if ~isempty(regexpi(header.image.psd_iname,'3drad')) || ...
        ~isempty(regexpi(header.image.psd_iname,'silen')) || ...
        ~isempty(regexpi(header.image.psd_iname,'burzte'))
    if ~isfield(header.rdb_hdr,'nspokes_lowres')
        header.rdb_hdr.nspokes_lowres = header.rdb_hdr.user7; 
    end
    if ~isfield(header.rdb_hdr,'nspokes_highres')
        header.rdb_hdr.nspokes_highres = header.rdb_hdr.user8;
        if (header.rdb_hdr.nspokes_highres==0)
            warning('header.rdb_hdr.user8==0');
            fprintf('setting header.rdb_hdr.nspokes_highres to\n');
            rhnframes = header.rdb_hdr.nframes;
            opslquant = ceil(header.image.dfov/header.image.slthick/2)*2;
            header.rdb_hdr.nspokes_highres = rhnframes*opslquant;
            fprintf('\trhnframes*opslquant (=%g)\n',...
                header.rdb_hdr.nspokes_highres);
        end
    end
    if ~isfield(header.rdb_hdr,'noncart_dual_traj')
        header.rdb_hdr.noncart_dual_traj = header.rdb_hdr.user5;
    end
    if ~isfield(header.rdb_hdr,'noncart_traj_kmax_ratio')
        header.rdb_hdr.noncart_traj_kmax_ratio = header.rdb_hdr.user6;
    end
    if ~isfield(header.rdb_hdr,'oversamplingfactor')
        header.rdb_hdr.oversamplingfactor = header.rdb_hdr.user33;
    end
    if ~isfield(header.image,'spokesPerSeg')
        header.image.spokesPerSeg = header.rdb_hdr.user34;
    end
    if correction==3
        if header.rdb_hdr.nechoes>1
            if isfield(header,'format')
                warning('Bug fix for ZTE-BURST p-files: time reversing even acqs\n');
                data(:,:,:,2:2:end,:) = data(:,end:-1:1,:,2:2:end,:);
            end
        end
    end
    if header.image.user8>1       % multi-phase
        fprintf('Multi-phase 3dradial data: reshaping slices into dim3\n');
        np = header.image.user8;
        if n3~=1
            warning('n3(=%g)~=1, skipping reshaping',n3);
        else
            if abs(n5/np-floor(n5/np))>1d-10
                warning('n5(=%g) not divisble by np(=%g); skipping',n5,np);
            else
                dtmp = complex(zeros(n1,n2,np,n4,n5/np,n6,class(data)));
                for lp=1:np
                    ii = (1:(n5/np))+(lp-1)*(n5/np);
                    dtmp(:,:,lp,:,:,:) = data(:,:,1,:,ii,:);
                end
                data = dtmp;
                clear dtmp;
                % data = reshape(data,[n1 n2 np n4 n5/np n6]);
            end
        end
    end
end


%% checking data
if header.rdb_hdr.nframes~=n1
    warning('h.rdb_hdr.nframes(=%d)~=n1(=%d)',header.rdb_hdr.nframes,n1);
end
if header.rdb_hdr.frame_size~=n2
    warning('h.rdb_hdr.frame_size(=%d)~=n2(=%d)',header.rdb_hdr.frame_size,n2);
end
if header.image.numecho~=n4
    warning('h.image.numecho(=%d)~=n4(=%d)',header.image.numecho,n4);
end
if header.rdb_hdr.nslices~=n5
    warning('h.rdb_hdr.nslices(=%d)~=n5(=%d)',header.rdb_hdr.nslices,n5);
end
ncoils = (header.rdb_hdr.dab(2)-header.rdb_hdr.dab(1))+1; % number of Rx coil elements
if ncoils~=n6
    warning('ncoils(=%d)~=n6(=%d)',ncoils,n6);
end


%% print execution time
if verb>0
    fprintf('read_p.m finished: ');
    fprintf('total reading time = %.2f [s]\n',toc(timerVal));
end

end    % main function read_p.m
