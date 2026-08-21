function [data,header] = read_p(fname,do_single,use_ox,correction,verb,...
    anon,invmns,excoil)
% READ_P Read GE rawdata: P-files (*.7 or *.7.gz) or hdf5 files (*.h5)
%[data,header] = read_p(fname,do_single,use_ox,correction,verb,anon,...
%                       invmns,excoil)
%        fname   Filename (full name with/without directory,'*.h5','*.7')
%    do_single   Read data in single precision                       (true)
%       use_ox   Use Orchestra for reading raw data                 (false)
%                if error: automatically tries other approach
%   correction   Apply correction: []->determine from header           ([])
%            0   none
%            1   ptx: if ((header.rdb_hdr.user2==13) && (f0>35d6) && 
%                         (header.rdb_hdr.rdbm_rev<15))
%                correct for updown convert: f0=160d6-f0 & data=conj(data)
%            2   mr901: if
%                f0=f0+1931067306
%            3   reverse odd burzte acquisitions
%         verb   Verbose mode                                        (true)
%         anon   Anonymise header                                   (false)
%                (true, if environment variable mnsrpanon set)
%       invmns   Invert freq axis for selected nuclei -> conj(data)    ([])
%                read from environment variable mnsrpinvmns if present
%       excoil   Exclude coil elements                                 ([])
%                read from environment variable mnsrpexcoil if present
%
%         data   Raw data: size=[YRES XRES PHASES ECHOES SLICES RECEIVERS]
%       header   Header structure
%
% 12/2024 Rolf Schulte
% See also read_MR_rawdata, read_MR_headers, read_archive, archive2header,
% read_arc.
if (nargout<1), help(mfilename); return; end


%% default input
if ~exist('fname','var'),       fname = []; end
if ~exist('do_single','var'),   do_single = []; end
if isempty(do_single),          do_single = true; end
if ~exist('use_ox','var'),      use_ox = []; end
if isempty(use_ox),             use_ox = false; end
if ~exist('corrections','var'), correction = []; end
if ~exist('verb','var'),        verb = []; end
if isempty(verb),               verb = true; end
if ~exist('anon','var'),        anon = []; end
if isempty(anon)
    if isempty(getenv('mnsrpanon'))
        anon = false;
    else
        anon = true;
    end
end
if ~exist('invmns','var'),      invmns = []; end
mnsrpinvmns = getenv('mnsrpinvmns');
if islogical(invmns)
    if invmns, invmns = -1; end
end
if isempty(invmns) && ~isempty(mnsrpinvmns)
    invmns = str2num(mnsrpinvmns);
end
if ~exist('excoil','var'),      excoil = []; end
mnsrpexcoil = getenv('mnsrpexcoil');
if isempty(excoil) && ~isempty(mnsrpexcoil)
    excoil = str2num(mnsrpexcoil);
end

f0corrmr901 = 193106730.6;        % correction freq for MR901
default_pfile = false;
if verb>0, timerVal = tic; end    % record execution time


%% check reading routines
if ~exist('GERecon','file') && ~exist('read_MR_headers','file')
    warning('GERecon && read_MR_headers not found');
end


%% look for suitable file to load
if isempty(fname)
    if default_pfile
        if ~isempty(dir('*.7'))   % load p-file by default
            fname = '*.7';
        else
            fname = 'ScanArchive_*.h5';
        end
    else
        if ~isempty(dir('ScanArchive_*.h5'))  % load ScanArchives by default
            fname = 'ScanArchive_*.h5';
        else
            fname = '*.7';
        end
    end
end            
if exist(fname,'dir')   % if directory -> look for p-file or ScanArchive
    if default_pfile
        fname = [fname '/*.7'];
        if isempty(dir(fname))
            fname = [fname(1:end-1) 'h5'];
        end
    else
        fname = [fname '/ScanArchive_*.h5'];
        if isempty(dir(fname))
            fname = [fname(1:end-2) '7'];
        end
    end
end  


if ~isempty(regexpi(fname,'*',  'once'))
    xx = dir(fname);
    if isempty(xx), error('no data found: isempty(xx)'); end
    if ~isempty(regexpi(xx(end).name,'mrdc\.7$',  'once'))
        xx = xx(1:end-1);
    end
    if isempty(xx), error('no pfile found; dicom discarded'); end
    if length(xx)>1
        warning('multiple files found; chosing last');
        for l=1:length(xx), fprintf('\t%s\n',xx(l).name); end
    end
    if isempty(xx)
        error('file ''%s'' not found',fname);
    end
    if isfield(xx(end),'folder')
        folder = xx(end).folder;
    else
        folder = pwd;
    end
    fname = [folder filesep xx(end).name];
    if ~exist(fname,'file'), error('file %s not found',fname); end
    if verb>0, fprintf('Loading fname=\n\t%s\n',fname); end
end


%% print info
if verb
    fprintf('read_p: fname, do_single=%d, use_ox=%d, correction=%d, ',...
        do_single,use_ox,correction);
    fprintf('verb=%d, anon=%d, invmns=%s, excoil=%s\n',...
        verb,anon,num2str(invmns),num2str(excoil));
end


%% check if file found
if ~exist(fname,'file'), error('file %s not found',fname); end


%% if suffix=='.h5' read as hdf5 file from image archive
if ~isempty(regexpi(fname,'\.h5$', 'once'))
    % if ~use_ox && (get_mver<9)
    %     warning('~use_ox && (get_mver(=%g)<9): setting use_ox=true',get_mver);
    %     if ~exist('GERecon','file'), error('GERecon not found'); end
    %     use_ox = true;
    % end
    if use_ox
        try
            if verb, fprintf('read_archive (Ox based) ***** \n'); end
            [data,header] = LoadData.GE_Recon.fidall.read_archive(fname,verb,do_single);
        catch ME
            warning('read_archive Ox failed: reverting back to Matlab');
            fprintf('Error message: %s\n',ME.message);
            if verb, fprintf('read_arc (Matlab based) ***** \n'); end
            [data,header] = LoadData.GE_Recon.fidall.read_arc(fname,verb);
            if ~do_single, data = double(data); end
        end
    else
        try
            if verb, fprintf('read_arc (Matlab based) ***** \n'); end
            [data,header] = LoadData.GE_Recon.fidall.read_arc(fname,verb);
            if ~do_single, data = double(data); end
        catch ME
            warning('read_arc Matlab failed: reverting back to Ox');
            fprintf('Error message: %s\n',ME.message);
            if verb, fprintf('read_archive (Ox based) ***** \n'); end
            [data,header] = LoadData.GE_Recon.fidall.read_archive(fname,verb,do_single);
        end
    end

else           % if h5-file -> p-file

    %% enable gunzip for Pxxx.7.gz
    removefile = false;
    if ~isempty(regexpi(fname,'\.gz$', 'once'))
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
    if use_ox
        try
            if verb, fprintf('read_pfile (Ox based) ***** \n'); end
            [data,header] = LoadData.GE_Recon.fidall.read_pfile(fname,do_single,verb);
        catch ME
            warning('read_pfile Ox failed: reverting back to Matlab');
            fprintf('Error message: %s\n',ME.message);
            if verb, fprintf('read_MR_rawdata (Matlab based) ***** \n'); end
            [data,header] = LoadData.GE_Recon.fidall.read_MR_rawdata(fname);
            if do_single, data = single(data); end
        end
    else
        try
            if verb, fprintf('read_MR_rawdata (Matlab based) ***** \n'); end
            [data,header] = LoadData.GE_Recon.fidall.read_MR_rawdata(fname);
            if do_single, data = single(data); end
        catch ME
            warning('read_MR_rawdata Matlab failed: reverting back to Ox');
            fprintf('Error message: %s\n',ME.message);
            if verb, fprintf('read_pfile (Ox based) ***** \n'); end
            [data,header] = LoadData.GE_Recon.fidall.read_pfile(fname,do_single,verb);
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

end      % if h5-file
header.fname = fname;


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


%% apply miscellaneous corrections
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
if ~isempty(regexpi(header.image.psd_iname,'3drad', 'once')) || ...
        ~isempty(regexpi(header.image.psd_iname,'silen', 'once')) || ...
        ~isempty(regexpi(header.image.psd_iname,'burzte', 'once'))
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


%% convert int header to char (for Matlab reading routines)
if ~use_ox
    hconv.exam   = {'patnameff','patidff','uniq_sys_id','service_id'};
    hconv.series = {'pure_cfg_params','se_suid'};
    hconv.image  = {'anatomy','GEcname','im_suid'};
    header = sub_hdr_int2char(header,hconv);
end


%% reshape if multi-phase stored in slice dimension
if ~isempty(regexpi(header.image.psd_iname,'fidspiral', 'once')) 
    nrep = header.image.user22;
    nsl = header.image.slquant;
    if nrep > 1
        if n4 == 1
            if n5 == nrep*nsl
                fprintf('Reshaping data: slices([%d,%d,%d,%d,%d,%d]) ',...
                    n1,n2,n3,n4,n5,n6);
                fprintf('-> echoes&slices ([%d,%d,%d,%d,%d,%d])\n',...
                    n1,n2,n3,nrep,nsl,n6);
                data = permute(reshape(data,[n1 n2 n3 nsl nrep n6]),[1 2 3 5 4 6]);
            else
                warning('size(data,5)(=%d) ~= nrep(=%d)*nsl(=%d)',n5,nrep,nsl);
            end
        else
            warning('size(data,4)(=%d)>1 and nrep(=%d)>1',n4,nrep);
        end
    end
end


%% inverting freq axis
if ~isempty(invmns)
    if any(invmns==header.image.specnuc) || any(invmns==-1)
        fprintf('Attention: data=conj(data);\n');
        data = conj(data);
    end
end


%% excluding coil elements
if ~isempty(excoil)
    fprintf('Attention: excluding coil elements: %s\n',num2str(excoil));
    nc = size(data,6);
    iic = 1:nc;
    for l=1:length(excoil), iic = iic(iic~=excoil(l)); end
    data = data(:,:,:,:,:,iic);
end


%% checking data size
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
if isfield(header,'cttentry')
    if isfield(header.cttentry,'cttentry')
        if isfield(header.cttentry.cttentry,'numChannels')
            if ncoils~=header.cttentry.cttentry(1).numChannels
                warning('ncoils(=%d)~=numChannels(=%d)',...
                    ncoils,header.cttentry.cttentry(1).numChannels);
            end
        end
    end

    % cttentry not really used + large -> remove
    header = rmfield(header,'cttentry');
end


%% checking correct R2
R2 = header.rdb_hdr.ps_mps_r2;
if isfield(header.rdb_hdr,'dacq_data_type')
    switch header.rdb_hdr.dacq_data_type
        case 1     % 'int16'
            if ~any(R2==[13 14 15])
                warning('EDR off; R2(=%g)~=13, 14 or 15',R2); pause(1);
            end
        case 2     % 'int32'
            if ~any(R2==[28 29 30])
                warning('EDR on; R2(=%g)~=28, 29 or 30',R2); pause(1);
            end
        otherwise
            if ~any(R2==[13 14 15 28 29 30])
                warning('R2(=%g)~=14,15,29 or 30',R2); pause(1);
            end
    end
else
    if ~any(R2==[13 14 15 28 29 30])
        warning('R2(=%g)~=13,14,15,28,29 or 30',R2); pause(1);
    end
end


%% checking data
if isempty(data), warning('data is empty'); end
if any(any(any(any(any(any(isnan(data)))))))
    % if any(isnan(data(:)))
    warning('data contains NaN; replacing by 0');
    data(isnan(data)) = 0;
end
if any(any(any(any(any(any(isinf(data)))))))
    % if any(isinf(data(:)))
    warning('data contains inf; replacing by 0');
    data(isinf(data)) = 0;
end


%% anonymise header
if anon, header = header_anonymise(header); end


%% print execution time
if verb>0, fprintf('read_p: runtime = %.2f [s]\n',toc(timerVal)); end


end      % read_p.m


%% sub-functions
function h = sub_hdr_int2char(h,hconv)
flds = fieldnames(hconv);
for l1=1:length(flds)
    if ~isfield(h,flds{l1}), error('header.%s not existing',flds{l1}); end
    for l2=1:length(hconv.(flds{l1}))
        if isfield(h.(flds{l1}),hconv.(flds{l1}){l2})
            ent = h.(flds{l1}).(hconv.(flds{l1}){l2});
            if ~ischar(ent)
                if any(ent<0)
                    warning('h.%s.%s < 0',(flds{l1}),(hconv.(flds{l1}){l2}));
                end
                 if any(ent>256)
                    warning('h.%s.%s > 256',(flds{l1}),(hconv.(flds{l1}){l2}));
                end
                h.(flds{l1}).(hconv.(flds{l1}){l2}) = deblank(char(ent));
            else
                % fprintf('~ischar(h.%s.%s)\n',(flds{l1}),(hconv.(flds{l1}){l2}));
            end
        else
            warning('header.%s: no field ''%s''',flds{l1},hconv.(flds{l1}){l2});
        end
    end
end

end      % sub_hdr_int2char
