function write_lcmodel(data,hdr,basis,ref_data,fname,cwdir,dim_spec)
% WRITE_LCMODEL  Write data into ASCII file+header readable by LCModel
%    write_lcmodel(data,hdr,basis,ref_data,fname,cwdir,dim_spec)
%
%       data  time domain signal: pre-processed, but not be filtered
%        hdr  P-file header
%      basis  file+pathname for LDModel basis file (Linux)
%   ref_data  reference data (w/o water suppression)
%      fname  output filename; if empty->lcm_matlab/P<#>
%      cwdir  current working directory
%   dim_spec  Dimension along acquisition (default=largest dim)
%
% 3/2015 Rolf Schulte
if nargin<1, help(mfilename); return; end

%% inputs
if ~exist('basis','var'),    basis = []; end
if ~exist('ref_data','var'), ref_data = []; end
if ~exist('fname','var'),    fname = []; end
if ~exist('cwdir','var'),    cwdir = []; end
if ~exist('dim_spec','var'), dim_spec = []; end
if isempty(dim_spec),        [tmp,dim_spec] = max(size(data)); end
ns = size(data,dim_spec);

if ~isempty(ref_data), 
    if(any(size(data)~=size(ref_data))),
        disp('size(data)');     disp(size(data));
        disp('size(ref_data)'); disp(size(ref_data));
        warning('size mismatch between metabolite and reference data');
    end
end

%% adapt filename + create sub directories
if isempty(cwdir), cwdir = [pwd '/']; end
cwdir(regexpi(cwdir,'\\')) = '/';
if ~isunix,
    switch lower(cwdir(1:2)),
        case 'h:', cwdir = ['/home/schulte' cwdir(3:end)];
        case 's:', cwdir = ['/projects/spectro' cwdir(3:end)];
        case 'k:', cwdir = ['/projects/hyperpol_eu' cwdir(3:end)];
        otherwise, warning('No replacement for drive (%s) found',cwdir(1:3));
    end
end

if isempty(fname), 
    fname = sprintf('lcm_matlab/P%.5d/P%.5d',hdr.image.rawrunnum,hdr.image.rawrunnum);
else
    if ~isempty(regexpi(fname,'\.raw$')), fname = fname(1:end-4); end
    fname(regexpi(fname,'\\')) = '/';   % replace Win dir sep through Linux
end
dirsep = regexpi(fname,'/');
for l=1:length(dirsep),
    dname = fname(1:dirsep(l));
    if exist(dname,'dir')~=7, mkdir(dname); end    
end


%% determine parameters
vol = hdr.rdb_hdr.roilenx*hdr.rdb_hdr.roileny*hdr.rdb_hdr.roilenz*1d-3;
% volume [mL]
f0 = hdr.rdb_hdr.ps_mps_freq*1d-7;          % B0 frequency [MHz]
TE = hdr.image.te*1d-3;
id = sprintf('P%g.7; E%g-S%g',hdr.image.rawrunnum,hdr.image.im_exno,...
    hdr.image.im_seno);
dt = 1/hdr.rdb_hdr.user0;                   % dwell time=1/BW [s]

%% write metabolite data
fileid	= fopen([fname '.raw'],'w');
fprintf(fileid,' $SEQPAR\n');
fprintf(fileid,' HZPPPM=  %.4e\n',f0);      % B0 frequency [MHz]
fprintf(fileid,' ECHOT= %.2f\n',TE);        % echo time [ms]
fprintf(fileid,' SEQ= ''PRESS''\n');        % must be PRESS or STEAM
fprintf(fileid,' $END\n');
fprintf(fileid,' $NMID\n');
fprintf(fileid,' ID= ''%s''\n',id);         % 
fprintf(fileid,' BRUKER= F\n');             % complex conjugate data
fprintf(fileid,' FMTDAT= ''(2e15.6)''\n');  % Fortran format specification
fprintf(fileid,' VOLUME= %.3e\n',vol);      % volume [mL]
fprintf(fileid,' TRAMP= 1.0\n');            % scaling factor for absoluete quantification
fprintf(fileid,' $END\n');
for l=1:ns,
    %     fprintf(fileid,'  %e   %e\n',real(data(l)),imag(data(l)));
    str = sprintf('   %e   %e',real(data(l)),imag(data(l)));
    str = regexprep(str,'e\+0','e\+');
    str = regexprep(str,'e\-0','e\-');
    str = regexprep(str,' -','-');
    fprintf(fileid,'%s\n',str);    
end
fclose(fileid);

%% write reference (H2O) data
if ~isempty(ref_data),
    fileid	= fopen([fname '_h2o.raw'],'w');
    fprintf(fileid,' $NMID\n');
    fprintf(fileid,' ID= ''%s''\n',id);         %
    fprintf(fileid,' BRUKER= F\n');             % complex conjugate data
    fprintf(fileid,' FMTDAT= ''(2e15.6)''\n');  % Fortran format specification
    fprintf(fileid,' VOLUME= %.3e\n',vol);      % volume [mL]
    fprintf(fileid,' TRAMP= 1.0\n');            % scaling factor for absoluete quantification
    fprintf(fileid,' $END\n');
    for l=1:ns,
        % fprintf(fileid,'  %e   %e\n',real(ref_data(l)),imag(ref_data(l)));
        str = sprintf('   %e   %e',real(ref_data(l)),imag(ref_data(l)));
        str = regexprep(str,'e\+0','e\+');
        str = regexprep(str,'e\-0','e\-');
        str = regexprep(str,' -','-');
        fprintf(fileid,'%s\n',str);
    end
    fclose(fileid);
end


%% write control file
if ~isempty(basis),
    fileid	= fopen([fname '.control'],'w');
    fprintf(fileid,' $LCMODL\n');
    fprintf(fileid,' TITLE=''%s; %s %s; %s; TE/TR/NS=%g/%g/%g; %gmL; %s''\n',...
        id,hdr.rdb_hdr.scan_date,hdr.rdb_hdr.scan_time,...
        hdr.image.psd_iname,TE,hdr.image.tr*1d-3,hdr.image.user4,vol,...
        hdr.exam.hospname);
    fprintf(fileid,' FILBAS=''%s''\n',basis);
    fprintf(fileid,' FILRAW=''%s.raw''\n',[cwdir fname]);
    if ~isempty(ref_data), 
        fprintf(fileid,' FILH2O=''%s_h2o.raw''\n',[cwdir fname]); 
        doecc = 'T';
    else
        doecc = 'F';
    end
    fprintf(fileid,' FILPS=''%s.ps''\n',[cwdir fname]);
    fprintf(fileid,' PGNORM=''A4''\n');
    fprintf(fileid,' FILTAB=''%s.table''\n',[cwdir fname]);
    fprintf(fileid,' LTABLE=7\n');
    fprintf(fileid,' HZPPPM=%g\n',f0);      % B0 frequency [MHz]
    fprintf(fileid,' NUNFIL=%g\n',ns);      % acq points
    fprintf(fileid,' DELTAT=%g\n',dt);      % dwell time=1/BW
    fprintf(fileid,' DOECC=%s\n',doecc);    % eddy current correction
    fprintf(fileid,' NSIMUL=12\n');         % #simulated basis spectra (1-12=lipids)
    fprintf(fileid,' DOWS=F\n');            % scaling to water
    fprintf(fileid,' SDDEGZ=3.\n');
    fprintf(fileid,' SDDEGP=1.\n');
    fprintf(fileid,' $END\n');
    fclose(fileid);
else
    % warning('No basis file specified: skipping control file'); 
end