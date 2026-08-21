function pfile2lcmodel(fname,dname,logfile)
%PFILE2LCMODEL  Convert GE raw data to LCModel raw data
%pfile2lcmodel(fname,dname,logfile)
%  fname  P-file or ScanArchive filename
%  dname  output directory name;                              ('P<runnum>')
%logfile  write output to logfile
%
% 3/2024 Rolf Schulte
if nargin<1, help(mfilename); return; end

if ~exist('dname','var'),    dname = []; end
if ~exist('logfile','var'),  logfile = ''; end



%% redirect output to logile
if ~isempty(logfile)
    if exist(logfile,'file')
        fprintf('logfile(=''%s'') existing; deleting\n',logfile);
        delete(logfile);
    end
    diary(logfile);
end


%% read GE raw data
[d,h] = read_p(fname);
if isempty(dname), dname = sprintf('P%.5d',h.image.rawrunnum); end


%% write LCModel raw + control files
sub_write_lcmodel(d,h,dname);


%% close again
if ~isempty(logfile)
    diary off
end


end      % pfile2lcmodel.m



%% sub-function sub_write_lcmodel
function sub_write_lcmodel(data,hdr,dname)
% WRITE_LCMODEL  Write data into ASCII file+header readable by LCModel
%    write_lcmodel(data,hdr,ref_data,fname,dim_spec)
%
%       data  time domain signal: pre-processed, but not be filtered
%        hdr  P-file header



%% adapt filename + create sub directories
dname(regexpi(dname,'\\')) = '/';   % replace Win dir sep through Linux
if exist(dname,'dir')~=7, mkdir(dname); end    
if exist([dname '/met'],'dir')~=7, mkdir([dname '/met']); end    
if exist([dname '/h2o'],'dir')~=7, mkdir([dname '/h2o']); end    


%% determine parameters
vol = hdr.rdb_hdr.roilenx*hdr.rdb_hdr.roileny*hdr.rdb_hdr.roilenz*1d-3;
% volume [mL]
f0 = hdr.rdb_hdr.ps_mps_freq*1d-7;          % B0 frequency [MHz]
TE = hdr.image.te*1d-3;
id = sprintf('P%g.7; E%g-S%g',hdr.image.rawrunnum,hdr.image.im_exno,...
    hdr.image.im_seno);
dt = 1/hdr.rdb_hdr.user0;                   % dwell time=1/BW [s]


%% coil-combine data
[nf,ns,n3,n4,n5,nc] = size(data);
if n3>1, warning('n3(=%d)>1',n3); end
if n4>1, warning('n4(=%d)>1',n4); end
if n5>1, warning('n5(=%d)>1',n5); end
if sub_ischop(hdr), data(2:2:end,:) = -data(2:2:end,:); end
if nc>1
    data = mrs_coil_combine(data);
end
nref = hdr.rdb_hdr.user19;
dmet = mean(data((nref+1):nf,:),1);
dref = mean(data(1:nref,:),1);


%% write metabolite data
fileid	= fopen([dname '/met/RAW'],'w');
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
for l2=1:ns
    %     fprintf(fileid,'  %e   %e\n',real(data(l)),imag(data(l)));
    str = sprintf('   %e   %e',real(dmet(1,l2)),imag(dmet(1,l2)));
    str = regexprep(str,'e\+0','e\+');
    str = regexprep(str,'e\-0','e\-');
    str = regexprep(str,' -','-');
    fprintf(fileid,'%s\n',str);    
end
fclose(fileid);


%% write reference (H2O) data
if ~isempty(dref)
    fileid	= fopen([dname '/h2o/RAW'],'w');
    fprintf(fileid,' $NMID\n');
    fprintf(fileid,' ID= ''%s''\n',id);         %
    fprintf(fileid,' BRUKER= F\n');             % complex conjugate data
    fprintf(fileid,' FMTDAT= ''(2e15.6)''\n');  % Fortran format specification
    fprintf(fileid,' VOLUME= %.3e\n',vol);      % volume [mL]
    fprintf(fileid,' TRAMP= 1.0\n');            % scaling factor for absoluete quantification
    fprintf(fileid,' $END\n');
    for l2=1:ns
        % fprintf(fileid,'  %e   %e\n',real(ref_data(l)),imag(ref_data(l)));
        str = sprintf('   %e   %e',real(dref(1,l2)),imag(dref(1,l2)));
        str = regexprep(str,'e\+0','e\+');
        str = regexprep(str,'e\-0','e\-');
        str = regexprep(str,' -','-');
        fprintf(fileid,'%s\n',str);
    end
    fclose(fileid);
end

%% write cpStart info file
fileid	= fopen([dname '/met/cpStart'],'w');
fprintf(fileid,'title= ''P%05d Exam%d Series%d %s\n',...
            hdr.image.rawrunnum,hdr.exam.ex_no,hdr.series.se_no,hdr.image.psd_iname);
fprintf(fileid,'filraw= ''%s/met/RAW''\n',dname);
fprintf(fileid,'filps= ''%s/ps''\n',dname);
fprintf(fileid,'hzpppm= %.4e\n',f0);
fprintf(fileid,'deltat= %.4e\n',dt);
fprintf(fileid,'nunfil= %d\n',ns);
fprintf(fileid,'sddegp= 1.\n');
fprintf(fileid,'iaverg= 1\n');
fprintf(fileid,'ndcols= 1\n');
fprintf(fileid,'ndrows= 1\n');
fprintf(fileid,'icolen= 1\n');
fprintf(fileid,'irowen= 1\n');
fprintf(fileid,'ndslic= 1\n');
fprintf(fileid,'sddegz= 3.\n');
fprintf(fileid,'doecc= T\n');
fprintf(fileid,'filh2o= ''%s/h2o/RAW''\n',dname);
fclose(fileid);



end      % sub_write_lcmodel


%% ischop
function chp = sub_ischop(h)
if ~isfield(h,'rdb_hdr'), error('~isfield(h,''rdb_hdr'')'); end
if ~isfield(h.rdb_hdr,'data_collect_type')
    error('~isfield(h.rdb_hdr,''data_collect_type'')'); 
end

if sub_isodd(h.rdb_hdr.data_collect_type)        % chopping
    chp = false;
else
    chp = true;
end

end      % sub_ischop


%% isodd
function bool = sub_isodd(number)
if ~isfloat(number), number = double(number); end
if number-floor(number) ~= 0
  warning('Input number not integer')
end

bool = logical(number/2-floor(number/2) ~= 0);
end      % sub_isodd

