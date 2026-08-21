function header2readme(what,sv)
%HEADER2README Extract info from header for README file
% header2readme(what,sv,image_user)
%      what   Selection of different field: 
%             cell structure with header field names
%             1=3dradial, 2=fidall fields                               (2)
%        sv   Write output to README_hdr.txt file                    (true)
%
% 12/2024 Rolf Schulte
% if (nargin<1), help(mfilename); return; end

if ~exist('what','var'), what = []; end
if isempty(what),        what = 2; end
if ~exist('sv','var'),   sv = []; end
if isempty(sv),          sv = true; end


%% check filename first
if sv
    fname_rdm = [pwd filesep 'README_hdr.txt'];
    fname_old = [pwd filesep 'README_hdr_old.txt'];
    if exist(fname_rdm,'file')
        warning('file ''%s'' exists; renaming',fname_rdm);
        if exist(fname_old,'file')
            warning('deleting file ''%s''',fname_old);
            delete(fname_old);
        end
        movefile(fname_rdm,fname_old);
    end
end


%% selection of values to print
if iscell(what), flds = what; what = 0; end
switch what
    case 0
    case 1
        flds = {'pfname','series','psd_iname','slthick','fov','res',...
            'user10','user12','user13','etl','spokesPerSeg','tr',...
            'flip','bw','sctime','user36'};
    case 2
        flds = {'series','se_desc','runnum','slthick','fov','tr','te',...
            'flip','sctime','specnuc',...
            'user0','user1','user3','user4','user14','user19',...
            'user20','user21','nex','TG','R1','R2','data_type','psdname'};
    otherwise, error('what unknown');
end


%% read header
hdr = [dir('**/ScanArchive*.h5') ; dir('**/*.7')];
iinc = [];
for l=1:length(hdr)
    if isempty(regexpi(hdr(l).name,'noisestatistics','once'))
        iinc = [iinc,l];
    end
end
hdr = hdr(iinc);
inc = 0;
for l=1:length(hdr)
    htmp = [];
    fname = [hdr(l).folder filesep hdr(l).name];
    if exist(fname,'file')
        fprintf('Reading header %s\n',fname);
        if ~isempty(regexpi(fname,'\.7$','once'))
            if isempty(regexpi(fname,'\MRDC.7$','once'))   % excluding dcm
                htmp = read_MR_headers(fname);
            end
        else
            if ~isempty(regexpi(fname,'\.h5$','once'))
                % htmp = read_archive_header(fname);
                try
                    htmp = read_arc_header(fname);
                catch ME
                    warning(sprintf('read_arc_header failed\n\t%s',ME.message));
                end
            else
                warning('fname(=%s) not ending with .7 or .h5',fname);
            end
        end
    else
        warning('strange input'); disp(fname);
    end
    if ~isempty(htmp)
        inc = inc+1;
        h{inc} = htmp;
    end
end
if inc==0, error('no suitable files found'); end


%% check if fidall
if what==2
    do_warn = true;
    for l=1:length(h)
        if isempty(regexpi(h{l}.image.psd_iname,'fidall','once')) && ...
                isempty(regexpi(h{l}.image.psd_iname,'fidcsi','once')) && ...
                do_warn
            warning('not fidall/fidcsi: l=%d psd_iname=%s',l,h{l}.image.psd_iname);
            do_warn = false;
        end
    end
end


%% exclude double files
for l=1:length(h), tt(l) = h{l}.image.im_datetime; end
[~,i2] = unique(tt);
if length(tt)~=length(i2)
    fprintf('Raw data not unique; excluding duplicates\n');
    fprintf('length(tt) = %d; length(i2) = %d\n',length(tt),length(i2));
    for l=1:length(i2), tmp{l} = h{i2(l)}; end
    h = tmp;
    clear tmp
    tt = tt(i2);
end


%% sort according to time when data was acquired
[~,ii] = sort(tt);
for l=1:length(h), tmp{l} = h{ii(l)}; end
h = tmp;
clear tmp


%% exclude corrupted files (series==65537)
l2 = 0;
clear tmp
for l1=1:length(h)
    if h{l1}.series.se_no ~= 65537
        l2 = l2+1;
        tmp{l2} = h{l1};
    end
end
nexcl = length(h)-length(tmp);
if nexcl > 0, fprintf('Excluding %d corrupted headers\n',nexcl); end
h = tmp;
clear tmp


%% extract info
if exist('h','var')
    [se_str,ex_str] = sub_extractinfo(h,flds);
else
    warning('no files found');
    se_str{1} = '';
    return;
end


%% print README info
if sv, fid = fopen(fname_rdm,'wt'); end
switch what
    case 1
        fprintf('\nuser10: Dual drive channel selection: 0=off(normal),1=I,2=Q,3=I&Q,4=I+Q\n');
        fprintf('user12: burst_mode 0=off;1=GE;2=SE; 3=2SE;4=GESE;5=GE2SE\n');
        fprintf('user13: burst_rf7file 0=Block,1=SLR,2=HS1,3=Mitschang\n');
        fprintf('user36: TE [ms]\n\n');
end

for lc=1:length(se_str)
    % fprintf('length(ex_str{%d})=%d\n',lc,length(ex_str{lc}));
    if ~isempty(ex_str{lc})
        % exam info
        if length(ex_str{lc})>60
            fprintf('\n%s\n',ex_str{lc});
            if sv, fprintf(fid,'\n%s\n',ex_str{lc}); end
            % series info title
            for lf=1:length(flds)
                fprintf('%s ',flds{lf});
                if sv, fprintf(fid,'%s ',flds{lf}); end
            end
            fprintf('\n');
            if sv, fprintf(fid,'\n'); end
        else
            fprintf('%s',ex_str{lc});
            if sv, fprintf(fid,'%s',ex_str{lc}); end
        end
    end

    fprintf('%s\n',se_str{lc});
    if sv, fprintf(fid,'%s\n',se_str{lc}); end
end



%% closing
if sv, fclose(fid); end
if nargout==0, clear se_str ex_str; end

end      % header2readme.m


%% sub-functions
function [se_str,ex_str] = sub_extractinfo(hc,flds)

% f0 = h.rdb_hdr.ps_mps_freq/10;
% TE = h.rdb_hdr.te/1d3;
% TR = h.image.tr*1d-6;
% r1 = h.rdb_hdr.ps_mps_r1;
% r2 = h.rdb_hdr.ps_mps_r2;
% tg = h.rdb_hdr.ps_mps_tg;

%% get fields
len = struct;
for lc=1:length(hc)
    h = hc{lc};
    % exam info
    w{lc}.ex_no   = h.series.se_exno;
    w{lc}.scan_date = h.rdb_hdr.scan_date;
    w{lc}.patid   = deblank(char(h.exam.patidff));
    w{lc}.patname = deblank(char(h.exam.patnameff));
    w{lc}.cname   = deblank(char(h.image.cname));

    % series info
    w{lc}.pfname  = sprintf('P%05d.7',h.image.rawrunnum);
    w{lc}.runnum  = h.image.rawrunnum;
    w{lc}.series  = h.series.se_no;
    w{lc}.psd_iname = h.image.psd_iname;
    w{lc}.psdname   = h.image.psdname;
    if isfield(h.image,'specnuc')
        w{lc}.specnuc = h.image.specnuc;
    else
        warning('h.image.specnuc not existing; using h.rdb_hdr.user2');
        w{lc}.specnuc = h.rdb_hdr.user2;
    end
    w{lc}.slthick = h.image.slthick;
    w{lc}.fov = h.rdb_hdr.fov;
    w{lc}.res = h.image.dim_X;
    w{lc}.etl = h.rdb_hdr.etl;
    if isfield(h.image,'spokesPerSeg')
        w{lc}.spokesPerSeg = h.image.spokesPerSeg;
    else
        warning('h.image.spokesPerSeg not existing; setting to 0');
        w{lc}.spokesPerSeg = 0;
    end
    w{lc}.tr   = h.image.tr*1d-3;
    w{lc}.te   = h.rdb_hdr.te*1d-3;
    w{lc}.flip = h.image.mr_flip;
    w{lc}.bw   = h.image.vbw;
    sct = round(h.image.sctime*1d-6);
    w{lc}.sctime = sprintf('%d:%02d',floor(sct/60),round(rem(sct,60)));
    w{lc}.user36 = h.rdb_hdr.user36*1d-3;
    w{lc}.se_desc = h.series.se_desc;
    w{lc}.nex = h.image.nex;
    w{lc}.TG  = h.rdb_hdr.ps_mps_tg;
    w{lc}.R1 = h.rdb_hdr.ps_mps_r1;
    w{lc}.R2 = h.rdb_hdr.ps_mps_r2;
    if isfield(h.rdb_hdr,'dacq_data_type')
        w{lc}.data_type = h.rdb_hdr.dacq_data_type;
    else
        w{lc}.data_type = '';
    end
    for l=0:30
        % h.image = opuser
        % h.rdb_hdr = rhuser
        % w{lc}.(['user' num2str(l)]) = h.image.(['user' num2str(l)]);
        w{lc}.(['user' num2str(l)]) = h.rdb_hdr.(['user' num2str(l)]);
    end
    if ~isempty(regexpi(h.image.psd_iname,'fidall', 'once'))
        w{lc}.user3  = h.image.user3;
        w{lc}.user14 = h.image.user14;
        w{lc}.user20 = h.image.user20;
        w{lc}.user21 = h.image.user21;
    end
    for lf=1:length(flds)
        if ~ischar(w{lc}.(flds{lf}))
            w{lc}.(flds{lf}) = sprintf('%g',w{lc}.(flds{lf}));
        end
    end

    for lf=1:length(flds)
        lstr = length(w{lc}.(flds{lf}));
        if isfield(len,flds{lf})
            if (len.(flds{lf}) < lstr)
                len.(flds{lf}) = lstr;
            end
        else
            len.(flds{lf}) = lstr;
        end
    end
end


%% assemble string
ex_no_old = 0;
cname_old = '';
for lc=1:length(hc)
    str = '';
    ex_no = w{lc}.ex_no;
    if ex_no ~= ex_no_old
        str = sprintf('Exam number:  %d\n',w{lc}.ex_no);
        str = sprintf('%sScan date:    %s\n',str,w{lc}.scan_date);
        str = sprintf('%sPatient ID:   %s\n',str,w{lc}.patid);
        str = sprintf('%sPatient name: %s\n',str,w{lc}.patname);
        ex_no_old = ex_no;
        cname_old = '';
    end
    if ~strcmpi(cname_old,w{lc}.cname)
        str = sprintf('%sCoil name:    %s\n',str,w{lc}.cname);
        cname_old = w{lc}.cname;
    end

    ex_str{lc} = str;
    str = '';
    for lf=1:length(flds)
        ss = w{lc}.(flds{lf});
        if ~ischar(ss), warning('~ischar(ss)'); end
        str = sprintf('%s%s  ',str,pad(ss,len.(flds{lf})));
    end
    se_str{lc} = str;
end

end      % sub_extractinfo.m
