function dcm2readme(exam,sv,dict)
%DCM2README Extract info from dicom header for README file
% dcm2readme(exam,sv,dict)
%  exam/dir   Exam number to extract info from database (on MRI only)  ([])
%             Directory name
%        sv   Write output to README_hdr.txt file                    (true)
%      dict   Dicom dictionary                      ('gems-dicom-dict.txt')
%
%  2/2025 Rolf Schulte
% if (nargin<1), help(mfilename); return; end

if ~exist('exam','var'), exam = []; end
if ~exist('sv','var'),   sv = []; end
if isempty(sv),          sv = true; end
if ~exist('dict','var'), dict = []; end
if isempty(dict),        dict = 'gems-dicom-dict.txt'; end


%% checks
if ~exist(dict,'file')
    dcmdict = getenv('DCMDICT');
    if exist(dcmdict,'file')
        dict = dcmdict;
    else
        warning('dicom dictionary (%s) not found; using default',dict);
        dict = 'dicom-dict.mat';
    end
end
fprintf('Using dicom dictiory ''%s''\n',dict);


%% check filename first
if sv
    fname_rdm = [pwd filesep 'README_dcm.txt'];
    fname_old = [pwd filesep 'README_dcm_old.txt'];
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
flds = {'series','se_desc','runnum','slthick','fov','tr','te',...
    'flip','sctime','specnuc','nex','TG','R1','R2','data_type','psdname'};


%% extract info directly from scanner databse
dname = '';
if ~isempty(exam)
    if ischar(exam)
        if exist(exam,'dir')
            dname = exam;
        else
            warning('directory not found (''%s'')',exam);
        end
    else
        if isdeployed
            [status,cmdout] = system(['pathExtract ' num2str(exam)]);
            if status~=0, warning('status(=%g)~=0',status); end
            if ~isempty(regexpi(cmdout,'NOT FOUND',  'once'))
                warning('No exam found: ignoring');
            else
                dname = regexpi(cmdout,...
                    '/export/home1/sdc_image_pool/images/p\d*/e\d*/',...
                    'once','match');
                if ~exist(dname,'dir')
                    warning('Problem extracting dname(''%s''): ignoring',...
                        dname);
                    dname = '';
                end
            end
        else
            warning('exam(=%g): feature only available on MRI: ignoring',exam);
        end
    end
end
if ~isempty(dname)
    cwd = pwd;
    cd(dname);
else
    cwd = '';
end


%% read header
dcm = [dir('**/*.dcm') ; dir('**/*.dicom') dir('**/*.MRDC.*')];

ndcm = length(dcm);
inc = 0;
fprintf('Reading header for %d dicom files\n',ndcm);
for l=1:ndcm
    htmp = [];
    fname = [dcm(l).folder filesep dcm(l).name];
    if exist(fname,'file')
        fprintf('%d/%d: %s\n',l,ndcm,fname);
        htmp = dicominfo(fname,'dictionary',dict);
    else
        warning('strange input'); disp(fname);
    end
    if ~isempty(htmp)
        inc = inc+1;
        h{inc} = htmp;
    end
end
if inc==0, error('no suitable files found'); end



%% exclude double files
for l=1:length(h), tt(l) = str2num(h{l}.SeriesTime); end
[~,i2] = unique(tt);
if length(tt)~=length(i2)
    fprintf('Dicom not unique; excluding duplicates\n');
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

if ~isempty(cwd), cd(cwd); end

end      % dcm2readme.m


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
    % fprintf('%d/%d\n',lc,length(hc));
    h = hc{lc};
    % exam info
    w{lc}.ex_no   = str2double(h.StudyID);
    w{lc}.scan_date = h.StudyDate;
    w{lc}.patid   = h.PatientID;
    w{lc}.patname = h.PatientsName.FamilyName;
    w{lc}.cname   = h.ReceivingCoil;

    % series info
    w{lc}.runnum  = h.RawDataRunNumber;
    w{lc}.series  = h.SeriesNumber;
    % w{lc}.psd_iname = h.InternalPulseSeqName;
    % w{lc}.psdname   = h.PulseSeqName;
    w{lc}.psdname   = h.InternalPulseSeqName;
    w{lc}.specnuc = h.ImagedNucleus;
    w{lc}.slthick = h.SliceThickness;
    w{lc}.fov = h.DisplayFieldOfView;
    w{lc}.res = h.ImageDimensionX;
    w{lc}.etl = h.EchoTrainLength;
    w{lc}.tr   = h.RepetitionTime*1d-3;
    w{lc}.te   = h.EchoTime*1d-3;
    w{lc}.flip = h.FlipAngle;
    % w{lc}.bw   = h.image.vbw;
    sct = round(h.AcquisitionDuration*1d-6);
    w{lc}.sctime = sprintf('%d:%02d',floor(sct/60),round(rem(sct,60)));
    % w{lc}.user36 = h.rdb_hdr.user36*1d-3;
    w{lc}.se_desc = h.SeriesDescription;
    if isfield(h,'NumberOfExcitations')
        w{lc}.nex = h.NumberOfExcitations;
    else
        w{lc}.nex = [];
    end
    if isfield(h,'TransmissionGain')
        w{lc}.TG  = h.TransmissionGain;
    else
        w{lc}.TG  = [];
    end
    if isfield(h,'ActualReceiveGainAnalog')
        w{lc}.R1 = h.ActualReceiveGainAnalog;
    else
        w{lc}.R1 = [];
    end
    if isfield(h,'ActualReceiveGainDigital')
        w{lc}.R2 = h.ActualReceiveGainDigital;
    else
        w{lc}.R2 = [];
    end
    if isfield(h,'RawDataType_ImageType')
        w{lc}.data_type = h.RawDataType_ImageType;
    else
        w{lc}.data_type = [];
    end

    % for l=0:24
        % w{lc}.(['user' num2str(l)]) = h.(['UserData' num2str(l)]);
    % end

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
