function [data,header] = read_dicom(fname,dict,do_sort,read_all_hdr,verb)
%READ_DICOM  Read dicom data and header from single image or series
%[data,header] = read_dicom(fname,dict,do_sort,read_all_hdr,verb)
%        fname   File name 
%                single file -> read only single dicom
%                directory name -> read all dicoms in that directory
%                (default=current directory)
%         dict   Dicom dictionary                   ('gems-dicom-dict.txt')
%      do_sort   Read in dicom header and sort according to             (0) 
%                (1) .InstanceNumber; (2) .SliceLocation
%                Otherwise: sorting via dir or suffix number (if ending with num)
% read_all_hdr   Read headers from all images                       (false)
%         verb   Verbose mode                                        (true)
%
%         data   Dicom data
%       header   Dicom header
%
% 8/2020  Rolf Schulte
% See also DICOMREAD, DICOMINFO.


%% default input
if ~exist('fname','var'),  fname = []; end
if isempty(fname),         fname = '.'; end
if ~exist('dict','var'),   dict = []; end
if isempty(dict),          dict = 'gems-dicom-dict.txt'; end
if ~exist('do_sort','var'),do_sort = []; end
if isempty(do_sort),       do_sort = 0; end
if ~exist('read_all_hdr','var'),read_all_hdr = []; end
if isempty(read_all_hdr),  read_all_hdr = false; end
if ~exist('verb','var'),   verb = []; end
if isempty(verb),          verb = true; end


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
if ~exist(fname,'file'), error('file ''%s'' not found',fname); end
if verb, fprintf('Using dicom dictiory ''%s''\n',dict); end


%% reading single dicom file (if fname == file name) exit
if ~isdir(fname)
    if verb, fprintf('Reading single dicom file ''%s''\n',fname); end
    if ~isdicom(fname), warning('fname(=%s) not dicom',fname); end
    data = double(dicomread(fname));
    if nargout>1, header = dicominfo(fname,'dictionary',dict); end
    return
end


%% reading dicom series
% extract + check file names
fns = dir;                        % list all files in directory
ii_dcm = [];                      % index for including files
ii_sort = [];                     % index for sorting files
for l=1:length(fns)
    if ~fns(l).isdir              % exclude directory names
        if isdicom(fns(l).name)   % check if file names are dicoms
            ii_dcm = [ii_dcm, l]; % include file
            
            % sorting via file name suffix .<number>
            tmp = regexp(fns(l).name,'\.\d+$');
            if ~isempty(tmp)
                ii_sort = [ii_sort, str2double(fns(l).name(tmp+1:end))]; 
            end
        else
            fprintf('Warning: %s not a dicom; skipping\n',fns(l).name);
        end
    end
end
fns = fns(ii_dcm);

ii_sort = ii_sort(~isnan(ii_sort));    % exclude NaN's
if length(ii_dcm)==length(ii_sort)     % sort if all files end w/ number
    fns(ii_sort) = fns;                % pure file name list
    if verb, fprintf('Sorting dicom series via suffix\n'); end
end
nf = length(fns);                      % #dicom files to read
if verb, fprintf('Dicom files to read: nf = %g\n',nf); end


%% read dicom data
% assumes all data has same size
if verb, fprintf('fns(1).name = %s\n',fns(1).name); end
tmp = dicomread(fns(1).name);          % read single file -> get data size
data = single(zeros(size(tmp,1),size(tmp,2),nf));  % preallocate as single
data(:,:,1) = tmp;
for l=2:nf 
    if verb, fprintf('fns(%g).name = %s\n',l,fns(l).name); end
    data(:,:,l) = dicomread(fns(l).name); 
end


%% read dicom header
if (do_sort>0) || read_all_hdr
    nn = nf;                      % read all headers
else
    nn = 1;                       % read only 1 header (saves time)
end
for l=1:nn
    if l==2, warning('off'); end  % avoid multiple warnings
    header(l) = dicominfo(fns(l).name,'dictionary',dict); 
end
if do_sort>0, warning('on'); end


%% sort data
switch do_sort    % sort if not done so already by suffix
    case 0
    case 1
        if verb, fprintf('Sorting dicom series via header.InstanceNumber\n'); end
        if isfield(header(1),'InstanceNumber')
            i_in = zeros(1,nf);
            i_sn = zeros(1,nf);
            for l=1:nf
                i_in(1,l) = header(l).InstanceNumber;
                i_sn(1,l) = header(l).SeriesNumber;
            end
        else
            warning('field missing; skipping sorting');
        end
        if any(diff(sort(i_in))==0)
            warning('same InstanceNumber multiple times; data loss');
        end
        if any(diff(i_sn))
            warning('different SeriesNumber');
            fprintf('sort series into different directories first\n');
        end
        if any(diff(sort(i_in))>1)
            warning('delta(InstanceNumber)>1; compressing');
            [tmp,i_in] = sort(i_in);
        end
        if min(i_in)~=1
            warning('min(i_in)(=%g)~=1; subtracting min',min(i_in));
            i_in = i_in-min(i_in);
        end
        data(:,:,i_in) = data;
    case 2
        if verb, fprintf('Sorting dicom series via header.SliceLocation\n'); end
        if isfield(header(1),'SliceLocation')
            sl_loc = zeros(1,nf);
            i_sn = zeros(1,nf);
            for l=1:nf
                sl_loc(1,l) = header(l).SliceLocation;
                i_sn(1,l) = header(l).SeriesNumber;
            end
        else
            warning('field missing; skipping sorting');
        end
        if any(abs(diff(sort(sl_loc)))<1d-10)
            warning('same SliceLocation multiple times; data loss');
        end
        if any(diff(i_sn))
            warning('different SeriesNumber');
            fprintf('sort series into different directories first\n');
        end
        [tmp,i_in] = sort(sl_loc);
        data(:,:,i_in) = data;
    otherwise
end


end   % main function read_dicom.m
