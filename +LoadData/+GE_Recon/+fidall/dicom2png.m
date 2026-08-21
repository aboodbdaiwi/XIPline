function dicom2png(fname,verb)
%DICOM2PNG  Convert dicom data into png files
%
%  dicom2png(fname,verb)
%        fname  File name 
%               single file -> read only single dicom
%               directory name -> read all dicoms in that directory
%               (default=current directory)
%         verb  >0 verbose mode; >1 plotting
%
% 9/2016  Rolf Schulte
% See also DICOMREAD, DICOMINFO.


%% default input + checks
if ~exist('fname','var'),  fname = []; end
if isempty(fname),         fname = '.'; end
if ~exist('verb','var'),   verb = []; end
if isempty(verb),          verb = 2; end

if ~exist(fname,'file'), error('file ''%s'' not found',fname); end
if verb>1, clf; colormap gray; end


%% reading single dicom file (if fname == file name) exit
if ~isdir(fname),
    if verb, fprintf('Reading single dicom file ''%s''\n',fname); end
    if ~isdicom(fname), warning('fname(=%s) not dicom',fname); end
    read_dicom_save_png(fname,verb);
    return
end


%% reading dicom series
% extract + check file names
fns = dir;                        % list all files in directory
for l=1:length(fns),
    if ~fns(l).isdir              % exclude directory names
        if isdicom(fns(l).name),  % check if file names are dicoms
            if verb, fprintf('Reading dicom file ''%s''\n',fns(l).name); end
            read_dicom_save_png(fns(l).name,verb);
        else
            fprintf('Warning: %s not a dicom; skipping\n',fns(l).name);
        end
    end
end
end

%% subfunction
function read_dicom_save_png(fname,verb)
b = double(dicomread(fname));
h = dicominfo(fname); %,'dictionary',dict);
if isfield(h,'WindowCenter') && isfield(h,'WindowWidth')
    sca = h.WindowCenter + h.WindowWidth*[-1 1]/2;
else
    fprintf('Warning: header fields WindowCenter and/or WindowWidth missing\n');
    sca = [min(b(:)) max(b(:))];
end
if verb>1, imagesc(b,sca); axis image off; drawnow; end
if ~isempty(regexpi(fname,'\.dcm$')), fname = fname(1:end-4); end
% saveas(gcf,[fname '.png'],'png');

b = b/sca(2);  % only scale to max
imwrite(b,[fname '.png'],'png');
end
