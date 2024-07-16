function dicom2mat(varargin)
%DICOM2MAT Read in dicom file and save as matlab .mat
% Input cell array, string array or multiple file names
% Use dicom2mat(ls) for whole directory
if nargin<1, help(mfilename); return; end

for l1=1:nargin
    fname = varargin{l1};
    if iscell(fname),
        for l2=1:length(fname), dicom2mat_sub(fname{l2}); end
    else
        for l2=1:size(fname,1), dicom2mat_sub(fname(l2,:)); end
    end
end


end

function dicom2mat_sub(fname)
if exist(fname)==2,
    fprintf('Converting ''%s''\n',fname);
    b1 = dicomread(fname);
    b1h = dicominfo(fname);
    save([fname '.mat'],'b1','b1h');
else
    warning('''%s'' not a file',fname);
end

end
