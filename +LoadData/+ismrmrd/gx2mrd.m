function gx2mrd(Subj_ID,xe_file,cal_file,ute_file)

% Get data files if needed
if nargin < 1
    [xe_file,cal_file,ute_file] = uigetgxfiles();
    Subj_ID = inputdlg('Subject ID','Input Subject ID',[1 50]); % Or maybe better to prompt user to type in Subject ID?
    Subj_ID = Subj_ID{1};
elseif nargin < 4
    [xe_file,cal_file,ute_file] = uigetgxfiles();
end
%%
% get file extension and path from the xenon file - For GE, there are 3
% files, so xe-file will be a cell
if ~iscell(xe_file)
    [xe_path,~,fileext] = fileparts(xe_file);
else
    [xe_path,~,fileext] = fileparts(xe_file{1});
end

% Send to platform-specific MRD converters
if strcmp(fileext,'.dat')
    %Try both John's xpdixon as well as the normal way
    try
        Siemens.xpdixon_2303_2_mrd(xe_file,Subj_ID,fullfile(xe_path,[Subj_ID '_dixon.h5']));
    catch
        Siemens.dissolved_to_ismrmrd(xe_file,Subj_ID,fullfile(xe_path,[Subj_ID '_dixon.h5']));
    end
    Siemens.calibration_to_ismrmrd(cal_file,Subj_ID,fullfile(xe_path,[Subj_ID '_calibration.h5']));
    Siemens.ute_to_ismrmrd(ute_file,Subj_ID,fullfile(xe_path,[Subj_ID '_proton.h5']));
    
else
    GE.dissolved_to_ismrmrd(xe_file,Subj_ID,fullfile(xe_path,[Subj_ID '_dixon.h5']));
    GE.calibration_to_ismrmrd(cal_file,Subj_ID,fullfile(xe_path,[Subj_ID '_calibration.h5']))
    GE.ute_to_ismrmrd(ute_file,Subj_ID,fullfile(xe_path,[Subj_ID '_proton.h5']));
end