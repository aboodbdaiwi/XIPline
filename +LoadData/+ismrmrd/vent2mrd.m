function vent2mrd(Subj_ID,xe_file,anat_file)

if nargin < 1
    if ~ispc 
       menu('Select Ventilation, then Ventilation Anatomic Raw Data','OK'); 
    end
    [xe_file_a,xe_path] = uigetfile('*.dat','Select Ventilation Raw Data file');
    xe_file = fullfile(xe_path,xe_file_a);
    [anat_file_a,anat_path] = uigetfile('*.dat','Select Ventilation Anatomic Raw Data file');
    anat_file = fullfile(anat_path,anat_file_a);
    Subj_ID = inputdlg('Subject ID','Input Subject ID',[1 50]); % Or maybe better to prompt user to type in Subject ID?
    Subj_ID = Subj_ID{1};
elseif nargin < 2
    if ~ispc 
       menu('Select Ventilation, then Ventilation Anatomic Raw Data','OK'); 
    end
    [xe_file_a,xe_path] = uigetfile('*.dat','Select Ventilation Raw Data file');
    xe_file = fullfile(xe_path,xe_file_a);
    [anat_file_a,anat_path] = uigetfile('*.dat','Select Ventilation Anatomic Raw Data file');
    anat_file = fullfile(anat_path,anat_file_a);
end
    
[~,~,fileext] = fileparts(xe_file);

if strcmp(fileext,'.dat')
    Siemens.gre_to_ismrmrd(xe_file,Subj_ID,fullfile(xe_path,[Subj_ID '_vent.h5']));
    Siemens.gre_to_ismrmrd(anat_file,Subj_ID,fullfile(anat_path,[Subj_ID '_ventanat.h5']));
else
    GE.gre_to_ismrmrd(xe_file,Subj_ID,fullfile(xe_path,[Subj_ID '_vent.h5']));
    GE.gre_to_ismrmrd(anat_file,Subj_ID,fullfile(anat_path,[Subj_ID '_ventanat.h5']));
end