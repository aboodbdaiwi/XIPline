function [xe_file,cal_file,ute_file] = uigetgxfiles()

platform = questdlg('On what platform was your data acquired','Platform','Siemens','GE','Philips','Siemens');

switch platform
    case 'Siemens'
        if ~ispc 
            menu('Select Gas Exchange, then Calibration, then Anatomic Raw Data','OK'); 
        end
        [xe_file_a,xe_path] = uigetfile('*.dat','Select Gas Exchange Raw Data file');
        xe_file = fullfile(xe_path,xe_file_a);
        [cal_file_a,cal_path] = uigetfile('*.dat','Select Calibration Raw Data file');
        cal_file = fullfile(cal_path,cal_file_a);
        [ute_file_a,ute_path] = uigetfile('*.dat','Select Anatomic Raw Data file');
        ute_file = fullfile(ute_path,ute_file_a);
    case 'GE'
        [xe_file_a,xe_path] = uigetfile('*.h5','Select Gas Exchange Scan Archive file');
        xe_file{1} = fullfile(xe_path,xe_file_a);
        [xe_file_a,xe_path] = uigetfile('*.mat','Select Gas Exchange .mat file generated from prep_dixon_gx');
        xe_file{2} = fullfile(xe_path,xe_file_a);
        [xe_file_a,xe_path] = uigetfile('*.mat','Select Calibration .mat file generated from recon_calibration');
        xe_file{3} = fullfile(xe_path,xe_file_a);
        
        [cal_file_a,cal_path] = uigetfile('*.h5','Select Calibration Scan Archive file');
        cal_file = fullfile(cal_path,cal_file_a);

        [ute_file_a,ute_path] = uigetfile('*.h5','Select Anatomic Raw Data file');
        ute_file = fullfile(ute_path,ute_file_a);
    case 'Philips'
        [xe_file_a,xe_path] = uigetfile('*.raw','Select Gas Exchange .raw File');
        xe_file{1} = fullfile(xe_path,xe_file_a);
        [xe_file_a,xe_path] = uigetfile('*.data','Select Gas Exchange .data file generated from prep_dixon_gx');
        xe_file{2} = fullfile(xe_path,xe_file_a);
        
        [cal_file_a,cal_path] = uigetfile('*.raw','Select Calibration .raw file');
        cal_file{1} = fullfile(cal_path,cal_file_a);
        [cal_file_a,cal_path] = uigetfile('*.data','Select Calibration .data file');
        cal_file{2} = fullfile(cal_path,cal_file_a);

        [ute_file_a,ute_path] = uigetfile('*.raw','Select Anatomic .raw file');
        ute_file{1} = fullfile(ute_path,ute_file_a);
        [ute_file_a,ute_path] = uigetfile('*.data','Select Anatomic .data file');
        ute_file{2} = fullfile(ute_path,ute_file_a);
end