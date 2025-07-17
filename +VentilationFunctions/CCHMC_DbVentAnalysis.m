% ----------129Xe Ventilation Analysis (Database version)-----------
%   Inputs: 
%          
%   Outputs:
%         
%   Example: 
%   Package: 
%
%   Author: Abdullah Bdaiwi
%   Work email: abdullah.bdaiwi@cchmc.org
%   Personal email: abdaiwi89@gmail.com
%   Website: https://www.cincinnatichildrens.org/research/divisions/c/cpir
% ------------------------------------------------------------------------
clc; clear; close;
% inputs
MainInput.SubjectID = 'IRC186H-0468';
MainInput.Age = 11;
MainInput.Sex = 'F';
MainInput.Disease = 'BPD';
MainInput.ScanDate = '20230120';
MainInput.Scanner = 'Philips';
MainInput.ScannerSoftware = '5.9.0';
MainInput.SequenceType = '2DGRE';
MainInput.ReconType = 'Online';
MainInput.denoiseXe = 'no';
MainInput.Analyst = 'Database';
MainInput.VoxelSize = '[3,3,15]';
MainInput.SliceOrientation = 'transversal'; % 'coronal' ||'transversal' || 'sagittal' ||'isotropic'
MainInput.vent_file = 'D:\OneDrive - cchmc\Lab\Random Subject analysis\sub-186H_0468\ses-20230120\vent\sub-0468_ses-20230120_acq-2.dcm';
MainInput.anat_file = 'D:\OneDrive - cchmc\Lab\Random Subject analysis\sub-186H_0468\ses-20230120\anatvent\sub-0468_ses-20230120_acq-1.dcm';
MainInput.analysisFolder = 'D:\OneDrive - cchmc\Lab\Random Subject analysis\sub-186H_0468\ses-20230120\analysis\acq_1_rec_2';
MainInput.mask_file_name = '';
print_output = 1; % 0 for stop printing 


% PatientID = 'test';
% SequenceType = 'GRE';
% ReconType = 'online';
% VoxelSize = '[3,3,15]';
% SliceOrientation = 'coronal';  
% vent_file = 'D:\OneDrive - cchmc\Lab\Random Subject analysis\CTC003-020-01\Xe\XeimageData.nii.gz';

% [~,~,ext] = fileparts(vent_file);
% if contains(ext,'.dcm') % only for dicom
%     if isempty(PatientID) || isempty(SequenceType) || isempty(ReconType) || isempty(VoxelSize)
%         VDP_inputs = LoadData.load_VDP_inputs(vent_file);
%         if isempty(PatientID)
%             PatientID = VDP_inputs.PatientID;
%         end
%         if isempty(SequenceType)
%             SequenceType = VDP_inputs.SequenceType;
%         end
%         if isempty(ReconType)
%             ReconType = VDP_inputs.ReconType;
%         end
%         if isempty(VoxelSize)
%             VoxelSize = VDP_inputs.VoxelSize;
%         end
%     end
% end 

if print_output == 0
    output = evalc(['VentilationFunctions.batchrun_vdp(MainInput)']); % Captures and suppresses output
else
    VentilationFunctions.CCHMC_Db_Vent_Pipeline(MainInput);
end


    %% 

