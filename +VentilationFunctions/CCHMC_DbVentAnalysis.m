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
PatientID = 'test';
SequenceType = 'GRE';
ReconType = 'online';
VoxelSize = '[1,1,15]';
SliceOrientation = 'coronal';  % 'coronal' ||'transversal' || 'sagittal' ||'isotropic'
% vent_file = 'D:\OneDrive - cchmc\Lab\Random Subject analysis\CTC003-020-01\Xe\XeimageData.nii.gz';
vent_file = 'D:\OneDrive - cchmc\Lab\Random Subject analysis\CTC003-020-01\XeDicom\sub-020_ses-01_acq-01temp.dcm';
anat_file = ''; %D:\OneDrive - cchmc\Lab\Random Subject analysis\CTC003-020-01\H\HimageRegistered.nii.gz';
analysisFolder = 'D:\OneDrive - cchmc\Lab\Random Subject analysis\CTC003-020-01';
mask_file_name = 'D:\OneDrive - cchmc\Lab\Random Subject analysis\CTC003-020-01\mask\automask.nii.gz';
% mask_file_name = '';
print_output = 1; % 0 for stop printing 

[~,~,ext] = fileparts(vent_file);
if contains(ext,'.dcm') % only for dicom
    if isempty(PatientID) || isempty(SequenceType) || isempty(ReconType) || isempty(VoxelSize)
        VDP_inputs = LoadData.load_VDP_inputs(vent_file);
        if isempty(PatientID)
            PatientID = VDP_inputs.PatientID;
        end
        if isempty(SequenceType)
            SequenceType = VDP_inputs.SequenceType;
        end
        if isempty(ReconType)
            ReconType = VDP_inputs.ReconType;
        end
        if isempty(VoxelSize)
            VoxelSize = VDP_inputs.VoxelSize;
        end
    end
end 

if print_output == 0
    output = evalc(['VentilationFunctions.batchrun_vdp(PatientID, ' ...
        'SequenceType, ReconType, VoxelSize, vent_file, anat_file, analysisFolder, mask_file_name)']); % Captures and suppresses output
else
    VentilationFunctions.batchrun_vdp(PatientID,...
                                        SequenceType,...
                                        ReconType,...
                                        VoxelSize,...
                                        SliceOrientation,...
                                        vent_file,...
                                        anat_file,...
                                        analysisFolder,...
                                        mask_file_name);
end


    %% 

