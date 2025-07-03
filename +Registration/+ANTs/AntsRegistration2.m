function [MainInput, Proton] = AntsRegistration2(MainInput, Proton, Ventilation, GasExchange)
% AntsRegistration - Registers proton image to fixed image using ANTs
% and optionally applies transform to a second image (e.g., lung mask).
%
% Inputs:
%   MainInput.AnalysisType = 'Ventilation' or 'GasExchange'
%   Proton.Image = proton image to be registered (moving)
%   Ventilation.Image = fixed image for registration
%   Ventilation.LungMask or GasExchange.LungMask (optional second image)
%
% Outputs:
%   MainInput.HRegisteredPath = path to registered proton image
%   Proton.ProtonRegistered = registered volume
%   Proton.ProtonRegisteredColored = overlay with fixed image

clc;
MainInput.transform_reg = false;  % set to true if applying to 2nd image

% Determine fixed and optional second moving image
moving1 = double(Proton.Image);
switch lower(MainInput.AnalysisType)
    case 'ventilation'
        fixed1 = try_get_field(Ventilation, 'Image', Ventilation);
        if MainInput.transform_reg
            moving2 = try_get_field(Ventilation, 'LungMask');
        end
    case 'gasexchange'
        fixed1 = try_get_field(GasExchange, 'VentImage', Ventilation);
        if MainInput.transform_reg
            moving2 = try_get_field(GasExchange, 'LungMask');
        end
    otherwise
        error('Unknown AnalysisType: %s', MainInput.AnalysisType);
end

% Resize if needed
if ~isequal(size(moving1), size(fixed1))
    moving1 = imresize3(moving1, size(fixed1));
end

% Prepare paths
ANTSPath = fileparts(mfilename('fullpath'));
tmp_path = fullfile(MainInput.OutputPath, 'tmp');
if ~exist(tmp_path, 'dir'), mkdir(tmp_path); end

% File paths
fixed_path = fullfile(tmp_path, 'image_static.nii');
moving1_path = fullfile(tmp_path, 'image_moving1.nii');
output_prefix = fullfile(tmp_path, 'thisTransform_');
moving1_out = fullfile(tmp_path, 'moving_reg.nii.gz');
transform_mat = [output_prefix, '0GenericAffine.mat'];
niftiwrite(fixed1, fixed_path);
niftiwrite(moving1, moving1_path);

% Build registration command
ants_reg = fullfile(ANTSPath, 'antsRegistration.exe');
cmd_register = sprintf(['"%s" --dimensionality 3 --float 0 --interpolation BSpline ' ...
    '--metric MI[%s,%s,1,32,Regular,1] --transform Affine[0.1] ' ...
    '--convergence [20x20x20,1e-6,20] --shrink-factors 4x2x1 ' ...
    '--smoothing-sigmas 0x0x0 --output [%s,%s] --verbose 1'], ...
    ants_reg, fixed_path, moving1_path, output_prefix, moving1_out);

% Execute
fprintf('Running ANTs registration...\n');
[status, cmdout] = system(cmd_register);
if status ~= 0
    error('ANTs registration failed:\n%s', cmdout);
end

% Read registered image
MainInput.HRegisteredPath = moving1_out;
try
    A1 = LoadData.load_nii(moving1_out);
catch
    A1 = LoadData.load_untouch_nii(moving1_out);
end
A = double(squeeze(A1.img));
A = imrotate(A, 90);  % adjust orientation
Proton.ProtonRegistered = A;

% Overlay visualization
fixedVolume = fixed1;
for s = 1:size(fixedVolume, 3)
    A = Proton.ProtonRegistered(:,:,s);
    B = fixedVolume(:,:,s);
    Proton.ProtonRegisteredColored(:,:,:,s) = imfuse(A,B,'falsecolor','ColorChannels','green-magenta');
end

% Optional: transform second image (e.g., lung mask)
if MainInput.transform_reg && exist('moving2', 'var')
    moving2_path = fullfile(tmp_path, 'image_moving2.nii');
    output2_path = fullfile(tmp_path, 'mask_reg.nii.gz');
    niftiwrite(double(moving2), moving2_path);
    
    ants_apply = fullfile(ANTSPath, 'antsApplyTransforms.exe');
    cmd_apply = sprintf(['"%s" -d 3 -e 0 -i %s -r %s -o %s -t %s'], ...
        ants_apply, moving2_path, fixed_path, output2_path, transform_mat);
    system(cmd_apply);
    fprintf('Second image transformed and saved: %s\n', output2_path);
end
end

function out = try_get_field(structure, field, default)
% Try to get a field from a struct, else return default
if isfield(structure, field)
    out = structure.(field);
else
    out = default;
end
end
