function frac_vent = computeFracVent(image, mask_including_trachea, bag_volume, voxelSize_mm)
% computeFracVent  Fractional ventilation (unitless) from HP 129Xe magnitude image.
%
%   frac_vent = computeFracVent(image, bag_volume, mask_including_trachea)
%
% Inputs
%   image                  : numeric array (magnitude image; any real/complex)
%   bag_volume             : positive scalar, liters (L)
%   mask_including_trachea : logical/numeric mask (lungs + trachea), same size as image
%
% Output
%   frac_vent : double array, same size as image (unitless fractional ventilation)
%
% Notes
% - Voxel volume assumed isotropic with side = 0.3125 cm (i.e., 3.125 mm).
% - DEAD_SPACE_ML fixed at 10 mL.
% - Attempts to save NIfTI as 'frac_vent_output.nii' if niftiwrite is available.

    %---------------------- Input checks ----------------------%
    narginchk(4, 4);

    if ~isscalar(bag_volume) || ~isnumeric(bag_volume) || ~(bag_volume > 0) 
        error('FRAC_VENT requires a positive numeric bag_volume in liters.');
    end

    if nargin < 4 || isempty(voxelSize_mm)
        error('FRAC_VENT requires voxel Size.');
    end
    if nargin < 3 || isempty(bag_volume)
        error('FRAC_VENT requires bag volume.');
    end    
    if nargin < 2 || isempty(mask_including_trachea)
        error('FRAC_VENT requires mask_including_trachea (lungs + trachea).');
    end
    if ~isequal(size(mask_including_trachea), size(image))
        error('mask_including_trachea must have the same size as image.');
    end

    if ~any(mask_including_trachea(:))
        error('mask_including_trachea is empty.');
    end

    %---------------------- Constants -------------------------%
    % voxel_side_cm = 0.3125;                 % cm
    % voxel_vol_ml  = voxel_side_cm^3;        % cm^3 == mL

    % voxelSize_mm = [3 3 15];          
    % voxel_vol_ml = prod(voxelSize_mm) / 1000;  % mm3 → mL
    voxel_vol_ml = prod(voxelSize_mm);  % mm3 → mL

    DEAD_SPACE_ML = 10.0;                   % mL

    tcv_gas_vol_ml = bag_volume - DEAD_SPACE_ML; % total delivered minus dead space
    if ~(tcv_gas_vol_ml > 0)
        error('Computed TCV gas volume <= 0 mL; check bag_volume / dead-space.');
    end

    %---------------------- Computation -----------------------%
    mag = abs(double(image));                         % ensure double; handle complex
    big = logical(mask_including_trachea);            % ensure logical
    signal_total = sum(mag(big), 'all');

    if ~(signal_total > 0)
        error('Total signal in mask_including_trachea is zero; cannot compute FV.');
    end

    sig_to_vol = tcv_gas_vol_ml / signal_total;       % mL per signal unit
    frac_vent  = (mag * sig_to_vol) / voxel_vol_ml;   % unitless FV

    %---------------------- Optional NIfTI save ---------------%
    try
        if exist('niftiwrite', 'file') == 2
            % Create output directory if needed
            if ~exist('tmp', 'dir'); mkdir('tmp'); end
            % Minimal metadata; adjust if you have spatial info to embed
            niftiwrite(frac_vent, fullfile('tmp', 'frac_vent_output.nii'));
        end
    catch %#ok<CTCH>
        % Silent if writing fails (matches Python try/except behavior)
    end
end
