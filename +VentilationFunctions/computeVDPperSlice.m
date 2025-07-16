function [vdp_per_slice_local, vdp_per_slice_global] = computeVDPperSlice(defectMap)
% computeVDPperSlice calculates slice-wise ventilation defect percentage (VDP)
% using both local (slice) and global (volume-wide) lung volume.
%
% Inputs:
%   defectMap - 3D binary mask:
%               0 = background
%               1 = ventilated lung
%               2 = defect
%
% Outputs:
%   vdp_per_slice_local  - VDP per slice, using only the lung area in that slice
%   vdp_per_slice_global - VDP per slice, using total lung volume across all slices

    % Get dimensions
    [~, ~, numSlices] = size(defectMap);
    vdp_per_slice_local = zeros(numSlices, 1);
    vdp_per_slice_global = zeros(numSlices, 1);

    % Total lung mask across all slices (lung + defect)
    global_lung_mask = (defectMap == 1 | defectMap == 2);
    total_lung_volume = sum(global_lung_mask(:));

    for sl = 1:numSlices
        slice = defectMap(:, :, sl);

        % Local lung + defect volume
        local_lung_volume = sum(slice(:) == 1) + sum(slice(:) == 2);
        local_defect_volume = sum(slice(:) == 2);

        % Local VDP (per-slice lung volume)
        if local_lung_volume > 0
            vdp_per_slice_local(sl) = 100 * local_defect_volume / local_lung_volume;
        else
            vdp_per_slice_local(sl) = NaN;
        end

        % Global VDP (using total lung volume)
        vdp_per_slice_global(sl) = 100 * local_defect_volume / total_lung_volume;
    end
end
