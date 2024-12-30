function [convex_boundary, concave_boundary, region_inbetween_lungs, mask_spline] = split_lung_boundaries(mask_struct, slice_range, binary_flag)
% SPLIT_LUNG_BOUNDARIES separates the lung boundaries into convex and concave portions.
%
% Input:
%   - mask_struct: structure containing binary lung masks.
%   - slice_range: list of slices for boundary splitting.
%   - binary_flag: flag to control binary or labeled output.
%
% Output:
%   - convex_boundary: 3D matrix containing convex ("outer") boundaries.
%   - concave_boundary: 3D matrix containing concave ("inner") boundaries.
%   - region_inbetween_lungs: 3D mask of the region between lungs.
%   - mask_spline: 3D matrix of combined convex and concave regions.
%
% W. Zha, modified with fixes on [current_date]

% Default inputs
lungsmask = mask_struct.lungsmask;
if nargin < 2
    slice_range = 1:size(lungsmask, 3);
    binary_flag = true;
elseif nargin < 3
    binary_flag = true;
end

% Initialize outputs
[nRows, nCols, nSlices] = size(lungsmask);
convex_boundary = zeros(nRows, nCols, nSlices);
concave_boundary = zeros(nRows, nCols, nSlices);
region_inbetween_lungs = zeros(nRows, nCols, nSlices);
mask_spline = zeros(nRows, nCols, nSlices);

% Identify lobes and initialize lobe masks
[mask_struct, lobe_label_values] = VentilationFunctions.Akmeans.identify_lobes(mask_struct);
lobemask = initialize_lobes(mask_struct, lobe_label_values, lungsmask);

% Define lobe names and values
lobe_list = {'RUL', 'LUL'};
lobe_values = [2, 5];

% Process each slice in the slice range
for n_ind = 1:length(slice_range)
    n_sl = slice_range(n_ind);
    slname = sprintf('sl%d', n_sl);
    concave_pts.RUL = [];
    concave_pts.LUL = [];

    % Process each lobe
    for n_lobe = 1:numel(lobe_list)
        lobe_name = lobe_list{n_lobe};
        this_lobe = lobemask.(lobe_name);
        slice_mask = this_lobe(:, :, n_sl);

        if sum(slice_mask(:)) > 0
            % Label connected components
            [labeled_mask, nObj] = bwlabel(slice_mask);
            for n = 1:nObj
                % Extract object and contours
                object_mask = labeled_mask == n;
                lobe_contours = VentilationFunctions.Akmeans.mask2poly(object_mask, 'Exact', 'CW');

                % Remove invalid contour points
                invalids = lobe_contours(:, 1) < 1 | lobe_contours(:, 2) < 1 | ...
                           lobe_contours(:, 1) > nCols | lobe_contours(:, 2) > nRows;
                lobe_contours(invalids, :) = [];

                if isempty(lobe_contours)
                    warning('Empty contours for lobe "%s" on slice %d.', lobe_name, n_sl);
                    continue;
                end

                % Split into convex and concave regions
                [convex_list, concave_list] = VentilationFunctions.Akmeans.find_conv_region(lobe_contours);

                % Convex boundary
                convexBoundary = lobe_contours(convex_list, :);
                valid_convex = validate_indices(convexBoundary, nRows, nCols);
                convexBoundary = convexBoundary(valid_convex, :);
                if ~isempty(convexBoundary)
                    convex_indices = sub2ind([nRows, nCols], convexBoundary(:, 2), convexBoundary(:, 1));
                    convex_edge = zeros(nRows, nCols);
                    convex_edge(convex_indices) = 1;
                    if binary_flag
                        convex_boundary(:, :, n_sl) = convex_boundary(:, :, n_sl) + convex_edge;
                    else
                        convex_boundary(:, :, n_sl) = convex_boundary(:, :, n_sl) + convex_edge * lobe_values(n_lobe);
                    end
                end

                % Concave boundary
                concaveBoundary = lobe_contours(concave_list, :);
                valid_concave = validate_indices(concaveBoundary, nRows, nCols);
                concaveBoundary = concaveBoundary(valid_concave, :);
                if ~isempty(concaveBoundary)
                    concave_indices = sub2ind([nRows, nCols], concaveBoundary(:, 2), concaveBoundary(:, 1));
                    concave_edge = zeros(nRows, nCols);
                    concave_edge(concave_indices) = 1;
                    concave_pts.(lobe_name).(slname) = concaveBoundary;
                    if binary_flag
                        concave_boundary(:, :, n_sl) = concave_boundary(:, :, n_sl) + concave_edge;
                    else
                        concave_boundary(:, :, n_sl) = concave_boundary(:, :, n_sl) + concave_edge * lobe_values(n_lobe);
                    end
                end
            end
        end
    end

    % Combine concave regions between lungs
    if isfield(concave_pts.RUL, slname) && isfield(concave_pts.LUL, slname)
        region_inbetween_lungs(:, :, n_sl) = poly2mask([concave_pts.RUL.(slname)(:, 1); concave_pts.LUL.(slname)(end:-1:1, 1)], ...
                                                       [concave_pts.RUL.(slname)(:, 2); concave_pts.LUL.(slname)(end:-1:1, 2)], ...
                                                       nRows, nCols);
    end
end
end

%% Helper Function: Initialize Lobes
function lobemask = initialize_lobes(mask_struct, lobe_label_values, lungsmask)
[nRows, nCols, nSlices] = size(lungsmask);
lobemask = struct();

if length(lobe_label_values) == 5
    % Initialize masks for RLL, RML, and LLL
    se = strel('disk', 2);
    RLL = imdilate(mask_struct.lobemaskRLL, se);
    RML = imdilate(mask_struct.lobemaskRML, se);
    LLL = imdilate(mask_struct.lobemaskLLL, se);
    lobemask.RUL = (mask_struct.lobemaskRUL + RML + RLL) .* lungsmask;
    lobemask.RUL(lobemask.RUL > 0) = 1;
    lobemask.LUL = (mask_struct.lobemaskLUL + LLL) .* lungsmask;
    lobemask.LUL(lobemask.LUL > 0) = 1;
elseif length(lobe_label_values) == 2
    lobemask.RUL = mask_struct.lobemaskRUL;
    lobemask.LUL = mask_struct.lobemaskLUL;
end
end

%% Helper Function: Validate Indices
function valid = validate_indices(contours, nRows, nCols)
valid = contours(:, 1) >= 1 & contours(:, 1) <= nCols & ...
        contours(:, 2) >= 1 & contours(:, 2) <= nRows;
end
