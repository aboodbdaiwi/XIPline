function idx = pickSlices(mask3d, plane, nKeep, margin)
% mask3d : 3D logical/numeric mask [X Y Z]
% plane  : 'ax' (axial, vary Z), 'co' (coronal, vary Y), 'sa' (sagittal, vary X)
% nKeep  : number of evenly spaced slices to keep (default 7)
% margin : exclude this many slices at each end (default 5)

    if nargin < 3 || isempty(nKeep),  nKeep  = 7; end
    if nargin < 4 || isempty(margin), margin = 5; end

    switch lower(plane)
        case 'ax'  % axial → vary Z; collapse X and Y
            nz = find(squeeze(any(any(mask3d, 1), 2)));   % 1x1xZ -> Z-vector
        case 'co'  % coronal → vary Y; collapse X and Z
            nz = find(squeeze(any(any(mask3d, 1), 3)));   % 1xYx1 -> Y-vector
        case 'sa'  % sagittal → vary X; collapse Y and Z
            nz = find(squeeze(any(any(mask3d, 2), 3)));   % Xx1x1 -> X-vector
        otherwise
            error('plane must be ''ax'', ''co'', or ''sa''.');
    end

    % Exclude first/last 'margin' slices if enough exist
    if numel(nz) > 2*margin
        nz = nz((margin+1):(end-margin));
    end

    % Choose up to nKeep evenly spaced slices
    if numel(nz) > nKeep
        idx = round(linspace(nz(1), nz(end), nKeep));
    else
        idx = nz(:).'; % keep all remaining
    end
end