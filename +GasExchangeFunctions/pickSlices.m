function idx = pickSlices(mask3d, plane, nKeep, margin)
% mask3d : 3D logical/numeric mask [X Y Z]
% plane  : 'ax' (axial, vary Z), 'co' (coronal, vary Y), 'sa' (sagittal, vary X)
% nKeep  : number of evenly spaced slices to keep (default 7)
% margin : exclude this many slices at each end (default 5)

    if nargin < 3 || isempty(nKeep),  nKeep  = 7; end
    if nargin < 4 || isempty(margin), margin = 5; end
    
    switch lower(plane)
        case 'ax'  % axial → vary Y
            sliceIndx = zeros(1, size(mask3d,1));
            for i = 1:size(mask3d,1)
                mask_slice = squeeze(mask3d(i,:,:));
                sliceIndx(i) = i;
                if sum(mask_slice(:)) <= 25
                    mask3d(i,:,:) = 0;
                    sliceIndx(i) = 0;
                end
            end
            nz = sliceIndx(sliceIndx ~= 0);
        case 'co'  % coronal → vary Z
            sliceIndx = zeros(1, size(mask3d,3));
            for i = 1:size(mask3d,3)
                mask_slice = mask3d(:,:,i);
                sliceIndx(i) = i;
                if sum(mask_slice(:)) <= 25
                    mask3d(:,:,i) = 0;
                    sliceIndx(i) = 0;
                end
            end
            nz = sliceIndx(sliceIndx ~= 0);
    
        case 'sa'  % sagittal → vary X
            sliceIndx = zeros(1, size(mask3d,2));
            for i = 1:size(mask3d,2)
                mask_slice = squeeze(mask3d(:,i,:));
                sliceIndx(i) = i;
                if sum(mask_slice(:)) <= 25
                    mask3d(:,i,:) = 0;
                    sliceIndx(i) = 0;
                end
            end
            nz = sliceIndx(sliceIndx ~= 0);
    
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