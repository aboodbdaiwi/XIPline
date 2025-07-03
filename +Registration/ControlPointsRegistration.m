
function [MainInput, Proton, Ventilation, GasExchange] = ControlPointsRegistration(MainInput, Proton, Ventilation, GasExchange)
% Register proton images to xenon images using manual control points (cpselect)
% Uses 2D slice-by-slice registration and optionally applies transform to second image (e.g., lung mask)
% Supported transform types: 'similarity', 'affine', 'projective', 'polynomial', etc.

    clc;
    
    if ~isfield(MainInput, 'TransformType')
        %error('MainInput must include a valid TransformType.');
        MainInput.TransformType = "affine";
    end
    TransformType = lower(MainInput.TransformType);
    
    % if ~ismember(TransformType, ["similarity", "reflectivesimilarity", "affine", "projective", "polynomial"])
    %     error('Unsupported TransformType: %s', TransformType);
    % end
    MainInput.transform_reg = false;
    switch MainInput.AnalysisType
        case 'Ventilation'
            fixed1 = Ventilation.Image;
            moving1 = double(Proton.Image);
            if MainInput.transform_reg
            end
        case 'GasExchange'
            fixed1 = GasExchange.VentImage;
            moving1 = double(Proton.ProtonImageHR);
            if MainInput.transform_reg
            end
        otherwise
            error('Unknown AnalysisType: %s', MainInput.AnalysisType);
    end

    if strcmp(MainInput.AnalysisType,'GasExchange') && size(moving1,1) > 130
        % Get the size of the arrays
        size1 = size(fixed1);
        size2 = size(moving1);
        
        % Calculate the starting and ending indices for cropping
        startIdx = (size2 - size1) / 2 + 1;
        endIdx = startIdx + size1 - 1;
        
        % Crop the center portion of array2
        moving1 = moving1(startIdx(1):endIdx(1), startIdx(2):endIdx(2), startIdx(3):endIdx(3));
        Proton.ProtonImageHR = Proton.ProtonImageHR(startIdx(1):endIdx(1), startIdx(2):endIdx(2), startIdx(3):endIdx(3));
        % Verify the size of the cropped array
        disp(size(moving1));
    end
    % Resize if dimensions mismatch
    if ~isequal(size(moving1), size(fixed1))
        moving1 = imresize3(moving1, size(fixed1));
    end

    % Manual slice-by-slice registration
    numSlices = size(fixed1, 3);
    ProtonRegistered = zeros(size(fixed1), 'like', moving1);

    for slice = 1:numSlices
        A = moving1(:, :, slice);
        B = fixed1(:, :, slice);
        fprintf('Select control points for slice %d/%d\n', slice, numSlices);
        [mp, fp] = cpselect(A, B, 'Wait', true);

        if size(mp, 1) < 3 || size(fp, 1) < 3
            warning('Insufficient control points for slice %d. Skipping registration.', slice);
            ProtonRegistered(:, :, slice) = A;
            continue;
        end

        try
            tform = fitgeotform2d(mp, fp, TransformType);
        catch ME
            warning('Transformation failed for slice %d: %s', slice, ME.message);
            ProtonRegistered(:, :, slice) = A;
            continue;
        end
        
        Rfixed = imref2d(size(B));
        ProtonRegistered(:, :, slice) = imwarp(A, tform, 'OutputView', Rfixed);
    end

    % Visualization
    for slice = 1:numSlices
        A = ProtonRegistered(:, :, slice);
        B = fixed1(:, :, slice);
        ProtonRegisteredColored(:, :, :, slice) = imfuse(A, B, 'falsecolor', 'ColorChannels', 'green-magenta');
    end
    Global.write_imshowpair(ProtonRegistered,fixed1,MainInput.OutputPath)

    % Store results
    Proton.ProtonRegistered = ProtonRegistered;
    Proton.ProtonRegisteredColored = permute(ProtonRegisteredColored, [1 2 4 3]);

    % Save montage if available
    if isfield(MainInput, 'OutputPath') && isfolder(MainInput.OutputPath)
        Global.write_imshowpair(ProtonRegistered, fixed1, MainInput.OutputPath);
    end


end