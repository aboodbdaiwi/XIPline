function [Ventilation] = calculate_VHI(Ventilation)
    % Function to calculate Coefficient of Variation (CV) maps and the Ventilation
    % Heterogeneity Index (VHI) for a 3D image using a 3x3 sliding window.
    % The calculations consider only the pixels where the binary mask equals 1.
    % Additionally, the neighboring pixels in the window are filtered to exclude zeros.
    %
    % Input:
    %   Ventilation: Struct with fields:
    %       - Image: A 3D matrix of ventilation data.
    %       - LungMask: Binary mask for the lung region.
    %       - VesselMask: Binary mask for vessels, excluded from analysis.
    %
    % Output:
    %   Ventilation: Struct with additional fields:
    %       - CV_maps: A 3D array containing the CV maps for all slices.
    %       - sliceMeanCV: Mean CV for each slice.
    %       - sliceVHI: Ventilation Heterogeneity Index (IQR of CV) for each slice.
    %       - overallMeanCV: Mean CV across the entire 3D image.
    %       - overallVHI: Ventilation Heterogeneity Index for the entire 3D array.
    %
    % Written by Abdullah S. Bdaiwi 11/23/2024

    % Extract input data and define the effective binary mask
    image = Ventilation.Image;
    maskarray = double(Ventilation.LungMask);
    
    if size(image, 1) > 128
        imageresized = zeros(128, 128, size(image, 3));
        maskresized = zeros(128, 128, size(image, 3));
        cleanMask = zeros(128, 128, size(image, 3)); % To store cleaned masks
        
        % Define a small structuring element for morphological operations
        se = strel('disk', 1); % Disk with radius 1
        
        for i = 1:size(image, 3)
            % Resize image and mask
            imageresized(:, :, i) = imresize(image(:, :, i), [128, 128]);
            maskresized(:, :, i) = imresize(maskarray(:, :, i), [128, 128]);
            
            % Convert the resized mask to binary
            binaryMask = maskresized(:, :, i) > 0.5; 
            
            % Apply binary opening to remove isolated pixels
            cleanedMask = imopen(binaryMask, se);
            
            % Store the cleaned mask
            cleanMask(:, :, i) = cleanedMask;
        end
    end

    image = imageresized;
    maskarray = double(cleanMask > 0);
    [rows, cols, slices] = size(image);

    NormMR = image.*(maskarray>0);
    incomplete = 0.33;
    defectmask = VentilationFunctions.medFilter(double(NormMR<(mean(NormMR(NormMR>0))*incomplete)*(maskarray>0)));

    % Initialize outputs
    CV_maps = zeros(rows, cols, slices);  % To store CV maps for all slices
    sliceMeanCV = zeros(slices, 1);       % Mean CV for each slice
    sliceVHI = zeros(slices, 1);          % VHI for each slice
    allCV_values = [];                    % To collect all CV values across 3D

    % Loop through each slice
    for sliceIdx = 1:slices
        % Extract the current slice and mask
        slice = image(:, :, sliceIdx);
        maskSlice = maskarray(:, :, sliceIdx) - defectmask(:, :, sliceIdx);

        % Initialize the CV map for the current slice
        CV_map = NaN(rows, cols);
        CV_values = [];

        % Slide the 3x3 window across the slice
        for i = 2:rows-1
            for j = 2:cols-1
                % Extract 3x3 window
                window = slice(i-1:i+1, j-1:j+1);
                windowMask = maskSlice(i-1:i+1, j-1:j+1);

                % Filter out zero pixels using the mask
                validValues = window(windowMask == 1);
                validValues(validValues == 0) = [];

                % Calculate CV only if there are valid neighbors
                if ~isempty(validValues) && windowMask(2, 2) == 1
                    mean_val = mean(validValues);
                    std_val = std(validValues);

                    % Calculate CV
                    if mean_val ~= 0
                        CV = std_val / mean_val;
                    else
                        CV = 0;
                    end

                    % Store the CV value in the CV map and the list
                    CV_map(i, j) = CV;
                    CV_values = [CV_values, CV];
                end
            end
        end

        % Store the CV map for the current slice
        CV_maps(:, :, sliceIdx) = CV_map;

        % Calculate mean CV and VHI (IQR of CV) for the current slice
        if ~isempty(CV_values)
            sliceMeanCV(sliceIdx) = mean(CV_values);
            sliceVHI(sliceIdx) = iqr(CV_values).*100;  % IQR = Q3 - Q1
            allCV_values = [allCV_values, CV_values];  % Append to global list
        else
            sliceMeanCV(sliceIdx) = NaN;
            sliceVHI(sliceIdx) = NaN;
        end
    end

    % Calculate overall metrics
    if ~isempty(allCV_values)
        overallMeanCV = mean(allCV_values);
        overallVHI = iqr(allCV_values).*100;
    else
        overallMeanCV = NaN;
        overallVHI = NaN;
    end

    % Display results
    disp('Slice-wise Mean CV and VHI:');
    disp(table((1:slices)', sliceMeanCV, sliceVHI, ...
        'VariableNames', {'Slice', 'Mean_CV', 'VHI'}));
    disp(['Overall Mean CV: ', num2str(overallMeanCV)]);
    disp(['Overall VHI: ', num2str(overallVHI)]);

    % Store results in the Ventilation struct
    Ventilation.CV_maps = CV_maps;
    Ventilation.sliceMeanCV = sliceMeanCV;
    Ventilation.sliceVHI = sliceVHI;
    Ventilation.overallMeanCV = overallMeanCV;
    Ventilation.overallVHI = overallVHI;
    Ventilation.CV_maps(isnan(Ventilation.CV_maps)) = 0; % Replace NaN with 0
    Ventilation.CV_maps(isinf(Ventilation.CV_maps)) = 0; % Replace Inf with 0    
end
