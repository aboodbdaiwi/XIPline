function [MainInput, Proton, Ventilation, GasExchange] = AntsRegistration(MainInput, Proton, Ventilation, GasExchange)
    % Register proton images to xenon using ANTs
    % Optionally applies the transform to a second image (e.g., lung mask)
    
    clc;
    MainInput.transform_reg = false;  % set true to apply to second image (e.g., lung mask)
    TransformType = MainInput.TransformType;
    % Determine fixed image and optional second moving image
    switch MainInput.AnalysisType
        case 'Ventilation'
            fixed1 = Ventilation.Image;
            moving1 = double(Proton.Image);
            if MainInput.transform_reg
                moving2 = double(Ventilation.LungMask);
            end
        case 'GasExchange'
            fixed1 = GasExchange.VentImage;
            moving1 = double(Proton.ProtonImageHR);
            if MainInput.transform_reg
                moving2 = double(GasExchange.LungMask);
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
    
    % Determine required subfolder name
    if strcmp(MainInput.AnalysisType,'GasExchange')
        H_RecMatrix = Proton.H_RecMatrix;
        subFolder = 'GasExchange_Analysis';
    else
        subFolder = 'Ventilation_Analysis';
    end
    
    % Check if OutputPath already ends with the subfolder
    [~, lastFolder] = fileparts(MainInput.OutputPath);
    
    if ~strcmp(lastFolder, subFolder)
        DataLocation = fullfile(MainInput.OutputPath, subFolder);
    else
        DataLocation = MainInput.OutputPath;
    end
    
    % % Create directory if it does not exist
    % if ~exist(DataLocation, 'dir')
    %     mkdir(DataLocation);
    % end
    
    % Change to output directory
    cd(DataLocation);

    % create a lung mask
    LungMask = double(Segmentation.SegmentLungthresh(fixed1,1,1));
    se = strel('sphere', 3);  % Structuring element (sphere with radius 1)
    
    % Perform erosion on the 3D binary mask
    LungMask = imerode(LungMask, se);
    registration_path = 'C:\XIPline\Registration';
    if MainInput.SkipRegistration == 0
        % Prepare paths
        ANTSPath = fileparts(mfilename('fullpath'));
        tmp_path = fullfile(registration_path, 'tmp');
        if ~exist(tmp_path, 'dir'); mkdir(tmp_path); end
    
        pathStatic       = fullfile(tmp_path, 'image_static.nii');
        pathMoving1      = fullfile(tmp_path, 'image_moving1.nii');
        pathMoving2      = fullfile(tmp_path, 'image_moving2.nii');
        pathOutputPrefix = fullfile(tmp_path, 'thisTransform_');
        pathOutputMoving = fullfile(tmp_path, 'moving_reg.nii.gz');
        pathReg          = fullfile(ANTSPath, 'antsRegistration.exe');
        pathApply        = fullfile(ANTSPath, 'antsApplyTransforms.exe');
        tdata            = fullfile(tmp_path, 'thisTransform_0GenericAffine.mat');
    
        % Write NIfTI files
        niftiwrite(abs(fixed1), pathStatic);
        niftiwrite(abs(moving1), pathMoving1);
        if MainInput.transform_reg
            niftiwrite(abs(moving2), pathMoving2);
        end
    
        % Registration command
        fprintf('*** Running ANTs registration...\n');
        output_prefix = sprintf('[%s,%s]', pathOutputPrefix, pathOutputMoving);
        switch lower(TransformType)  % ensure case-insensitivity
            case 'translation'
                transformCmd = 'Translation[0.1]';
            case 'rigid'
                transformCmd = 'Rigid[0.1]';
            case 'similarity'
                transformCmd = 'Similarity[0.1]';
            case 'affine'
                transformCmd = 'Affine[0.1]';
            otherwise
                error('Unsupported TransformType: %s. Use translation, rigid, similarity, or affine.', TransformType);
        end
        
        % Assemble registration command using the selected transform
        cmd_register = sprintf(...
            ['"%s" --dimensionality 3 --float 0 --interpolation BSpline ' ...
            '--metric MI[%s,%s,1,32,Regular,1] --transform %s ' ...
            '--convergence [20x20x20,1e-6,20] --shrink-factors 4x2x1 ' ...
            '--smoothing-sigmas 0x0x0 --output %s --verbose 1'], ...
            pathReg, pathStatic, pathMoving1, transformCmd, output_prefix);

        % Execute registration
        [status, cmdout] = system(cmd_register);
        if status ~= 0
            error('ANTs registration failed:\n%s', cmdout);
        end
    
        MainInput.pathOutputMoving = pathOutputMoving;
    
        % Load registered image
        try
            A1 = LoadData.load_nii(MainInput.pathOutputMoving);
        catch
            A1 = LoadData.load_untouch_nii(MainInput.pathOutputMoving);
        end
    
        A = double(squeeze(A1.img));
        A = flip(A, 1);
        A = imrotate(A, 180);
        A = flip(A, 1);
        ProtonRegistered = A;    
        %imslice(A)

        Proton.tdata = load(tdata);
        % Optionally apply transform to second moving image
        if MainInput.transform_reg
            pathOutputMoving2 = fullfile(tmp_path, 'transform_reg.nii.gz');
            cmd_applyTransform = sprintf('"%s" -d 3 -e 0 -i "%s" -r "%s" -o "%s" -t "%s"', ...
                pathApply, pathMoving2, pathStatic, pathOutputMoving2, tdata);
            system(cmd_applyTransform);

             MainInput.pathOutputMoving2 = pathOutputMoving2;        
            % Load registered image
            try
                A1 = LoadData.load_nii(MainInput.pathOutputMoving2);
            catch
                A1 = LoadData.load_untouch_nii(MainInput.pathOutputMoving2);
            end        
            A = double(squeeze(A1.img));
            A = flip(A, 1);
            A = imrotate(A, 180);
            A = flip(A, 1);
            ProtonMaskRegistered = A;                 
        end
    else
        ProtonRegistered = moving1;
        if MainInput.transform_reg
            ProtonMaskRegistered = moving2;
            Proton.ProtonMaskRegistered = ProtonMaskRegistered;
        end
    end

    % Visualization
    fixedVolume = fixed1;
    for slice = 1:size(fixedVolume, 3)
        A = ProtonRegistered(:, :, slice);
        B = fixedVolume(:, :, slice);
        ProtonRegisteredColored(:, :, :, slice) = imfuse(A, B, 'falsecolor', 'ColorChannels', 'green-magenta');
    end
        
    % Determine required subfolder
    if strcmp(MainInput.AnalysisType, 'Ventilation')
        subFolder = 'Ventilation_Analysis';
    elseif strcmp(MainInput.AnalysisType, 'GasExchange')
        subFolder = 'GasExchange_Analysis';
    else
        error('Unsupported AnalysisType: %s', MainInput.AnalysisType);
    end
    
    % Check if OutputPath already ends with the subfolder
    [~, lastFolder] = fileparts(MainInput.OutputPath);
    
    if ~strcmp(lastFolder, subFolder)
        OutputPath = fullfile(MainInput.OutputPath, subFolder);
    else
        OutputPath = MainInput.OutputPath;
    end
    
    % % Ensure directory exists
    % if ~exist(OutputPath, 'dir')
    %     mkdir(OutputPath);
    % end
    
    % Create the output directory if it does not exist
    if ~exist(OutputPath, 'dir')
        mkdir(OutputPath);
    end

    Global.write_imshowpair(ProtonRegistered,fixedVolume,OutputPath)
    ProtonRegisteredColored = permute(ProtonRegisteredColored,[1 2 4 3]);
    Proton.ProtonRegistered = ProtonRegistered;
    Proton.ProtonRegisteredColored = ProtonRegisteredColored;
    if MainInput.transform_reg
        Proton.ProtonMaskRegistered = ProtonMaskRegistered;
    end

    if strcmp(MainInput.AnalysisType,'GasExchange')
        %Determine Slices to Plot
        NumPlotSlices = 7;
        if size(fixed1,3) < 7
            NumPlotSlices = size(fixed1,3);
        end
        %Coronal
        Co_MaskVoxels = squeeze(sum(LungMask,[1,2]));
        Co_StartIndex = find(Co_MaskVoxels,1,'first');
        Co_EndIndex = find(Co_MaskVoxels,1,'last');
        Co_MiddleIndex = round((Co_StartIndex+Co_EndIndex)/2);
        Co_Step = floor((Co_EndIndex-Co_StartIndex)/NumPlotSlices);
        Slices_Co = (Co_MiddleIndex-Co_Step*(NumPlotSlices-1)/2):Co_Step:(Co_MiddleIndex+Co_Step*(NumPlotSlices-1)/2);
        %Axial
        Ax_MaskVoxels = squeeze(sum(LungMask,[2,3]));
        Ax_StartIndex = find(Ax_MaskVoxels,1,'first');
        Ax_EndIndex = find(Ax_MaskVoxels,1,'last');
        Ax_MiddleIndex = round((Ax_StartIndex+Ax_EndIndex)/2);
        Ax_Step = floor((Ax_EndIndex-Ax_StartIndex)/NumPlotSlices);
        Slices_Ax = (Ax_MiddleIndex-Ax_Step*(NumPlotSlices-1)/2):Ax_Step:(Ax_MiddleIndex+Ax_Step*(NumPlotSlices-1)/2);
        fixedMasked = fixed1.*LungMask;
        %View
        Registrationfig = figure('Name','Registration','units','normalized','outerposition',...
            [0 0 1 4/NumPlotSlices]);set(Registrationfig,'WindowState','minimized');
        set(Registrationfig,'color','white','Units','inches','Position',[0.25 0.25 2*NumPlotSlices 4*2+1])
        tiledlayout(4,NumPlotSlices,'TileSpacing','none','Padding','compact');
        for slice=1:NumPlotSlices
            nexttile
            imshowpair(abs(moving1(:,:,Slices_Co(slice)))/max(abs(moving1(:))),...
                abs(fixed1(:,:,Slices_Co(slice)))/max(abs(fixedMasked(:))),...
                'Scaling','none','ColorChannels','red-cyan')
            if slice == round(NumPlotSlices/2)
                title('Before Registration','FontSize',24)
            end
        end
        for slice=1:NumPlotSlices
            nexttile
            imshowpair(fliplr(rot90(squeeze(abs(moving1(Slices_Ax(slice),:,:))),-1))/...
                max(abs(moving1(:))),fliplr(rot90(squeeze(abs(fixed1(Slices_Ax(slice),:,:))),-1))/...
                max(abs(fixedMasked(:))),'Scaling','none','ColorChannels','red-cyan')
        end
        for slice=1:NumPlotSlices
            nexttile
            imshowpair(abs(ProtonRegistered(:,:,Slices_Co(slice)))/...
                max(abs(ProtonRegistered(:))),abs(fixed1(:,:,Slices_Co(slice)))/...
                max(abs(fixedMasked(:))),'Scaling','none','ColorChannels','red-cyan')
            if slice == round(NumPlotSlices/2)
                title('After Registration','FontSize',24)
            end
        end
        for slice=1:NumPlotSlices
            nexttile
            imshowpair(fliplr(rot90(squeeze(abs(ProtonRegistered(Slices_Ax(slice),:,:))),-1))/...
                max(abs(ProtonRegistered(:))),fliplr(rot90(squeeze(abs(fixed1(Slices_Ax(slice),:,:))),-1))/...
                max(abs(fixedMasked(:))),'Scaling','none','ColorChannels','red-cyan')
        end
        savefig('Registerationfig.fig')
        close(gcf)
        
        
        Proton.LungMask = LungMask;
        try
            Proton.tform = tform;
            Proton.optimizer = optimizer;
        catch
            Proton.tform = [];
            Proton.optimizer = [];
        end
        Proton.Slices_Co = Slices_Co; 
        Proton.Slices_Ax = Slices_Ax; 
        
        Proton.ProtonRegistered = ProtonRegistered;
        if strcmp(MainInput.AnalysisType,'GasExchange')
            Proton.ProtonHRRegistered = ProtonRegistered;
        end
        Proton.ProtonRegisteredColored = permute(ProtonRegisteredColored,[1 2 4 3]);
        GasExchange.ProtonRegistered = ProtonRegistered;
        GasExchange.ProtonRegisteredColored = permute(ProtonRegisteredColored,[1 2 4 3]);
    
    end

end

