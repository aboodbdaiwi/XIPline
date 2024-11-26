

% function [moving1_reg, moving2_reg] = register_ants(Ventilation,Proton,MainInput)
    % Register images using ANTs executables
    %
    % Args:
    %   image_static: 3D array static image
    %   image_moving1: 3D array moving image 1
    %   image_moving2: 3D array moving image 2, transformed using 
    %                  the calculated transform between image_static and image_moving1.
    %
    % Returns:
    %   Registered images moving1_reg and moving2_reg
    clc;
    image_static = Ventilation.Image; 
    image_moving1 = Proton.Image;  
    % image_moving2 = Ventilation.LungMask;

    ANTSPath = mfilename('fullpath');
    idcs = strfind(ANTSPath,filesep);%determine location of file separators
    ANTSPath = [ANTSPath(1:idcs(end)-1),filesep];%remove file
    ANTSPath = 'D:\Github\XIPline\+Registration\ANTs';
    % Set paths

    tmp_path = fullfile(MainInput.XeDataLocation, 'tmp');
    if ~exist(tmp_path, 'dir')
        mkdir(tmp_path);
    end

    % Define file paths
    pathInputStatic = fullfile(tmp_path, 'image_static.nii');
    pathInputMoving1 = fullfile(tmp_path, 'image_moving1.nii');
    pathInputMoving2 = fullfile(tmp_path, 'image_moving2.nii');
    pathOutputPrefix = fullfile(tmp_path, 'thisTransform_');
    pathOutputMoving1 = fullfile(tmp_path, 'moving_reg.nii.gz');
    pathOutputMoving2 = fullfile(tmp_path, 'transform_reg.nii.gz');
    pathReg = fullfile(ANTSPath, 'antsRegistration.exe');
    pathApply = fullfile(ANTSPath, 'antsApplyTransforms.exe');
    tdata = fullfile(tmp_path, 'thisTransform_0GenericAffine.mat');

    % Save the inputs as NIfTI files
    niftiwrite(abs(image_static), pathInputStatic);
    niftiwrite(abs(image_moving1), pathInputMoving1);
    % niftiwrite(abs(image_moving2), pathInputMoving2);

    % Registration command
    fprintf('*** Using ANTs executables to register images ...\n');
    output_prefix = sprintf('[%s, %s]', pathOutputPrefix, pathOutputMoving1);
    cmd_register = sprintf(...
        ['%s --dimensionality 3 --float 0 --interpolation BSpline ' ...
        '--metric MI[%s,%s,1,32,Regular,1] --transform Affine[0.1] ' ... % Rigid or  Affine
        '--convergence [20x20x20,1e-6,20] --shrink-factors 4x2x1 ' ...
        '--smoothing-sigmas 0x0x0 --output %s --verbose 1'], ...
        pathReg, pathInputStatic, pathInputMoving1, output_prefix);

    % Execute registration
    system(cmd_register);

    % % Apply transformation to the second moving image
    % cmd_applyTransform = sprintf('%s -d 3 -e 0 -i %s -r %s -o %s -t %s', ...
    %     pathApply, pathInputMoving2, pathOutputMoving1, pathOutputMoving2, tdata);
    % system(cmd_applyTransform);

% figure; imslice(ProtonRegistered)
%% stor data
for slice =1:size(fixedVolume,3)
    A = ProtonRegistered(:,:,slice);
    B = fixedVolume(:,:,slice);
    ProtonRegisteredColored(:,:,:,slice) = imfuse(A,B,'falsecolor','ColorChannels','green-magenta'); %'red-cyan'|'green-magenta'
end
ProtonRegistered = ProtonRegistered(:,:,1:nSlice);
ProtonRegisteredColored = ProtonRegisteredColored(:,:,:,1:nSlice);
Proton.ProtonRegistered = ProtonRegistered;
Proton.Rmoving = Rmoving;
Proton.Rfixed = Rfixed;
Proton.optimizer = optimizer;
Proton.geomtform = geomtform;
Proton.ProtonRegisteredColored = permute(ProtonRegisteredColored,[1 2 4 3]);
% figure; orthosliceViewer(Proton.ProtonRegisteredColored)
%% view

if strcmp(MainInput.AnalysisType,'GasExchange')
    %Determine Slices to Plot
    NumPlotSlices = 7;
    if size(fixed,3) < 7
        NumPlotSlices = size(fixed,3);
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
    fixedMasked = fixed.*LungMask;
    %View
    Registrationfig = figure('Name','Registration','units','normalized','outerposition',[0 0 1 4/NumPlotSlices]);set(Registrationfig,'WindowState','minimized');
    set(Registrationfig,'color','white','Units','inches','Position',[0.25 0.25 2*NumPlotSlices 4*2+1])
    tiledlayout(4,NumPlotSlices,'TileSpacing','none','Padding','compact');
    for slice=1:NumPlotSlices
        nexttile
        imshowpair(abs(moving(:,:,Slices_Co(slice)))/max(abs(moving(:))),abs(fixed(:,:,Slices_Co(slice)))/max(abs(fixedMasked(:))),'Scaling','none','ColorChannels','red-cyan')
        if slice == round(NumPlotSlices/2)
            title('Before Registration','FontSize',24)
        end
    end
    for slice=1:NumPlotSlices
        nexttile
        imshowpair(fliplr(rot90(squeeze(abs(moving(Slices_Ax(slice),:,:))),-1))/max(abs(moving(:))),fliplr(rot90(squeeze(abs(fixed(Slices_Ax(slice),:,:))),-1))/max(abs(fixedMasked(:))),'Scaling','none','ColorChannels','red-cyan')
    end
    for slice=1:NumPlotSlices
        nexttile
        imshowpair(abs(ProtonRegistered(:,:,Slices_Co(slice)))/max(abs(ProtonRegistered(:))),abs(fixed(:,:,Slices_Co(slice)))/max(abs(fixedMasked(:))),'Scaling','none','ColorChannels','red-cyan')
        if slice == round(NumPlotSlices/2)
            title('After Registration','FontSize',24)
        end
    end
    for slice=1:NumPlotSlices
        nexttile
        imshowpair(fliplr(rot90(squeeze(abs(ProtonRegistered(Slices_Ax(slice),:,:))),-1))/max(abs(ProtonRegistered(:))),fliplr(rot90(squeeze(abs(fixed(Slices_Ax(slice),:,:))),-1))/max(abs(fixedMasked(:))),'Scaling','none','ColorChannels','red-cyan')
    end
    savefig('Registerationfig.fig')
    close(gcf)
    
    
    Proton.LungMask = LungMask;
    Proton.tform = tform;
    Proton.Slices_Co = Slices_Co; 
    Proton.Slices_Ax = Slices_Ax; 
    Proton.optimizer = optimizer;
    Proton.ProtonRegistered = ProtonRegistered;
    if strcmp(MainInput.AnalysisType,'GasExchange')
        Proton.ProtonHRRegistered = ProtonHRRegistered;
    end
    Proton.ProtonMaskRegistred = LungMask;
    Proton.ProtonRegisteredColored = permute(ProtonRegisteredColored,[1 2 4 3]);
    GasExchange.LungMask = LungMask;
    GasExchange.ProtonRegistered = ProtonRegistered;
    GasExchange.ProtonMaskRegistred = LungMask;
    GasExchange.ProtonRegisteredColored = permute(ProtonRegisteredColored,[1 2 4 3]);

end
% save images
%Tiffs (movingRegisteredVolume)

Global.write_imshowpair(ProtonRegistered,fixed(:,:,1:nSlice),DataLocation)

%     % Clean up temporary files
%     delete(pathInputStatic, pathInputMoving1, pathInputMoving2, ...
%         [pathOutputPrefix, '0GenericAffine.mat'], pathOutputMoving1, pathOutputMoving2);
% % end
