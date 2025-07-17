
function [Proton,Ventilation,Diffusion, GasExchange] = PerformSegmentation(Proton,Ventilation,Diffusion,GasExchange,MainInput)
%   Inputs:
%          
%   Outputs:
%      LungMask
%
%   Example: 
%   Package: 
%
%   Author: Abdullah Bdaiwi 
%   Work email: abdullah.bdaiwi@cchmc.org
%   Personal email: abdaiwi89@gmail.com
%   Website: https://www.cincinnatichildrens.org/research/divisions/c/cpir
%% perform segmention
switch MainInput.SegmentationMethod
    case 'Manual' % ===========================================Manual=================================================  
            [Proton,Ventilation,Diffusion,GasExchange] = Segmentation.PerformManualThresholdSegmentation(Proton,Ventilation,Diffusion,GasExchange,MainInput);      
    case 'Threshold' % ========================================Threshold====================================================

        switch MainInput.SegmentAnatomy
            case 'Parenchyma' 
                %1) segment lungs
            case 'Airway'
                % 2) segment airways
                MainInput.Segment = 'airway'; % 'Airway'; || 'Parenchyma'
                MainInput.SegmentMethod = 'Manual'; % 'threshold' || 'manual'
                MainInput.SegmentManual = 'Freehand';
        end        
        
        [Proton,Ventilation,Diffusion,GasExchange] = Segmentation.PerformManualThresholdSegmentation(Proton,Ventilation,Diffusion,GasExchange,MainInput);

    case 'Auto' % ============================================Auto================================================
        automasking_folder = 'XIPline';
        destinationFolderPath = join(['C:\',automasking_folder]);
        cd('C:\');
        if ~exist(automasking_folder, 'dir')
            mkdir(destinationFolderPath);
        else        
            disp('HPXeAnalysisApp folder already exists in the destination folder.');
        end 
        cd(destinationFolderPath)
        models_folder = '\models';
        modelsFolderPath = join([destinationFolderPath,models_folder]);
        if ~exist(modelsFolderPath, 'dir')
            mkdir(modelsFolderPath);
        else
            disp('models folder already exists in the destination folder.');
        end
    
        % copy model to the HPXeAnalysisApp folder
        FunctionDirectory = which('XIPline');
        idcs = strfind(FunctionDirectory,filesep);%determine location of file separators
        FunctionDirectory = FunctionDirectory(1:idcs(end)-1);%remove file

        sourcemodel1Path = [FunctionDirectory,'\+Segmentation\AutoSegmentation.py'];
        sourcemodel2Path = [FunctionDirectory,'\+Segmentation\2DVent_Xe_axial_1000e_20250509.hdf5'];
        sourcemodel3Path = [FunctionDirectory,'\+Segmentation\2DVent_Xe_coronal_1000e_20250509.hdf5'];
        sourcemodel4Path = [FunctionDirectory,'\+Segmentation\2DVent_Xe_H_coronal_1000e_20230528.hdf5'];
        sourcemodel5Path = [FunctionDirectory,'\+Segmentation\3DGasExchange_Xe_100e_20250324.hdf5'];
        sourcemodel6Path = [FunctionDirectory,'\+Segmentation\3DGasExchange_Xe_HLR_100e_20250324.hdf5'];
        sourcemodel7Path = [FunctionDirectory,'\+Segmentation\2DDiff_Xe_axial_2000e_20240118.hdf5'];
        try
            copyfile(sourcemodel1Path, destinationFolderPath); % force to copy
        catch
            disp('Python script is not there')
        end

        cd(modelsFolderPath);
        if ~exist(fullfile(modelsFolderPath, '2DVent_Xe_axial_1000e_20250509.hdf5'), 'file') ||...
            ~exist(fullfile(modelsFolderPath, '2DVent_Xe_coronal_1000e_20250509.hdf5'), 'file') ||...
            ~exist(fullfile(modelsFolderPath, '2DVent_Xe_H_coronal_1000e_20230528.hdf5'), 'file') ||...
            ~exist(fullfile(modelsFolderPath, '3DGasExchange_Xe_100e_20250324.hdf5'), 'file') ||...
            ~exist(fullfile(modelsFolderPath, '3DGasExchange_Xe_HLR_100e_20250324.hdf5'), 'file') ||...
            ~exist(fullfile(modelsFolderPath, '2DDiff_Xe_axial_2000e_20240118.hdf5'), 'file') ||...~exist(fullfile(modelsFolderPath, 'AutoSegment_3DGasExchange_Xe_H_1000e.hdf5'), 'file') ||...
            ~exist(fullfile(destinationFolderPath, 'AutoSegmentation.py'), 'file')

            copyfile(sourcemodel1Path, destinationFolderPath);
            % copyfile(sourcemodel2Path, modelsFolderPath);
            % copyfile(sourcemodel3Path, modelsFolderPath);
            % copyfile(sourcemodel4Path, modelsFolderPath);
            % copyfile(sourcemodel5Path, modelsFolderPath);
            % copyfile(sourcemodel6Path, modelsFolderPath);
            % copyfile(sourcemodel7Path, modelsFolderPath);

            disp('File copied successfully.');
        else
            disp('models files already exist in the models folder.');
        end
        % delete old files
        fileName1 = 'AutoMask.mat';
        fileName2 = 'AutoMask.nii.gz';
        fullFilePath1 = fullfile(destinationFolderPath, fileName1);
        fullFilePath2 = fullfile(destinationFolderPath, fileName2);
        if exist(fullFilePath1, 'file') == 2
            delete(fullFilePath1);
            disp(['File ' fileName1 ' has been deleted successfully.']);
        else
            disp(['File ' fileName1 ' does not exist in the specified folder.']);
        end
        if exist(fullFilePath2, 'file') == 2
            delete(fullFilePath2);
            disp(['File ' fileName2 ' has been deleted successfully.']);
        else
            disp(['File ' fileName2 ' does not exist in the specified folder.']);
        end
        fileName = 'InputImage.mat';
        fullFilePath = fullfile(destinationFolderPath, fileName);
        if exist(fullFilePath, 'file') == 2
            delete(fullFilePath);
            disp(['File ' fileName ' has been deleted successfully.']);
        else
            disp(['File ' fileName ' does not exist in the specified folder.']);
        end

        %============================Making Lung Mask==========================
        disp('Making Lung Mask (trained ML model method)...');

%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %Need a raw CPython installation (NOT anaconda).  These commands below
%         % shouldn't need to be called if everthing is set up properly.
        
        % Define the folder path
        folderPath = destinationFolderPath;        
        % Define the file names
        satisfiedFile = fullfile(folderPath, '\python_requirement_satisfied.txt');
        pathFile = fullfile(folderPath, '\python_path.txt');
        
        % Check if 'python_requirement_satisfied.txt' exists
        if exist(satisfiedFile, 'file')
            disp('Python requirements already satisfied.');
        else
            % If 'python_requirement_satisfied.txt' does not exist, check for 'python_path.txt'
            if exist(pathFile, 'file')
                % Read the first line from 'python_path.txt'
                fid = fopen(pathFile, 'r');
                if fid == -1
                    error('Cannot open file: %s', pathFile);
                end
                pythonPath = fgetl(fid);
                fclose(fid);
                
                % Display the Python path (this is where you can add additional steps if needed)
                fprintf('Python path: %s\n', pythonPath);
                
                terminate(pyenv)
                pyenv('Version', pythonPath);
                system('pip install numpy')
                system('pip install keras==2.10.0'); % Specific version of Keras
                system('pip install tensorflow==2.10.1'); % Specific version of TensorFlow
                system('pip install nibabel')
                system('pip install scipy') 
                terminate(pyenv)
        
                % Create an empty 'python_requirement_satisfied.txt' file
                fid = fopen(satisfiedFile, 'w');
                if fid == -1
                    error('Cannot create file: %s', satisfiedFile);
                end
                fclose(fid);
                
                disp('Python requirements are now marked as satisfied.');
            else
                error('Neither python_requirement_satisfied.txt nor python_path.txt found.');
            end
        end

        cd(FunctionDirectory)
        pathToMod = fileparts(fullfile(destinationFolderPath, 'AutoSegmentation.py'));
        if count(py.sys.path,pathToMod)==0
            insert(py.sys.path,int32(0),pathToMod)
            disp('Uploaded AutoSegmentation')
        end
        % segmentation type
        % Check if MainInput.NoProtonImage is 'yes' or 'no', then convert to 1 or 0
        if ~isnumeric(MainInput.NoProtonImage)
            if strcmp(MainInput.NoProtonImage, 'yes')
                MainInput.NoProtonImage = 1;
            elseif strcmp(MainInput.NoProtonImage, 'no')
                MainInput.NoProtonImage = 0;
            end        
        end
        % segmentation type
        switch MainInput.AnalysisType
            case 'Ventilation'                            
               switch MainInput.SliceOrientation
                    case 'coronal' 
                        if MainInput.NoProtonImage == 0 && strcmp(MainInput.Imagestosegment, 'Proton & Xe Registered') == 1
                            SegmentType = 'vent_2D_2ch_cor'; %system('predict_mask_2DVent_w_H_coronal.exe')
                        elseif strcmp(MainInput.Imagestosegment, 'Xenon') == 1
                            SegmentType = 'vent_2D_1ch_cor'; 
                        else
                            SegmentType = 'vent_2D_1ch_cor'; %system('predict_mask_2DVent_wout_H_coronal.exe')           
                        end                       
                    case 'transversal'
                        SegmentType = 'vent_2D_1ch_axi'; % 2ch is not supported yet
                    case 'sagittal'
                        SegmentType = 'not_supported'; % not supported yet
                    case 'isotropic'
                        if (MainInput.NoProtonImage == 0) && (strcmp(MainInput.Imagestosegment, 'Proton & Xe Registered'))
                            SegmentType = 'gx_3D_2ch_iso'; %system('predict_mask_3DGasExchange_w_H_isotropic.exe')
                        elseif strcmp(MainInput.Imagestosegment, 'Xenon')
                            SegmentType = 'gx_3D_1ch_iso'; 
                        else 
                            SegmentType = 'gx_3D_1ch_iso'; %system('predict_mask_3DGasExchange_wout_H_isotropic.exe')           
                        end  
               end 
            case 'Diffusion'
               switch MainInput.SliceOrientation
                    case 'coronal' 
                        SegmentType = 'not_supported'; % not supported yet
                    case 'transversal'
                        SegmentType = 'diff_2D_1ch_axi'; % 2ch is not supported yet
                   case 'sagittal'
                        SegmentType = 'not_supported'; % not supported yet
                    case 'isotropic'
                        SegmentType = 'not_supported'; % not supported yet
               end 
            case 'GasExchange'
               switch MainInput.SliceOrientation
                    case 'coronal' 
                        SegmentType = 'not_supported'; % not supported yet
                    case 'transversal'
                        SegmentType = 'not_supported'; % not supported yet
                    case 'sagittal'
                        SegmentType = 'not_supported'; % not supported yet                        
                    case 'isotropic'
                        if (MainInput.NoProtonImage == 0) && (strcmp(MainInput.Imagestosegment, 'Proton & Xe Registered'))
                            SegmentType = 'gx_3D_2ch_iso'; %system('predict_mask_3DGasExchange_w_H_isotropic.exe')
                        elseif strcmp(MainInput.Imagestosegment, 'Xenon')
                            SegmentType = 'gx_3D_1ch_iso'; 
                        else 
                            SegmentType = 'gx_3D_1ch_iso'; %system('predict_mask_3DGasExchange_wout_H_isotropic.exe')           
                        end  
                end                                     
        end 
        MainInput.SegmentType = SegmentType;
        Segmentation.write_configSeg_file(MainInput);
        [~, MainInput] = Segmentation.preprocess_images_for_auto_segmentation(Proton,Ventilation,Diffusion,GasExchange,MainInput);
        cd(MainInput.AutoSegmentPath)  

        %Step 2: Run external recon executable
        exePath = 'C:\XIPline\segmentation\AutoSegmentation.exe';
        [status, cmdout] = system(['"', exePath, '"']);

        if status ~= 0
            error('Failed to run AutoSegmentation.exe:\n%s', cmdout);
        end

        if strcmp(SegmentType, 'not_supported') == 0
            % run python script 
            cd(destinationFolderPath)
            command = string(strcat('python AutoSegmentation.py',{' '},SegmentType));
            % LungMask = system(command);
            % load auto mask
            cd(destinationFolderPath);
            if exist('AutoMask.mat')
                load([destinationFolderPath,'\AutoMask.mat']);
                Mask = AutoMask > 0;
                disp('auto mask process Completed.')
                %imslice(AutoMask)
                switch MainInput.AnalysisType
                    case 'Ventilation'
                        if size(Mask,1) ~= size(Ventilation.Image,1) || size(Mask,2) ~= size(Ventilation.Image,2)
                            temp_mask = zeros(size(Ventilation.Image));
                            if strcmp(MainInput.SegmentType, 'gx_3D_1ch_iso') || strcmp(MainInput.SegmentType, 'gx_3D_2ch_iso')
                                temp_mask = imresize3(Mask, size(Ventilation.Image));
                            else
                                for i = 1:size(Ventilation.Image,3)
                                    temp_mask(:,:,i) = imresize(Mask(:,:,i),[size(Ventilation.Image,1),size(Ventilation.Image,2)]);
                                end
                            end
                            temp_mask = temp_mask > 0.5;
                        else
                            temp_mask = Mask;
                        end
                        Mask = double(temp_mask);
                        mask_existing = 1;
                    case 'Diffusion'    
                        if size(Mask,1) ~= size(Diffusion.Image,1)
                            temp_mask = zeros(size(Diffusion.Image,1),size(Diffusion.Image,2),size(Diffusion.Image,3));
                            for i = 1:size(Diffusion.Image,3)
                                temp_mask(:,:,i) = imresize(Mask(:,:,i),[size(Diffusion.Image,1),size(Diffusion.Image,2)]);
                            end
                            temp_mask = temp_mask > 0.5;
                        else
                            temp_mask = Mask;
                        end
                        for i = 1:size(temp_mask,3)
                            temp_mask(:,:,i) = imerode(temp_mask(:,:,i), strel('square', 3));
                        end
                        Mask = double(temp_mask);
                        mask_existing = 1;                        
                    case 'GasExchange'
                        if size(Mask,1) ~= size(GasExchange.VentImage,1)
                            temp_mask = imresize3(Mask, [size(GasExchange.VentImage,1),size(GasExchange.VentImage,2),size(GasExchange.VentImage,3)]);
                            temp_mask = temp_mask > 0.5;
                        else
                            temp_mask = Mask;
                        end
                        Mask = double(temp_mask);
                        mask_existing = 1;    
                end                             
                % store mask
                if mask_existing == 1
                    % assign lung and airway mask
                    Mask(isinf(Mask)|isnan(Mask)) = 0;
                    lungmask = zeros(size(Mask));
                    lungmask(Mask == 1)= 1;
                    airwaymask = zeros(size(Mask));
                    airwaymask(Mask == 2)= 1;  
                    
                    % store mask
                    switch MainInput.AnalysisType
                        case 'Ventilation'
                           Ventilation.Mask = Mask;
                           Ventilation.LungMask = lungmask;
                           Ventilation.AirwayMask = airwaymask;                             
                        case 'Diffusion'
                           Diffusion.Mask = Mask;
                           Diffusion.LungMask = lungmask;
                           Diffusion.AirwayMask = airwaymask;
                        case 'GasExchange'
                           GasExchange.Mask = Mask;
                           GasExchange.LungMask = lungmask;
                           Proton.LungMask = lungmask;
                           Proton.ProtonMaskRegistred = lungmask;
                           GasExchange.AirwayMask = airwaymask;  
%                             %View
%                             cd(GasExchange.outputpath);
%                             ProtonMaskMontage = figure('Name','Lung Mask');set(ProtonMaskMontage,'WindowState','minimized');
%                             montage(GasExchange.Mask,'DisplayRange',[0 1])%unregistered for these
%                             disp('Lung Mask Completed.')
%                             savefig('ProtonMaskMontage.fig')
%                             close(gcf)                           
                    end    
                end
            else
                disp('Auto mask does not exist in the destination folder.');
            end
        else
            disp('Auto segmentation is not supported for the selected settings')
        end
        cd(MainInput.XeDataLocation)
end 
try
   Ventilation.UncorrectedImage = Ventilation.Image; 
catch

end
  
end

