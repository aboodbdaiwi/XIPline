
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
   
        if ispc
            % Windows
            XIPlineRoot = 'C:\XIPline';
        elseif ismac
            % macOS
            XIPlineRoot = fullfile(getenv('HOME'), 'XIPline');
        else
            % Linux
            XIPlineRoot = fullfile(getenv('HOME'), 'XIPline');
        end 

        destinationFolderPath = XIPlineRoot;

        % RH: removed hardcoded cd('C:\') — crashed on macOS ("Unable to change
        % current folder to '.../C:\'"). It only reset the cwd before creating
        % destinationFolderPath below, which mkdir/cd don't need.
        if ~exist(destinationFolderPath, 'dir')
            mkdir(destinationFolderPath);
        else        
            disp('HPXeAnalysisApp folder already exists in the destination folder.');
        end 
        cd(destinationFolderPath)
        % RH: was string-concatenated with a literal '\models', which produced a
        % bogus path (embedded backslash) on macOS/Linux; use fullfile instead.
        modelsFolderPath = fullfile(destinationFolderPath, 'models');
        if ~exist(modelsFolderPath, 'dir')
            mkdir(modelsFolderPath);
        else
            disp('models folder already exists in the destination folder.');
        end
    
        % Current function/script full path
        fullFilePath = mfilename('fullpath');
        
        % Function directory
        FunctionDirectory = fileparts(fullFilePath);
        [FunctionDirectory,~] = fileparts(FunctionDirectory);

        % RH: was hardcoded with backslash separators, breaking on macOS/Linux
        sourcemodel1Path = fullfile(FunctionDirectory, '+Segmentation', 'AutoSegmentation.py');

        try
            copyfile(sourcemodel1Path, destinationFolderPath); % force to copy
        catch
            disp('Python script is not there')
        end

        cd(modelsFolderPath);
        if ~exist(fullfile(destinationFolderPath, 'AutoSegmentation.py'), 'file')
            copyfile(sourcemodel1Path, destinationFolderPath);
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
        % RH: leading '\' made fullfile emit a bogus path on macOS/Linux
        % (fullfile('/a/b','\c.txt') -> '/a/b/\c.txt', not '/a/b/c.txt'),
        % so these files were never found there.
        satisfiedFile = fullfile(folderPath, 'python_requirement_satisfied.txt');
        pathFile = fullfile(folderPath, 'python_path.txt');

        % RH: read pythonPath unconditionally (not just when installing requirements)
        % so it's available below to actually invoke AutoSegmentation.py — MATLAB's
        % system() shell doesn't inherit the interactive PATH, so a bare "python"
        % command can silently fail to resolve (seen on macOS: "command not found").
        if exist(pathFile, 'file')
            fid = fopen(pathFile, 'r');
            if fid == -1
                error('Cannot open file: %s', pathFile);
            end
            pythonPath = strtrim(fgetl(fid));
            fclose(fid);
        else
            pythonPath = 'python3';
        end

        % Check if 'python_requirement_satisfied.txt' exists
        if exist(satisfiedFile, 'file')
            disp('Python requirements already satisfied.');
        else
            % If 'python_requirement_satisfied.txt' does not exist, check for 'python_path.txt'
            if exist(pathFile, 'file')
                % Display the Python path (this is where you can add additional steps if needed)
                fprintf('Python path: %s\n', pythonPath);

                terminate(pyenv)
                pyenv('Version', pythonPath);
                % RH: use "<pythonPath> -m pip ..." instead of bare "pip" — same PATH
                % issue as above; this also guarantees pip belongs to this interpreter.
                % RH: tensorflow==2.10.1/keras==2.10.0 (the Windows pins) have no
                % Apple Silicon wheel at all, so they can never install on macOS —
                % use the newer combo verified to load/run XIPline's existing legacy
                % .hdf5 models correctly under Keras 3 (tested against a real
                % production model). Same versions used for Linux since untested
                % there; tensorflow-macos is Apple-only so plain tensorflow is used
                % on Linux instead.
                pipStatus = zeros(1,5);
                if ispc
                    pipStatus(1) = system(sprintf('"%s" -m pip install numpy==1.26.4', pythonPath));
                    pipStatus(2) = system(sprintf('"%s" -m pip install keras==2.10.0', pythonPath)); % Specific version of Keras
                    pipStatus(3) = system(sprintf('"%s" -m pip install tensorflow==2.10.1', pythonPath)); % Specific version of TensorFlow
                    pipStatus(4) = system(sprintf('"%s" -m pip install nibabel', pythonPath));
                    pipStatus(5) = system(sprintf('"%s" -m pip install scipy', pythonPath));
                elseif ismac
                    pipStatus(1) = system(sprintf('"%s" -m pip install numpy==1.26.4', pythonPath));
                    pipStatus(2) = system(sprintf('"%s" -m pip install keras==3.15.1', pythonPath));
                    pipStatus(3) = system(sprintf('"%s" -m pip install tensorflow-macos==2.16.2', pythonPath));
                    pipStatus(4) = system(sprintf('"%s" -m pip install nibabel==5.3.0', pythonPath));
                    pipStatus(5) = system(sprintf('"%s" -m pip install scipy==1.13.1', pythonPath));
                else
                    pipStatus(1) = system(sprintf('"%s" -m pip install numpy==1.26.4', pythonPath));
                    pipStatus(2) = system(sprintf('"%s" -m pip install keras==3.15.1', pythonPath));
                    pipStatus(3) = system(sprintf('"%s" -m pip install tensorflow==2.16.2', pythonPath));
                    pipStatus(4) = system(sprintf('"%s" -m pip install nibabel==5.3.0', pythonPath));
                    pipStatus(5) = system(sprintf('"%s" -m pip install scipy==1.13.1', pythonPath));
                end
                %system('python -m pip uninstall -y protobuf')
                %system('python -m pip install "protobuf>=3.20.0,<3.21.0"')
                terminate(pyenv)

                % RH: only mark requirements satisfied if every install actually
                % succeeded — previously this was written unconditionally, so a
                % silently failed pip install (e.g. "pip: command not found") got
                % permanently marked as satisfied and never retried.
                if any(pipStatus ~= 0)
                    error(['One or more pip installs failed (exit codes: %s). ' ...
                        'Not marking Python requirements as satisfied — fix the ' ...
                        'issue and try again.'], mat2str(pipStatus));
                end

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
                        elseif strcmp(MainInput.Imagestosegment, 'Registered Proton')
                            SegmentType = 'vent_anat_2D_1ch_cor';
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
                        SegmentType = 'diff_2D_1ch'; % not supported yet
                    case 'transversal'
                        SegmentType = 'diff_2D_1ch'; % diff_2D_1ch_axi 2ch is not supported yet
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
       
        if strcmp(SegmentType, 'not_supported') == 0
            switch  MainInput.AIScript
                case 'Python'
                    % run python script
                    cd(destinationFolderPath)
                    % RH: was bare "python", which isn't on MATLAB's system() PATH on
                    % macOS; use the resolved interpreter from python_path.txt instead.
                    command = string(strcat(sprintf('"%s" AutoSegmentation.py', pythonPath),{' '},SegmentType));
                    LungMask = system(command);
                case 'Executable'
                    %Step 2: Run external recon executable
                    exePath = fullfile(XIPlineRoot, 'segmentation\AutoSegmentation.exe'); 
                    [status, cmdout] = system(['"', exePath, '"']);
                    if status ~= 0
                        error('Failed to run AutoSegmentation.exe:\n%s', cmdout);
                    end
            end            

            % load auto mask
            cd(destinationFolderPath);
            if exist('AutoMask.mat', 'file')
                % RH: hardcoded backslash separator, broke on macOS/Linux
                load(fullfile(destinationFolderPath, 'AutoMask.mat'));
                Mask = AutoMask > 0;

                % clean up the mask
                Mask = logical(Mask);
                
                % Remove slices containing <10% nonzero pixels
                minPixels = round(0.1 * size(Mask,1));
                nzCount = squeeze(sum(Mask,[1 2]));
                Mask(:,:,nzCount < minPixels) = 0;
                
                % Find and retain the longest continuous sequence of nonzero slices
                NZ = squeeze(any(Mask,[1 2]));
                
                d = diff([false; NZ; false]);
                runStart = find(d == 1);
                runEnd   = find(d == -1) - 1;
                
                if ~isempty(runStart)
                    [~,idx] = max(runEnd - runStart + 1);
                    keepRange = runStart(idx):runEnd(idx);
                
                    removeSlices = true(size(Mask,3),1);
                    removeSlices(keepRange) = false;
                    Mask(:,:,removeSlices) = 0;
                end              
                Mask = double(Mask);
                disp('auto mask process Completed.')
                %imslice(AutoMask)

                % resize mask to xenon image
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
                        % for i = 1:size(temp_mask,3)
                        %     temp_mask(:,:,i) = imerode(temp_mask(:,:,i), strel('square', 3));
                        % end
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

% if strcmp(SegmentType, 'not_supported') 
%     switch MainInput.AnalysisType
%         case 'Ventilation'
%            Ventilation.Mask = Mask;
%            Ventilation.LungMask = lungmask;
%            Ventilation.AirwayMask = airwaymask;                             
%         case 'Diffusion'
%            Diffusion.Mask = Mask;
%            Diffusion.LungMask = lungmask;
%            Diffusion.AirwayMask = airwaymask;
%         case 'GasExchange'
%     end
% 
% end

try
   Ventilation.UncorrectedImage = Ventilation.Image; 
catch

end
  
end

