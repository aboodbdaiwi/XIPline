
function [Images, MainInput] = preprocess_images_for_auto_segmentation(Proton,Ventilation,Diffusion,GasExchange,MainInput)
%% MATLAB script to preprocess_images_for_auto_segmentation
%
% Authors: Abdullah S. Bdaiwi
%         
% Date: 06/21/2023.

    switch MainInput.AnalysisType
        case'Ventilation'            
            % resize xenon images
            if size(Ventilation.Image,1) ~= 256
                Xe_Img = zeros(256,256,size(Ventilation.Image,3));
                for i = 1:size(Ventilation.Image,3)
                    Xe_Img(:,:,i) = imresize(Ventilation.Image(:,:,i),[256,256]);
                end
            end
            % normalize
            Xe_Img = Xe_Img./max(Xe_Img(:));
            % resize proton images
            if strcmp(MainInput.NoProtonImage, 'no') == 1
                % resize
                if size(Proton.Image,1) ~= 256
                    H_Img = zeros(256,256,size(Proton.Image,3));
                    for i = 1:size(Proton.Image,3)
                        H_Img(:,:,i) = imresize(Proton.ProtonRegistered(:,:,i),[256,256]);
                    end
                end
                % normalize
                H_Img = H_Img./max(H_Img(:));
            end
            % make stacks 
            Images = zeros(size(Ventilation.Image,3),256,256,3);
            for j = 1:size(Ventilation.Image,3)
                Images(j,:,:,1) = Xe_Img(:,:,j);
                Images(j,:,:,2) = Xe_Img(:,:,j);
                Images(j,:,:,3) = Xe_Img(:,:,j);
                if strcmp(MainInput.NoProtonImage, 'no') == 1
                    Images(j,:,:,2) = H_Img(:,:,j);
                end
            end
            
        
        case 'Diffusion'
        case 'GasExchange'
            % resize xenon images
            if size(GasExchange.VentImage,1) ~= 112
                Xe_Img = zeros(112,112,112);
                for i = 1:size(GasExchange.VentImage,3)
                    Xe_Img(:,:,i) = imresize(GasExchange.VentImage(:,:,i),[112,112]);
                end
            else
                Xe_Img = GasExchange.VentImage;
            end
            Xe_Img = Xe_Img./max(Xe_Img(:));
            % resize proton images
            if MainInput.NoProtonImage == 0
                % resize
                if size(Proton.ProtonRegistered,1) ~= 112
                    H_Img = zeros(112,112,112);
                    for i = 1:size(Proton.ProtonRegistered,3)
                        H_Img(:,:,i) = imresize(Proton.ProtonRegistered(:,:,i),[112,112]);
                    end
                else
                    H_Img = Proton.ProtonRegistered;
                end
                % normalize
                H_Img = H_Img./max(H_Img(:));
            end
            % make stacks 
            if MainInput.NoProtonImage == 0
                Images = zeros(1,112,112,112,2);
            else
                Images = zeros(1,112,112,112,1);
            end
            
            for j = 1:size(Ventilation.Image,3)
                Images(1,:,:,:,1) = Xe_Img;
                if MainInput.NoProtonImage == 0
                    Images(1,:,:,:,2) = H_Img;
                end
            end
            % imslice(squeeze(Images(1,:,:,:,:)))
    end
    save([MainInput.XeDataLocation '\Mat2Py_preprocessing.mat'],'Images');
    scriptLocation = fileparts(mfilename('fullpath'));
    MainInput.AutoSegmentPath = scriptLocation;
    %disp(scriptLocation)
disp('preprocessing completed')
end