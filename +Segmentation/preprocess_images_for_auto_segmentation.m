
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
    end
    save([MainInput.XeDataLocation '\preprocessed_Images.mat'],'Images');
    scriptLocation = fileparts(mfilename('fullpath'));
    MainInput.AutoSegmentPath = scriptLocation;
    %disp(scriptLocation)
disp('preprocessing completed')
end