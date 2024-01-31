
function [Images, MainInput] = preprocess_images_for_auto_segmentation(Proton,Ventilation,Diffusion,GasExchange,MainInput)
%% MATLAB script to preprocess_images_for_auto_segmentation
%
%   Inputs:
%      
%   Outputs:
%                   
%   Package: https://github.com/aboodbdaiwi/HP129Xe_Analysis_App
%
%   Author: Abdullah S. Bdaiwi
%   Work email: abdullah.bdaiwi@cchmc.org
%   Personal email: abdaiwi89@gmail.com
%   Website: https://www.cincinnatichildrens.org/research/divisions/c/cpir
%
%   Please add updates at the end. Ex: 3/10/24 - ASB: update ....   
% Date: 06/21/2023.

    switch MainInput.AnalysisType
        case'Ventilation'            
            % resize xenon images
            Im_size = 256;
            if size(Ventilation.Image,1) ~= Im_size
                Xe_Img = zeros(Im_size,Im_size,size(Ventilation.Image,3));
                for i = 1:size(Ventilation.Image,3)
                    Xe_Img(:,:,i) = imresize(Ventilation.Image(:,:,i),[Im_size,Im_size]);
                end
            end
            % normalize
            Xe_Img = Xe_Img./max(Xe_Img(:));
            % resize proton images
            if strcmp(MainInput.NoProtonImage, 'no') == 1
                % resize
                if size(Proton.Image,1) ~= Im_size
                    H_Img = zeros(Im_size,Im_size,size(Proton.Image,3));
                    for i = 1:size(Proton.Image,3)
                        H_Img(:,:,i) = imresize(Proton.ProtonRegistered(:,:,i),[Im_size,Im_size]);
                    end
                end
                % normalize
                H_Img = H_Img./max(H_Img(:));
            end
            % make stacks 
            Images = zeros(size(Ventilation.Image,3),Im_size,Im_size,3);
            for j = 1:size(Ventilation.Image,3)
                Images(j,:,:,1) = Xe_Img(:,:,j);
                Images(j,:,:,2) = Xe_Img(:,:,j);
                Images(j,:,:,3) = Xe_Img(:,:,j);
                if strcmp(MainInput.NoProtonImage, 'no') == 1
                    Images(j,:,:,2) = H_Img(:,:,j);
                end
            end     
        case 'Diffusion'
            % resize xenon images
            Im_size = 64;
            if size(Diffusion.Image,1) ~= Im_size
                Xe_Img = zeros(Im_size,Im_size,size(Diffusion.Image,3));
                for i = 1:size(Diffusion.Image,3)
                    Xe_Img(:,:,i) = imresize(Diffusion.Image(:,:,i,1),[Im_size,Im_size]);
                end
            end
            % normalize
            Xe_Img = Xe_Img./max(Xe_Img(:));

            % make stacks 
            Images = zeros(size(Diffusion.Image,3),Im_size,Im_size,3);
            for j = 1:size(Diffusion.Image,3)
                Images(j,:,:,1) = Xe_Img(:,:,j);
                Images(j,:,:,2) = Xe_Img(:,:,j);
                Images(j,:,:,3) = Xe_Img(:,:,j);
            end        
        case 'GasExchange'
            % resize xenon images
            if strcmp(MainInput.Institute, 'CCHMC')
                Im_size = [11, 112, 112];
            else
                Im_size = [112, 112, 112];
            end
            if size(GasExchange.VentImage,1) ~= Im_size
                Xe_Img = imresize3(GasExchange.VentImage, Im_size);
            else
                Xe_Img = GasExchange.VentImage;
            end
            Xe_Img = double(Xe_Img./max(Xe_Img(:)));
            % resize proton images
            if strcmp(MainInput.NoProtonImage, 'no')
                % resize
                if size(Proton.ProtonRegistered,1) ~= Im_size
                    H_Img = imresize3(Proton.ProtonRegistered, Im_size);
                else
                    H_Img = Proton.ProtonRegistered;
                end
                % normalize
                H_Img = double(H_Img./max(H_Img(:)));
            end
            % make stacks 
            if strcmp(MainInput.NoProtonImage, 'no')
                Images = zeros(1,Im_size(1),Im_size(2),Im_size(3),2);
                Images(1,:,:,:,1) = Xe_Img;
                Images(1,:,:,:,2) = H_Img;
            else
                Images = zeros(1,Im_size(1),Im_size(2),Im_size(3),1);
                Images(1,:,:,:,1) = Xe_Img;
            end
            % imslice(squeeze(Images(1,:,:,:,:)))
    end
%     save([MainInput.XeDataLocation '\Mat2Py_preprocessing.mat'],'Images');
    automasking_folder = 'HPXeAnalysisApp';
    destinationFolderPath = join(['C:\',automasking_folder]);
    cd(destinationFolderPath)
    save([destinationFolderPath '\InputImage.mat'],'Images');
%     scriptLocation = fileparts(mfilename('fullpath'));
    MainInput.AutoSegmentPath = destinationFolderPath;
    %disp(scriptLocation)
disp('preprocessing completed')
end