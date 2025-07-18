
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
            if strcmp(MainInput.SegmentType, 'gx_3D_1ch_iso') || strcmp(MainInput.SegmentType, 'gx_3D_2ch_iso')
                Im_size = [112, 112, 112];
                if ~isequal(size(Ventilation.Image), Im_size)
                    Xe_Img = imresize3(Ventilation.Image, Im_size);
                else
                    Xe_Img = Ventilation.Image;
                end
            else
                Im_size = 256;
                if size(Ventilation.Image,1) ~= Im_size || size(Ventilation.Image,2) ~= Im_size
                    Xe_Img = zeros(Im_size,Im_size,size(Ventilation.Image,3));
                    for i = 1:size(Ventilation.Image,3)
                         resImg = imresize(Ventilation.Image(:,:,i),[Im_size,Im_size]);
                         %resImg = (resImg - min(resImg(:)))./(max(resImg(:)) - min(resImg(:)));
                         Xe_Img(:,:,i) = resImg;
                    end
                else
                    Xe_Img = Ventilation.Image;               
                end
            end

            % run N4
            [Xe_Img, ~] = Segmentation.N4_bias_correction(Xe_Img, MainInput.XeDataLocation);

            % Normalize each slice independently
            Xe_Img = Xe_Img/max(Xe_Img(:));
            % for sl = 1:size(Xe_Img, 3)
            %     max_val = max(Xe_Img(:, :, sl), [], 'all');
            %     if max_val > 0
            %         Xe_Img(:, :, sl) = Xe_Img(:, :, sl) ./ max_val;
            %     end
            % end

            % resize proton images
            if (strcmp(MainInput.NoProtonImage, 'no') || MainInput.NoProtonImage == 0) && strcmp(MainInput.Imagestosegment, 'Proton & Xe Registered') 
                % resize
                if strcmp(MainInput.SegmentType, 'gx_3D_1ch_iso') || strcmp(MainInput.SegmentType, 'gx_3D_2ch_iso')
                    Im_size = [112, 112, 112];
                    if ~isequal(size(Proton.ProtonRegistered), Im_size)
                        H_Img = imresize3(Proton.ProtonRegistered, Im_size);
                    else
                        H_Img = Proton.ProtonRegistered;
                    end
                else
                    if size(Proton.ProtonRegistered,1) ~= Im_size || size(Proton.ProtonRegistered,2) ~= Im_size
                        H_Img = zeros(Im_size,Im_size,size(Proton.ProtonRegistered,3));
                        for i = 1:size(Proton.ProtonRegistered,3)
                            H_Img(:,:,i) = imresize(Proton.ProtonRegistered(:,:,i),[Im_size,Im_size]);
                        end
                    else
                        H_Img = Proton.ProtonRegistered;
                    end
                end
                % normalize
                H_Img = H_Img./max(H_Img(:));
            end
            % make stacks 
            if (strcmp(MainInput.NoProtonImage, 'no') || MainInput.NoProtonImage == 0) &&...
                    strcmp(MainInput.Imagestosegment, 'Proton & Xe Registered') &&...
                    strcmp(MainInput.SegmentType, 'gx_3D_1ch_iso')
                    Images = zeros(1,Im_size(1),Im_size(2),Im_size(3),2);
                    Images(1,:,:,:,1) = Xe_Img;
                    Images(1,:,:,:,2) = H_Img;
            elseif strcmp(MainInput.SegmentType, 'gx_3D_1ch_iso')
                Images = zeros(1,Im_size(1),Im_size(2),Im_size(3),1);
                Images(1,:,:,:,1) = Xe_Img;
            else
                Images = zeros(size(Ventilation.Image,3),Im_size,Im_size,3);
                for j = 1:size(Ventilation.Image,3)
                    Images(j,:,:,1) = Xe_Img(:,:,j);
                    Images(j,:,:,2) = Xe_Img(:,:,j);
                    Images(j,:,:,3) = Xe_Img(:,:,j);
                    if (strcmp(MainInput.NoProtonImage, 'no') || MainInput.NoProtonImage == 0) && strcmp(MainInput.SliceOrientation, 'coronal') && strcmp(MainInput.Imagestosegment, 'Proton & Xe Registered') 
                        Images(j,:,:,2) = H_Img(:,:,j);
                    end
                end 
            end

        case 'Diffusion'
            % resize xenon images
            Im_size = 64;
            if size(Diffusion.Image,1) ~= Im_size || size(Diffusion.Image,2) ~= Im_size
                Xe_Img = zeros(Im_size,Im_size,size(Diffusion.Image,3));
                for i = 1:size(Diffusion.Image,3)
                    Xe_Img(:,:,i) = imresize(Diffusion.Image(:,:,i,1),[Im_size,Im_size]);
                end
            else
                Xe_Img = Diffusion.Image(:,:,:,1);
            end
            % run N4
            % [Xe_Img, ~] = Segmentation.N4_bias_correction(Xe_Img, MainInput.XeDataLocation);
            
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
            % if strcmp(MainInput.Institute, 'CCHMC')
            %     Im_size = [112, 112, 112];
            % else
            %     Im_size = [112, 112, 112];
            % end
            Im_size = [112, 112, 112];
            if ~isequal(size(GasExchange.VentImage), Im_size)
                Xe_Img = imresize3(GasExchange.VentImage, Im_size);
            else
                Xe_Img = GasExchange.VentImage;
            end
            % run N4
            %[Xe_Img, ~] = Segmentation.N4_bias_correction(Xe_Img, MainInput.XeDataLocation);            
            Xe_Img = double(Xe_Img./max(Xe_Img(:)));
            
            % resize proton images
            if strcmp(MainInput.NoProtonImage, 'no') || MainInput.NoProtonImage == 0
                % resize
                try
                    ProtonImage = Proton.ProtonLRRegistered;
                catch
                    ProtonImage = Proton.ProtonRegistered;
                end
                if ~isequal(size(ProtonImage), Im_size)
                    H_Img = imresize3(ProtonImage, Im_size);
                else
                    H_Img = ProtonImage;
                end
                % normalize
                H_Img = double(H_Img./max(H_Img(:)));
            end
            % make stacks 
            if (strcmp(MainInput.NoProtonImage, 'no') || MainInput.NoProtonImage == 0) && strcmp(MainInput.SegmentType, 'gx_3D_2ch_iso')
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
    automasking_folder = 'XIPline';
    destinationFolderPath = join(['C:\',automasking_folder]);
    cd(destinationFolderPath)
    save([destinationFolderPath '\InputImage.mat'],'Images');
%     scriptLocation = fileparts(mfilename('fullpath'));
    MainInput.AutoSegmentPath = destinationFolderPath;
    %disp(scriptLocation)
disp('preprocessing completed')
end