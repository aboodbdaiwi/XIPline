
function [Image, parentPath, filename] = LoadData_Gas_VentDiff_Philips_GRE(MainInput)
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

    % DataFiles = dir([FullFilePath,'\*.data']);
    % FileNames = DataFiles.name;

    if strcmp(MainInput.ReconImageMode,'proton')
        parentPath = MainInput.HDataLocation;
        cd(MainInput.HDataLocation)
        filename = MainInput.H_name;
    else
        parentPath = MainInput.XeDataLocation;
        cd(MainInput.XeDataLocation)
        filename = MainInput.Xe_name;
    end
  
    ch_range = [1 1]; % assume one channel for data import by default
    filterim = 'Fermi';

    try
        [ImgK,NoiK,kx_oversample_factor] = LoadData.load_philips_extr1_2D(filename,ch_range);
    catch
        [ImgK,NoiK,kx_oversample_factor] = LoadData.load_philips_extr1_2D_R32(filename,ch_range);
    end
    diffK = squeeze(ImgK); 
    % figure; imslice(abs(diffK)); 

    if  length(size(diffK)) == 4
        data_size = size(diffK);
        ky = data_size(1);
        kx =data_size(2);
        slices = data_size(3);
        nbs = data_size(4);
        diffimgcmplx = zeros(ky,kx/kx_oversample_factor,slices, nbs);
        if strcmp(filterim,'Fermi') == 1
                filter= LoadData.fermi_filter_2D_AB(ky,kx,10);
        else
                filter = 1;
        end
        if ky <= 75
            imsize = 80;
            Image = zeros(imsize,imsize,slices,nbs);
        else
            Image = zeros(ky,ky,slices,nbs);
        end
        for sl =1:slices
            for b = 1:nbs
                data_slice = diffK(:,:,sl,b).*filter;
                recon_slice = fftshift(fft2(data_slice),2);
                recon_slice(:,1:kx/(2*kx_oversample_factor))=[]; % crop extra pixels from oversample factor
                recon_slice(:, (1+kx/kx_oversample_factor):end)=[]; 
                diffimgcmplx(:,:,sl,b) = recon_slice;
                diffimg(:,:,sl,b) = abs(diffimgcmplx(:,:,sl,b)); % make magnitude data
               
                % Get original image size
                [origY, origX] = size(diffimg(:,:,sl,b));
                
                % Compute uniform scaling factor
                if max(origY, origX) <= 75
                    imsize = 80;
                
                    scale = imsize / max(origY, origX);
                    newY  = round(origY * scale);
                    newX  = round(origX * scale);
                
                    % Resize with uniform scale
                    resized_img = imresize(diffimg(:,:,sl,b), [newY, newX], 'bicubic');
                
                    % Pad to imsize x imsize
                    padded_img = zeros(imsize, imsize);
                
                    startY = floor((imsize - newY)/2) + 1;
                    startX = floor((imsize - newX)/2) + 1;
                
                    padded_img(startY:startY+newY-1, startX:startX+newX-1) = resized_img;
                else
                    padded_img = diffimg(:,:,sl,b);
                end

                Image(:,:,sl,b) = padded_img;

            end
        end
        Image = rot90(rot90(Image));    
        % figure; imslice(Image);       
    elseif length(size(diffK)) <= 3
            data_size = size(diffK);
            ky = data_size(1);
            kx =data_size(2);
            if length(data_size) < 3
                slices = 1;
            else
                slices = data_size(3); 
            end
            diffimgcmplx = zeros(ky,round(kx/kx_oversample_factor),slices);
            if strcmp(filterim,'Fermi') == 1
                    filter= LoadData.fermi_filter_2D_AB(ky,kx,10);
            else
                    filter = 1;
            end
            Image = zeros(128,128,slices);
            for sl =1:slices
                    data_slice = diffK(:,:,sl).*filter;
                    recon_slice = fftshift(fft2(data_slice),2);
                    recon_slice(:,1:kx/(2*kx_oversample_factor))=[]; % crop extra pixels from oversample factor
                    recon_slice(:, (1+kx/kx_oversample_factor):end)=[]; 
                    diffimgcmplx(:,:,sl) = recon_slice;
                    diffimg(:,:,sl) = abs(diffimgcmplx(:,:,sl)); % make magnitude data
                    
                    [origY, origX] = size(diffimg(:,:,sl));
                    
                    % Compute uniform scaling factor
                    scale = 128 / max(origY, origX);
                    newY  = round(origY * scale);
                    newX  = round(origX * scale);
                    
                    % Resize using uniform scale
                    resized_img = imresize(diffimg(:,:,sl), [newY, newX]);
                    
                    % ---- FORCE BORDER PIXELS TO ZERO ----
                    resized_img(1, :)   = 0;
                    resized_img(end, :) = 0;
                    resized_img(:, 1)   = 0;
                    resized_img(:, end) = 0;
                    
                    % Pad to 128×128 with zeros
                    padded_img = zeros(128, 128);
                    startY = floor((128 - newY)/2) + 1;
                    startX = floor((128 - newX)/2) + 1;
                    
                    padded_img(startY:startY+newY-1, startX:startX+newX-1) = resized_img;
                    
                    Image(:,:,sl) = padded_img;

            end
            Image = rot90(rot90(Image));  
    end
end