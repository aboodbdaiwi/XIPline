function [Image, parentPath, filename] = LoadData_Gas_VentDiff_Philips_GRE_R590(MainInput)
%   Modified to perform zero-padding in k-space to change image size
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
%   12/24/24 - Modified to use k-space zero padding to change resolution

    if strcmp(MainInput.ReconImageMode,'proton')
        parentPath = MainInput.HDataLocation;
        cd(MainInput.HDataLocation)
        filename = MainInput.H_name;
    else
        parentPath = MainInput.XeDataLocation;
        cd(MainInput.XeDataLocation)
        try
            filename = MainInput.Xe_name;
        catch
            filename = MainInput.XeFileName;
        end
    end

    ch_range = [1 1]; % assume one channel for data import by default
    filterim = 'Fermi'; % 'Fermi' or 'no'

    [ImgK,NoiK,kx_oversample_factor] = LoadData.load_philips_extr1_2D(filename,ch_range);
    diffK = ImgK(:,:,:,:);    
    
    if length(size(diffK)) == 4
        data_size = size(diffK);
        ky = data_size(1);
        kx = data_size(2);
        slices = data_size(3);
        nbs = data_size(4);
        
        if strcmp(filterim,'Fermi')
            filter = LoadData.fermi_filter_2D_AB(ky,kx,10);
        else
            filter = 1;
        end
        
        kx_size = kx/kx_oversample_factor;
        
        % Determine target size based on ky
        if ky <= 75
            target_size = 80;
        else
            target_size = 128;
        end
        
        Image = zeros(target_size, target_size, slices, nbs);
        
        for sl = 1:slices
            for b = 1:nbs
                % Apply filter to k-space
                data_slice = diffK(:,:,sl,b) .* filter;
                
                % Remove oversampling in k-space
                kspace_cropped = data_slice(:, (kx/2)-(kx_size/2)+1:(kx/2)+(kx_size/2));
                
                % Go to image space first
                img_temp = fftshift(fft2(kspace_cropped), 2);
                
                % Back to k-space (unshifted for proper padding)
                kspace_unshifted = ifft2(ifftshift(img_temp, 2));
                
                % Zero-pad k-space to target size
                kspace_padded = zeros(target_size, target_size);
                
                % Calculate starting positions for centering
                ky_start = floor((target_size - ky)/2) + 1;
                kx_start = floor((target_size - kx_size)/2) + 1;
                
                % Place original k-space in center
                kspace_padded(ky_start:ky_start+ky-1, kx_start:kx_start+kx_size-1) = kspace_unshifted;
                
                % Reconstruct from padded k-space
                recon_padded = fftshift(fft2(kspace_padded), 2);
                
                % Take magnitude and rotate
                diffimg = abs(recon_padded);
                diffimg = rot90(rot90(diffimg));
                
                Image(:,:,sl,b) = diffimg;
            end
        end
        
    elseif length(size(diffK)) <= 3
        data_size = size(diffK);
        ky = data_size(1);
        kx = data_size(2);
        if length(data_size) < 3
            slices = 1;
        else
            slices = data_size(3); 
        end
        
        if strcmp(filterim,'Fermi')
            filter = LoadData.fermi_filter_2D_AB(ky,kx,10);
        else
            filter = 1;
        end
        
        kx_size = kx/kx_oversample_factor;
        
        % Set target size to 128x128
        target_size = 128;
        
        Image = zeros(target_size, target_size, slices);
        
        for sl = 1:slices
            % Apply filter to k-space
            data_slice = diffK(:,:,sl) .* filter;
            
            % Remove oversampling in k-space
            kspace_cropped = data_slice(:, (kx/2)-(kx_size/2)+1:(kx/2)+(kx_size/2));
            
            % Go to image space first
            img_temp = fftshift(fft2(kspace_cropped), 2);
            
            % Back to k-space (unshifted for proper padding)
            kspace_unshifted = ifft2(ifftshift(img_temp, 2));
            
            % Zero-pad k-space to target size
            kspace_padded = zeros(target_size, target_size);
            
            % Calculate starting positions for centering
            ky_start = floor((target_size - ky)/2) + 1;
            kx_start = floor((target_size - kx_size)/2) + 1;
            
            % Place original k-space in center
            kspace_padded(ky_start:ky_start+ky-1, kx_start:kx_start+kx_size-1) = kspace_unshifted;
            
            % Reconstruct from padded k-space
            recon_padded = fftshift(fft2(kspace_padded), 2);
            
            % Take magnitude and rotate
            diffimg = abs(recon_padded);
            diffimg = rot90(rot90(diffimg));
            
            Image(:,:,sl) = diffimg;
        end
    end
    
end