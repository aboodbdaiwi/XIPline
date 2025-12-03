function [VDP, VDPmap] = VDPAKms(image, mask, SliceOrientation, ImageResolution)
    Ventilation.Image = image;
    maskarray = double(mask);
    Ventilation.SliceOrientation = SliceOrientation;
    MainInput.SliceOrientation = SliceOrientation;
    MainInput.Institute = 'XeCTC';
    Ventilation.LungMask = maskarray;
    temp_handles.proton_used = 'SSFSE'; % hard code for now
    temp_handles.proton = zeros(size(image));
    if size(maskarray, 1) <= 128
        se = strel('disk', 2);
    elseif size(maskarray, 1) <= 192
        se = strel('disk', 4);
    else
        se = strel('disk', 6);    
    end
    dilated_mask = zeros(size(maskarray));
    for i =1:size(maskarray,3)
        dilated_mask(:,:,i) = double(imdilate(maskarray(:,:,i), se) > 0);
    end
    temp_handles.he = Ventilation.Image; %xe.img;
    temp_handles.lungsmask = dilated_mask; % maskarray
    temp_handles.vessel_mask = zeros(size(image));
    temp_handles.lobemask = maskarray;
    temp_handles.ventmaskbinary = maskarray;
    temp_handles.hist_mu = 385;
    temp_handles.frangi_thresh = 0.25;%0.15;
    temp_handles.frangi_sigmas = 1;%1.5;%
    temp_handles.pid = 'XeCTC';
    temp_handles.vessel_removal = false;
    temp_handles.vessel_thresh = 80;
    temp_handles.view = MainInput.SliceOrientation;
    temp_handles.helium_voxelsize = prod(ImageResolution)*1e-6;
    [VDPmap,~,~,~,~,~,ventilation_mask] = VentilationFunctions.Akmeans.ventilation_defect_by_kmeans(temp_handles);
    VDPmap = VDPmap.*maskarray;
    ventilation_mask = ventilation_mask.*maskarray;
    ventP_raw_cor = VentilationFunctions.Akmeans.classify_ventilation_levles(ventilation_mask,temp_handles.lungsmask);
    VDP = ventP_raw_cor(1);
end
