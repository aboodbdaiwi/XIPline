function VV_mask = SegVentVolume(vent_img)

    Ventilation.Image = vent_img;
    Ventilation.SliceOrientation = 'coronal';
    MainInput.Institute = 'London';
    Ventilation.LungMask = ones(size(vent_img));

    [Ventilation] = VentilationFunctions.kmeans.clustering_3DASB(Ventilation, MainInput);
    clusters_3D = Ventilation.clusters_3D;
    VV_mask = double(clusters_3D > 0);

    [Nx, Ny, Nz] = size(VV_mask);

    % -------------------------------------------------------------
    % (1) Adaptive boundary erosion based on image size
    %     ~0.5–1.5% of in-plane dimension
    % -------------------------------------------------------------
    inPlaneSize = sqrt(Nx * Ny);              % characteristic length
    nErode = max(1, round(0.01 * inPlaneSize));  % adaptive pixels to remove

    se = strel('disk', nErode, 0);

    erodedMask = false(size(VV_mask));
    for z = 1:Nz
        erodedMask(:, :, z) = imerode(VV_mask(:, :, z), se);
    end

    VV_mask = double(erodedMask);

    % -------------------------------------------------------------
    % (2) Remove small isolated regions (slice-wise)
    % -------------------------------------------------------------
    minPixels = round(0.001 * Nx * Ny);   % 0.1% of slice area

    cleanMask3D = false(size(VV_mask));
    for z = 1:Nz
        slice = VV_mask(:, :, z);
        cleanMask3D(:, :, z) = bwareaopen(slice, minPixels, 8);
    end

    VV_mask = double(cleanMask3D);

end
