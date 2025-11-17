
function VV_mask = SegVentVolume(vent_img)

        Ventilation.Image = vent_img;
        Ventilation.SliceOrientation = 'coronal';
        MainInput.Institute = 'London';
        Ventilation.LungMask = ones(size(vent_img));
        [Ventilation] = VentilationFunctions.kmeans.clustering_3DASB(Ventilation, MainInput);
        clusters_3D = Ventilation.clusters_3D;
        VV_mask = double(clusters_3D > 0);

        [Nx, Ny, Nz] = size(VV_mask);

        % Adaptive threshold: 0.1% of slice area
        minPixels = round(0.001 * Nx * Ny);
        
        cleanMask3D = false(size(VV_mask));
        
        for z = 1:Nz
            slice = VV_mask(:, :, z);
            slice_clean = bwareaopen(slice, minPixels, 8);
            cleanMask3D(:, :, z) = slice_clean;
        end
        VV_mask = double(cleanMask3D);

        % figure; imslice(VV_mask)
end