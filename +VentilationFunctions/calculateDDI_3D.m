
function [Ventilation] = calculateDDI_3D(Ventilation,Proton,MainInput)
    % Inputs:
    %   defectMap (background = 0, lung mask = 1, defect mask = 2)    
    %   voxel_size: [x,y,slice thickness], ex: [3,3,15]
    % Outputs:
    %   D, A, R, DR, C, DDI, mean_DDI
    %   
    % Example:
    %   Package: 129Xe Analysis App
    %
    %   Author: Abdullah Bdaiwi
    %   Work email: abdullah.bdaiwi@cchmc.org
    %   Personal email: abdaiwi89@gmail.com
    %   Website: https://www.cincinnatichildrens.org/research/divisions/c/cpir
    %
    % this code is based on this paper
    % https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9292253/pdf/MRM-86-3224.pdf
    %=========================================================================

    disp('processing 3D DDI, please wait...')
    % Defect Map Size
    defectMap = Ventilation.defectMap_forDDI;
    [dimX, dimY, dimZ] = size(defectMap);
    
    T = 0.5; % sphere radius threshold
    
    % Initialize variables
    meanSlice_DDI = zeros(1, dimZ); 
    stdSlice_DDI = zeros(1, dimZ);
    D = zeros(dimX, dimY, dimZ); % number of defect voxels inside the sphere
    A = zeros(dimX, dimY, dimZ); % number of lung voxels inside the sphere
    DR = zeros(dimX, dimY, dimZ); % defect ratio
    C = zeros(dimX, dimY, dimZ); % cluster score
    DDI = zeros(dimX, dimY, dimZ); % Defect distribution index
    R = zeros(dimX, dimY, dimZ); % sphere radius
    lungmap = double(defectMap > 0);    
    Alung = sum(lungmap(:)) / dimZ;
    voxel_size = [Ventilation.DDI3Dx,Ventilation.DDI3Dy,Ventilation.DDI3Dz];

    z_dim = ceil(voxel_size(3)/voxel_size(1));

    % Get the indices of voxels with value 2 (lung defects)
    defectMap3D = zeros(size(defectMap,1),size(defectMap,2),size(defectMap,3)*z_dim);
    z_dim_med = median(1:1:z_dim);
    shift = z_dim;
    defectMap3D(:,:,z_dim_med) = defectMap(:,:,1);

    for k = 1:size(defectMap,3)
        if k == 1
            defectMap3D(:,:,z_dim_med) = defectMap(:,:,k);
        else
            defectMap3D(:,:,z_dim_med + shift) = defectMap(:,:,k);
            shift = shift + z_dim;
        end
    end

    
    [defectIndicesX, defectIndicesY, defectIndicesZ] = ind2sub(size(defectMap3D), find(defectMap3D == 2));
    defectMap3D = zeros(size(defectMap,1),size(defectMap,2),size(defectMap,3)*z_dim);
    original_z = zeros(1,size(defectMap3D,3));
    shift = 0;
    
    for k = 1:size(defectMap,3)
        defectMap3D(:,:,shift + 1:z_dim + shift) = repmat(defectMap(:,:,k), [1, 1, z_dim]);
        original_z(shift + 1:z_dim + shift) = k;
        shift = shift + z_dim;
    end
       
%     figure; imslice(defectMap3D)    

    % Initialize temporary variables 
    Rtemp = zeros(1,length(defectIndicesX));
    DDItemp= zeros(1,length(defectIndicesX));
    Dtemp= zeros(1,length(defectIndicesX));
    Atemp= zeros(1,length(defectIndicesX));
    DRtemp= zeros(1,length(defectIndicesX));
    Ctemp = zeros(1,length(defectIndicesX));
    
    delete(gcp('nocreate'));
    maxWorkers = maxNumCompThreads;
    parpool('local', maxWorkers); % Change this depending on your CPU

    parfor i = 1:length(defectIndicesX)
        x = defectIndicesX(i);
        y = defectIndicesY(i);
        z = defectIndicesZ(i);
%         disp(['processing slice= ',num2str(z),' defect=',num2str(i),'/',num2str(length(defectIndicesX))])
        
        r = 0; radius = 0; v = 0; d = 0; c = 0; dr = 0; n = 0; volume = 0;
        rad_incrm = 0.5;

        while true
            r = r + 1;
            radius = radius + rad_incrm;
            sphere_mask = sphereMask(defectMap3D, x, y, z, radius);
            v(r) = sum(sphere_mask(:));
            defect_inside_sphere = sphere_mask .* double(defectMap3D == 2);
            d(r) = sum(defect_inside_sphere(:));
            dr(r) = abs(d(r)) / abs(v(r));
            c(r) = ((dr(r) - T) / (1 - T)) * 100;
            n(r) = abs(v(r)) / Alung;
            if c(r) < 0 
                Rtemp(i) = radius - rad_incrm;
               % R(x, y, original_z(z)) = radius - rad_incrm;
                cc = squeeze(c(1:end-1));
                vv = squeeze(v(1:end-1));
                vn = squeeze(abs(vv) ./ Alung);
                p = 0;

                for t = 1:length(cc)
                    width = vn(t) - p;
                    volume(t) = width * cc(t);
                    p = vn(t);
                end

                if round(sum(volume), 1) > 0
                    DDItemp(i) = round(sum(volume), 1);
                    Dtemp(i) = d(end - 1);
                    Atemp(i) = v(end - 1);
                    DRtemp(i) = dr(end - 1);
                    Ctemp(i) = c(end - 1);
                end
                break;
            end
        end
    end

    for i = 1:length(defectIndicesX)
        x = defectIndicesX(i);
        y = defectIndicesY(i);
        z = defectIndicesZ(i);
        
        R(x, y, original_z(z)) = Rtemp(i);

        DDI(x, y, original_z(z)) = DDItemp(i);
        D(x, y, original_z(z)) = Dtemp(i);
        A(x, y, original_z(z)) = Atemp(i);
        DR(x, y, original_z(z)) = DRtemp(i);
        C(x, y, original_z(z)) = Ctemp(i);
    end

    parfor j = 1:dimZ
        ddi_slice = squeeze(DDI(:, :, j));
        meanSlice_DDI(j) = mean(ddi_slice(ddi_slice ~= 0));
        stdSlice_DDI(j) = std(ddi_slice(ddi_slice ~= 0));
    end
    delete(gcp('nocreate'))
    DDIStat = regionprops3(DDI > 0.01);
    Ventilation.DDI3D_DDImap = DDI;
    Ventilation.DDI3D_meanSlice = meanSlice_DDI;
    Ventilation.DDI3D_stdSlice = stdSlice_DDI;
    Ventilation.DDI3D_mean = mean(DDI(DDI ~=0));
    Ventilation.DDI3D_std = std(DDI(DDI ~=0)); 
    Ventilation.DDI3D_max = max(DDI(:));
    Ventilation.DDI3D_Stat = DDIStat;
        %% %% write tiff and read back BinnedVent maps
    DDI_outputpath = [Ventilation.outputpath, '\VDP Analysis\'];
    mkdir(DDI_outputpath);
    cd(DDI_outputpath);
    tiff = figure('MenuBar','none','ToolBar','none','DockControls','off','Resize','off','WindowState','minimized');%figure for tiffs
    ax1 = axes('Parent',tiff);ax2 = axes('Parent',tiff);%make axis for both images
    set(ax1,'Visible','off');set(ax2,'Visible','off');%turn off axis
    set(ax1,'units','inches');set(ax2,'units','inches');%make axis units inches
    set(ax1,'position',[0 0 2 2]);set(ax2,'position',[0 0 2 2]);%make axis same as image
    set(gcf,'units','inches'); % set the figure units to pixels
    set(gcf,'position',[1 1 2 2])% set the position of the figure to axes
    disp('Saving Vent Tiff...')
    %Vent Binned
    for slice=1:size(Ventilation.Image,3) %repeat for rest of slices
        [~,~] = Global.imoverlay(squeeze(abs(Ventilation.Image(:,:,slice))),squeeze(DDI(:,:,slice)),[1,6],[0,0.9*max(Ventilation.Image(:))],colormap('parula'),1,gca);
        colormap(gca,colormap('parula')); clim([0 25]);
         Xdata = getframe(gcf);
         X = Xdata.cdata;     
        if (slice == 1)
            imwrite(X,[DDI_outputpath,'DDI3Dmap.tif'],'Description',strcat('Package Version: ', '1','; Cohort: ', 'test'));%write new/ overwrite tiff
        else
            imwrite(X,[DDI_outputpath,'DDI3Dmap.tif'],'WriteMode','append','Description',strcat('Package Version: ', '1','; Cohort: ', 'test'));%append tiff
        end
    end
    disp('Saving DDI3Dmap Tiff Completed.')
    close all;
    % read tiff
    cd(DDI_outputpath)
    tiff_info = imfinfo('DDI2Dmap.tif'); % return tiff structure, one element per image
    % tiff_stack = imread('BinnedVent.tif', 1) ; % read in first image
    DDI3Dmap = uint8(zeros(tiff_info(1).Height ,tiff_info(1).Width ,3,length(tiff_info)));
    %concatenate each successive tiff to tiff_stack
    for ii = 2 : size(tiff_info, 1)
        temp_tiff = imread('DDI2Dmap.tif', ii);
        DDI3Dmap(:,:,:,ii) = temp_tiff;
    end
    Ventilation.DDI3Dmap = DDI3Dmap;
%     DDI3Dmaptoview = permute(DDI3Dmap,[1 2 4 3]);
%     S = orthosliceViewer((DDI3Dmaptoview)); 
    %% save result in a powerpoint file // Abdullah 1/27/2024
    %% 
    close all;
    NumSliceView = 16;
        % check images size, 2D or 3D
    if size(Ventilation.Image,3) > 30
        Image_3D = 1;
        nrow = ceil(sqrt(size(Ventilation.Image,3)));
        % Initialize an empty vector to store indices of slices with no zeros
        slices_with_no_zeros = zeros(1,size(Ventilation.LungMask,3));    
        % Loop through each slice
        for slice_idx = 1:size(Ventilation.LungMask,3)
            % Extract the current slice
            current_slice = Ventilation.LungMask(:, :, slice_idx);
            
            % Check if all elements in the current slice are non-zero
            if sum(current_slice(:)) > 0
                % Append the slice index to the vector if all elements are non-zero
                slices_with_no_zeros(slice_idx) = 1;
            end
        end
        nonZeroSlices = find(slices_with_no_zeros);
        nonZeroSlice_lng = nonZeroSlices(end) - nonZeroSlices(1);
        nonZeroSlice_space = floor(nonZeroSlice_lng/16);
        slice_indices = nonZeroSlices(1):nonZeroSlice_space:nonZeroSlices(end);
        DDI3Dmap_SS = DDI3Dmap(:,:,:,slice_indices);
    else
        DDI3Dmap_SS = DDI3Dmap;
    end 
    DDIMontage = figure('Name','Vent Image');set(DDIMontage,'WindowState','minimized');
    montage(DDI3Dmap_SS, 'Size',[1 size(DDI3Dmap_SS,4)],'DisplayRange',[]);
    set(gca,'units','pixels'); % set the axes units to pixels
    x = get(gca,'position'); % get the position of the axes
    set(gcf,'units','pixels'); % set the figure units to pixels
    y = get(gcf,'position'); % get the figure position
    set(gcf,'position',[y(1) y(2) x(3) x(4)])% set the position of the figure to the length and width of the axes
    set(gca,'units','normalized','position',[0 0 1 1]) % set the axes units to pixels
    DDIMontagePosition=get(gcf,'position'); 

    parentPath = Ventilation.outputpath;
    cd(parentPath)
    % save ppt 
    ReportTitle = 'Ventilation_Analysis';
    %Start new presentation
    isOpen  = Global.exportToPPTX();
    if ~isempty(isOpen) %If PowerPoint already started, then close first and then open a new one
        Global.exportToPPTX('close');
    end
    ppt_file_name = 'Ventilation_Analysis.pptx';
    if isfile(ppt_file_name)           
        disp('file existed');
        Global.exportToPPTX('open',ppt_file_name);
        Global.exportToPPTX('switchslide',1);
    else            
        Global.exportToPPTX('new','Dimensions',[16 9], ...
            'Title',ReportTitle, ...
            'Author','CPIR @ CCHMC');
    end 
    %Add slides
    Global.exportToPPTX('addslide'); % Image/mask/VDP
    Global.exportToPPTX('addtext',['DDI (3D) map for ', Ventilation.DDIDefectMap, ' method'],'Position',[5 3 5 1],'Color','b','FontSize',25);
    Global.exportToPPTX('addpicture',DDIMontage,'Position',...
        [0 3.5 NumSliceView NumSliceView*(DDIMontagePosition(4)/DDIMontagePosition(3))]);

    meanDDI_title = sprintf('Mean = %0.1f%%±%0.1f%%, Max = %0.1f%%',Ventilation.DDI3D_mean,Ventilation.DDI3D_std,Ventilation.DDI3D_max);
    Global.exportToPPTX('addtext',meanDDI_title,'Position',[0 3.5 5 0.5],'Color','r','FontSize',20,'BackgroundColor','k');
    
    Global.exportToPPTX('save',fullfile(parentPath, ReportTitle));
    Global.exportToPPTX('close');
    fprintf('PowerPoint file has been saved\n');               
    
    close all;
    disp('processing 3D DDI completed...')

%     % create a montage
%     VentMontage = figure('Name','DDI map'); set(VentMontage,'WindowState','minimized');
%     montage(reshape(DDI,[size(DDI,1), size(DDI,2), 1, size(DDI,3)]),...
%         'Size',[1 size(DDI,3)],'DisplayRange',[0 max(DDI(:))]); colormap("jet");
%     set(gca,'units','pixels'); % set the axes units to pixels
%     x = get(gca,'position'); % get the position of the axes
%     set(gcf,'units','pixels'); % set the figure units to pixels
%     y = get(gcf,'position'); % get the figure position
%     set(gcf,'position',[y(1) y(2) x(3) x(4)])% set the position of the figure to the length and width of the axes
%     set(gca,'units','normalized','position',[0 0 1 1]) % set the axes units to pixels
end

function sphere_mask = sphereMask(defectmap, centerX, centerY, centerZ, radius)
    % Inputs:
    %     defectMap (lung mask = 1, defect mask = 2)
    %     centerX 
    %     centerY
    %     centerZ
    %     radius
    % Outputs:
    %   sphere_mask
    % Example:
    %   Package:
    %
    %   Author: Abdullah Bdaiwi
    %   Work email: abdullah.bdaiwi@cchmc.org
    %   Personal email: abdaiwi89@gmail.com
    %   Website: https://www.cincinnatichildrens.org/research/divisions/c/cpir
    %=========================================================================
    % Define array size
    [arraySizeX, arraySizeY, arraySizeZ] = size(defectmap);        
    % Create an empty array
    sphere_mask = zeros(arraySizeX, arraySizeY, arraySizeZ);

    % Create a meshgrid
    [X, Y, Z] = meshgrid(1:arraySizeY, 1:arraySizeX, 1:arraySizeZ);

    % Calculate the distance of each point from the sphere center
    distance = sqrt((X - centerY).^2 + (Y - centerX).^2 + (Z - centerZ).^2);

    % Create a binary mask for the sphere
    sphereMask = distance <= radius;
    % Set the sphere volume to 1 inside the array
    sphere_mask(sphereMask) = 1;
end