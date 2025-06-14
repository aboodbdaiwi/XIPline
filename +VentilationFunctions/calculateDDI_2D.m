
function [Ventilation] = calculateDDI_2D(Ventilation,Proton,MainInput)
%   Inputs: 
%       defectMap (background = 0, lung mask = 1, defect mask = 2)    
%   Outputs:
%       DDI map and mean DDI for each slice
%       D,A,R,DR,C,DDI,mean_DDI
%   Example: 
%   Package: 
%
%   Author: Abdullah Bdaiwi 
%   Work email: abdullah.bdaiwi@cchmc.org
%   Personal email: abdaiwi89@gmail.com
%   Website: https://www.cincinnatichildrens.org/research/divisions/c/cpir
%
% this code is based on this paper
%https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9292253/pdf/MRM-86-3224.pdf
%=========================================================================

    disp('processing 2D DDI, please wait...')
    defectMap = Ventilation.defectMap_forDDI;
    
    % Defect Map Size
    [dimX, dimY, dimZ] = size(defectMap);
    
    T = 0.5; % circal reduis threshold
    
    % Initialize variables
    meanSlice_DDI = zeros(1, dimZ);
    stdSlice_DDI = zeros(1, dimZ);
    D = zeros(dimX, dimY, dimZ);
    A = zeros(dimX, dimY, dimZ);
    DR = zeros(dimX, dimY, dimZ);
    C = zeros(dimX, dimY, dimZ);
    DDI = zeros(dimX, dimY, dimZ);
    R = zeros(dimX, dimY, dimZ);
    lungmap = double(defectMap > 0);
    
    Alung = sum(lungmap(:)) / dimZ;
    
    % Get the indices of voxels with value 2 (lung defects)
    delete(gcp('nocreate'));
    maxWorkers = maxNumCompThreads;
    parpool('local', maxWorkers);
    
    % Preallocate arrays for parfor
    parfor slice = 1:dimZ
        img = defectMap(:, :, slice);
        defectIndices = find(img == 2);
        defectonlymap = double(img > 1);
    
        % Initialize arrays for parfor
        D_slice = zeros(dimX, dimY);
        A_slice = zeros(dimX, dimY);
        DR_slice = zeros(dimX, dimY);
        C_slice = zeros(dimX, dimY);
        DDI_slice = zeros(dimX, dimY);
        R_slice = zeros(dimX, dimY);
    
        % Loop through the defect indices and place circles
        for i = 1:length(defectIndices)
            % Get the current defect voxel index
            currentDefectIndex = defectIndices(i);
    
            % Get the current defect voxel coordinates
            [x, y, ~] = ind2sub(size(defectMap(:, :, slice)), currentDefectIndex);
    
            % Initialize variables
            radius = 0; r = 0; a=0; d=0; c=0; dr=0; n=0; area=0;
    
            % Loop until defect cluster score < 0
            while true
                r = r + 1;
                radius = radius + 0.5;
                circal_mask = circalMask(defectMap, x, y, radius);
                a(r) = sum(circal_mask(:));
                defectinsidecircal = circal_mask .* defectonlymap;
                d(r) = sum(defectinsidecircal(:));
                dr(r) = abs(d(r)) / abs(a(r));
                c(r) = ((dr(r) - T) / (1 - T)) * 100;
                n(r) = abs(a(r)) / Alung;
    
                % Check if cluster score < 0
                if c(r) < 0
                    R_slice(x, y) = radius - 0.5;
                    cc = squeeze(c(1:end - 1));
                    aa = squeeze(a(1:end - 1));
                    nn = squeeze(abs(aa) / Alung);
                    p = 0;
                    for t = 1:length(cc)
                        width = nn(t) - p;
                        area(t) = width .* cc(t);
                        p = nn(t);
                    end
                    if round(sum(area), 1) > 0
                        D_slice(x, y) = d(end - 1);
                        A_slice(x, y) = a(end - 1);
                        DR_slice(x, y) = dr(end - 1);
                        C_slice(x, y) = c(end - 1);
                        DDI_slice(x, y) = round(sum(area), 1);
                    end
                    break;
                end
            end
        end
    
        % Store results for the slice
        D(:, :, slice) = D_slice;
        A(:, :, slice) = A_slice;
        DR(:, :, slice) = DR_slice;
        C(:, :, slice) = C_slice;
        DDI(:, :, slice) = DDI_slice;
        R(:, :, slice) = R_slice;
    
        % Calculate mean and std for DDI for the slice
        ddi_slice = squeeze(DDI(:, :, slice));
        meanSlice_DDI(slice) = mean(ddi_slice(ddi_slice ~= 0));
        stdSlice_DDI(slice) = std(ddi_slice(ddi_slice ~= 0));
    end
    
    % Close the parallel pool
    delete(gcp);
    DDIStat = regionprops3(DDI > 0.01);
    Ventilation.DDI2D_DDImap = DDI;
    Ventilation.DDI2D_meanSlice = meanSlice_DDI;
    Ventilation.DDI2D_stdSlice = stdSlice_DDI;
    Ventilation.DDI2D_mean = mean(DDI(DDI ~=0));
    Ventilation.DDI2D_std = std(DDI(DDI ~=0)); 
    Ventilation.DDI2D_max = max(DDI(:));
    Ventilation.DDI2D_Stat = DDIStat;
    %% %% write tiff and read it back 
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
            imwrite(X,[DDI_outputpath,'DDI2Dmap.tif'],'Description',strcat('Package Version: ', '1','; Cohort: ', 'test'));%write new/ overwrite tiff
        else
            imwrite(X,[DDI_outputpath,'DDI2Dmap.tif'],'WriteMode','append','Description',strcat('Package Version: ', '1','; Cohort: ', 'test'));%append tiff
        end
    end
    disp('Saving DDI2Dmap Tiff Completed.')
    close all;
    % read tiff
    cd(DDI_outputpath)
    tiff_info = imfinfo('DDI2Dmap.tif'); % return tiff structure, one element per image
    % tiff_stack = imread('BinnedVent.tif', 1) ; % read in first image
    DDI2Dmap = uint8(zeros(tiff_info(1).Height ,tiff_info(1).Width ,3,length(tiff_info)));
    %concatenate each successive tiff to tiff_stack
    for ii = 2 : size(tiff_info, 1)
        temp_tiff = imread('DDI2Dmap.tif', ii);
        DDI2Dmap(:,:,:,ii) = temp_tiff;
    end
    Ventilation.DDI2Dmap = DDI2Dmap;
%     DDI2Dmaptoview = permute(DDI2Dmap,[1 2 4 3]);
%     S = orthosliceViewer((DDI2Dmaptoview)); %colormap(SixBinMap);
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
        DDI2Dmap_SS = DDI2Dmap(:,:,:,slice_indices);
    else
        Image_3D = 0;
        DDI2Dmap_SS = DDI2Dmap;
    end 
    DDIMontage = figure('Name','Vent Image');set(DDIMontage,'WindowState','minimized');
    montage(DDI2Dmap_SS, 'Size',[1 size(DDI2Dmap_SS,4)],'DisplayRange',[]);
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

    %Start new presentation
    isOpen  = Global.exportToPPTX(); 
    if ~isempty(isOpen) %If PowerPoint already started, then close first and then open a new one
        Global.exportToPPTX('close');
    end
    % Generate filename with today's date
    today = datetime('today');
    date_str = datestr(today, 'yyyymmdd');
    ppt_file_name = ['Ventilation_Analysis_', date_str, '.pptx'];
    ReportTitle = ['Ventilation_Analysis_', date_str];
    % Create or open the presentation
    if isfile(ppt_file_name)
        disp('File existed')
        Global.exportToPPTX('open', ppt_file_name);
        Global.exportToPPTX('switchslide', 1);
    else            
        Global.exportToPPTX('new', 'Dimensions', [16 9], ...
            'Title', ReportTitle, ...
            'Author', 'CPIR @ CCHMC');
    end 
    
    %Add slides
    Global.exportToPPTX('addslide'); % Image/mask/VDP
    Global.exportToPPTX('addtext',['DDI (2D) map for ', Ventilation.DDIDefectMap, ' method'],'Position',[5 0 5 1],'Color','b','FontSize',25);
    Global.exportToPPTX('addpicture',DDIMontage,'Position',...
        [0 0.5 NumSliceView NumSliceView*(DDIMontagePosition(4)/DDIMontagePosition(3))]);

    meanDDI_title = sprintf('Mean = %0.1f%%±%0.1f%%, Max = %0.1f%%',Ventilation.DDI2D_mean,Ventilation.DDI2D_std,Ventilation.DDI2D_max);
    Global.exportToPPTX('addtext',meanDDI_title,'Position',[0 0.5 5 0.5],'Color','r','FontSize',20,'BackgroundColor','k');
    
    
    Global.exportToPPTX('save',fullfile(parentPath, ReportTitle));
    Global.exportToPPTX('close');
    fprintf('PowerPoint file has been saved\n');               
    
    close all;
end

function circal_mask = circalMask(defectmap,centerY,centerX,radius)
    % Inputs:
    %     defectMap (lung mask = 1, defect mask = 2)
    %     centerX 
    %     centerY
    %     radius
    % Outputs:
    %   circal_mask
    % Example:
    %   Package:
    %
    %   Author: Abdullah Bdaiwi
    %   Work email: abdullah.bdaiwi@cchmc.org
    %   Personal email: abdaiwi89@gmail.com
    %   Website: https://www.cincinnatichildrens.org/research/divisions/c/cpir
    %=========================================================================
    % % Define array size
    arraySizeX = size(defectmap,1);
    arraySizeY = size(defectmap,2);
    % Create an empty array
    circal_mask = zeros(arraySizeX, arraySizeY);
    % Create a meshgrid
    [X, Y] = meshgrid(1:arraySizeY, 1:arraySizeX);
    % Calculate the distance of each point from the circle center
    distance = sqrt((X - centerX).^2 + (Y - centerY).^2);
    % Create a binary mask for the circle
    circleMask = distance <= radius;
    % Set the circle area to 1 inside the array
    circal_mask(circleMask) = 1;
end