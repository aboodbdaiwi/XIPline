function [Ventilation] = clustering_3DASB(Ventilation)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function was designed to segment ventilated areas in hyperpolarized 
% helium-3 (3He) images, using K-means clustering. ***3D SEGMENTATION***

% heSlices = Ventilation.Image;
heSlices = Ventilation.Image.*Ventilation.LungMask;
heSlices (isnan(heSlices)) = 0;
% lungmask = Ventilation.LungMask;
%% initialization
% tic
k = 4; %number of clusters for clustering whole 3He MRI 
k2 = 4; %number of clusters for clustering Ventilation Defect Volume
heSlicesMask_3D = false(size(heSlices));
primaryMask = false(size(heSlices));
numberOfslices = size(heSlices,3);
clusterMin_3D = zeros(numberOfslices,k);
clusterArea_3D = zeros(numberOfslices,k);
clusterMin3D_All = zeros(1,k);
sliceIntensityMapping = zeros(numberOfslices,2);
width = size(heSlices,2);
height = size(heSlices,1);
if width <= 128
    if strcmp(Ventilation.SliceOrientation,'transversal') 
        removeSize = 2;  
        sideColumns = 5; 
    else
        removeSize = 3;  
        sideColumns = 10; 
    end
elseif width <= 256
    if strcmp(Ventilation.SliceOrientation,'transversal') 
        removeSize = 6; % def 6
        sideColumns = 10; % def 15 %  10 for axial slices works better 
    else
        removeSize = 6; % def 6
        sideColumns = 15; % def 15 %  12 for axial slices works better 
    end    
else
    removeSize = 12;
    sideColumns = 25;
end
for slice = 1:numberOfslices
    sliceIntensityMapping(slice,1) = max(max(heSlices(:,:,slice)));
end

%% First step of Hirarchical clustering
clc;
heSlices = round(((heSlices - min(heSlices(:)))/(max(heSlices(:))-min(heSlices(:))))*255);
[~,clusters_3D] = VentilationFunctions.kmeans.kmeans3D(heSlices,k); % kmeansModified kmeans3D
% [garbage,clusters_3D] = VentilationFunctions.kmeans.kmeansModified(heSlices.*lungmask,k); % kmeansModified kmeans3D
% removeSize = 6;
% sideColumns = 10;
% imslice(clusters_3D)
%scroll_matrix(clusters_3D)
for indx =1:k
    temp = clusters_3D == indx;
    if sum(temp(:)) == 0
        clusterMin3D_All(indx) = 0;
    else
        clusterMin3D_All(indx) = round(min(min(min(heSlices(temp)))));
    end
end

for slice = 1:numberOfslices
    disp(['----------------------------------------> ',num2str(slice)])
    mask = true(height,width);
    I = heSlices(:,:,slice);
    sliceIntensityMapping(slice,2) = max(I(:));
    clusterI = clusters_3D(:,:,slice);
    for indx =1:k
        temp = clusterI == indx;
        clusterArea_3D(slice,indx) = sum(temp(:));
        if clusterArea_3D(slice,indx) == 0
            clusterMin_3D(slice,indx) = 0;
        else
            clusterMin_3D(slice,indx) = round(min(min(I(temp))));
        end
        %removing areas of current cluster which their size is less or
        %equal to removeSize pixels in order to remove noise
        temp2 = VentilationFunctions.kmeans.miprmdebrisn(temp,removeSize);
        [~,col]=find(temp2);
        minCol = min(col);
        maxCol = max(col);
        % Delete the current clustre from ventilated area if the
        % cluster was scattered to the end left and right side
        if (minCol<sideColumns) & (maxCol>(width-sideColumns)) 
            mask(temp) = 0;
        end
        % figure; imslice(I)
    end
    
    % primaryMask is to show the difference after finding hypoVV
    primaryMask(:,:,slice) = VentilationFunctions.kmeans.miprmdebrisn(mask,removeSize); %removing the areas smaller than removeSize from ventilated area (assume to be noise)
    heSlicesMask_3D(:,:,slice) = mask;
end
% primaryMask = Ventilation.LungMask;
% imslice(primaryMask)
%% Second step of Hirarchical clustering (Clustering HypoVV)
% airwaymask = double(Segmentation.SegmentLungthresh(Ventilation.Image,1,0.5)); 
% backGroundMask = logical(imcomplement(primaryMask).*~airwaymask);

backGroundMask = ~primaryMask; % to find the nonventilated area and the background
J = heSlices(backGroundMask); % J is a vector contains t9he intensity of pixelclustersNews which were segmented asthe background at the first step of H-Kmeans
[garbage,clustersNew]= VentilationFunctions.kmeans.kmeansModified(J,k2);
% imslice(clustersNew)
temp = clustersNew == k2-1;
secondLastClusterMin = round(min(J(temp)));
secondLastClusterMax = round(max(J(temp)));
temp = clustersNew == k2;
if (sum(temp(:)))>0
     lastClusterMax = round(max(max(J(temp))));
else
     lastClusterMax = secondLastClusterMax;
end

threshold = secondLastClusterMin + ((lastClusterMax-secondLastClusterMin)/20);

% secondLastClusterMin
% secondLastClusterMax
lastClusterMin = round(min(min(J(temp))))
lastClusterMax = round(max(max(J(temp))))

HypoVent1 = (heSlices>threshold)  & (heSlices<=lastClusterMax);
% to remove a part of current segmented area (as HypoVV) if it was
% selected as VV at previous step. Usually is empty.
HypoVent2 = xor(HypoVent1,primaryMask); 
HypoVent = HypoVent1 & HypoVent2;
%if sum(mask(:))>0  % we assume that if there is not any VV in the current He slice, so no hypoVV as well in this slice
clusters_3D(HypoVent) = 1; % assining hypoVV to be named cluster 1
% imslice(clusters_3D)
heSlicesMask_3D(HypoVent) = 1; % Adding the HypoVV to the previously segmented VV
%end
% imslice(heSlicesMask_3D)
%% Removing small areas from 3D mask (slice by slice)
for slice = 1:numberOfslices
    mask = heSlicesMask_3D(:,:,slice);
    if (slice==1) | (slice==numberOfslices) | (slice==numberOfslices-1)
        mask = VentilationFunctions.kmeans.miprmdebrisn(mask,removeSize*5); %removing the areas smaller than removeSize from ventilated area (assume to be noise)
    else
        mask = VentilationFunctions.kmeans.miprmdebrisn(mask,removeSize); %removing the areas smaller than removeSize from ventilated area (assume to be noise)
    end
    negativeMask1 = ~mask;
    negativeMask2 = VentilationFunctions.kmeans.miprmdebrisn(negativeMask1,removeSize); %removing the areas smaller than removeSize from non-ventilated (assume to be noise)
    clusterI = clusters_3D(:,:,slice);
    clusterI(negativeMask2) = 0; % assining all non-ventilated area to the background cluster, named cluster 0
    mask = ~negativeMask2;
    temp3 = xor(negativeMask1,negativeMask2);
    clusterI(temp3) = 1; % assining small non-ventilated area to the hypoVV cluster, named cluster 1
    temp = clusterI == 1;
    clusterArea_3D(slice,1) = sum(temp(:));
    clusters_3D(:,:,slice) = clusterI;
    heSlicesMask_3D(:,:,slice) = mask;
end
Ventilation.clusters_3D = clusters_3D;
% imslice(clusters_3D)

end %function

