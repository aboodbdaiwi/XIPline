function clustering_3D(He, He_clusters)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function was designed to segment ventilated areas in hyperpolarized 
% helium-3 (3He) images, using K-means clustering. ***3D SEGMENTATION***

heSlices = metaImageRead(He);
heSlices (isnan(heSlices)) = 0;
hdr = metaImageInfo(He);

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
    removeSize = 3;
    sideColumns = 10;
elseif width <= 256
    removeSize = 6;
    sideColumns = 15;
else
    removeSize = 12;
    sideColumns = 25;
end
for slice = 1:numberOfslices
    sliceIntensityMapping(slice,1) = max(max(heSlices(:,:,slice)));
end

%% First step of Hirarchical clustering
heSlices = round(((heSlices - min(heSlices(:)))/(max(heSlices(:))-min(heSlices(:))))*255);
[garbage,clusters_3D]=kmeans3D(heSlices,k);
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
        temp2 = miprmdebrisn(temp,removeSize);
        [garbage,col]=find(temp2);
        minCol = min(col);
        maxCol = max(col);
        % Delete the current clustre from ventilated area if the
        % cluster was scattered to the end left and right side
        if (minCol<sideColumns) & (maxCol>(width-sideColumns)) 
            mask(temp) = 0;
        end
    end
    % primaryMask is to show the difference after finding hypoVV
    primaryMask(:,:,slice) = miprmdebrisn(mask,removeSize); %removing the areas smaller than removeSize from ventilated area (assume to be noise)
    heSlicesMask_3D(:,:,slice) = mask;
end

%% Second step of Hirarchical clustering (Clustering HypoVV)
backGroundMask = ~primaryMask; % to find the nonventilated area and the background
J = heSlices(backGroundMask); % J is a vector contains t9he intensity of pixelclustersNews which were segmented asthe background at the first step of H-Kmeans
[garbage,clustersNew]=kmeansModified(J,k2);
temp = clustersNew == k2-1;
secondLastClusterMin = round(min(J(temp)));
secondLastClusterMax = round(max(J(temp)));
temp = clustersNew == k2;
if (sum(temp(:)))>0
     lastClusterMax = round(max(max(J(temp))));
else
     lastClusterMax = secondLastClusterMax;
end

threshold = secondLastClusterMin + ((lastClusterMax-secondLastClusterMin)/7);

secondLastClusterMin
secondLastClusterMax
lastClusterMin = round(min(min(J(temp))))
lastClusterMax = round(max(max(J(temp))))

HypoVent1 = (heSlices>threshold)  & (heSlices<=lastClusterMax);
% to remove a part of current segmented area (as HypoVV) if it was
% selected as VV at previous step. Usually is empty.
HypoVent2 = xor(HypoVent1,primaryMask); 
HypoVent = HypoVent1 & HypoVent2;
%if sum(mask(:))>0  % we assume that if there is not any VV in the current He slice, so no hypoVV as well in this slice
clusters_3D(HypoVent) = 1; % assining hypoVV to be named cluster 1
heSlicesMask_3D(HypoVent) = 1; % Adding the HypoVV to the previously segmented VV
%end

%% Removing small areas from 3D mask (slice by slice)
for slice = 1:numberOfslices
    mask = heSlicesMask_3D(:,:,slice);
    if (slice==1) | (slice==numberOfslices) | (slice==numberOfslices-1)
        mask = miprmdebrisn(mask,removeSize*5); %removing the areas smaller than removeSize from ventilated area (assume to be noise)
    else
        mask = miprmdebrisn(mask,removeSize); %removing the areas smaller than removeSize from ventilated area (assume to be noise)
    end
    negativeMask1 = ~mask;
    negativeMask2 = miprmdebrisn(negativeMask1,removeSize); %removing the areas smaller than removeSize from non-ventilated (assume to be noise)
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

metaImageWrite(clusters_3D, He_clusters, hdr);
end %function

