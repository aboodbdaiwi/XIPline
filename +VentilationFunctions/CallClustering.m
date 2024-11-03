function [heSlicesClusters,heSlicesMask,clusterMin] = CallClustering(heSlices)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function CallClustering_3He(heSlices)
% This function segments Ventilated Volume (VV) in hyperpolarized 
% helium-3 (3He) images, using K-means clustering.
% 
% input: 
%         heSlices: All 3He slices in a 3D Matlab matrix
% 
% outputs:
%    The results will be saved in a MAT-file (HeClusteringResult.mat)
%    including following data: 
%         heSlicesClusters: cluster information for all slices 
%         heSlicesMask: mask information for all slices 
%         clusterMin: minimum intensity of pixels in each cluster for all
%               slices
% 
% Author: M.Reza Heydarian
% lab: Robarts Research Institute, London, ON, Canada
% http://www.imaging.robarts.ca/~gep/
% Date: October 30, 2009 
% Modified: Aug. 2010, June 2010 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%% initialization
tic
k = 5; %number of clusters for clustering whole 3He MRI 
k2 = 4; %number of clusters for clustering Ventilation Defect Volume
heSlicesClusters = zeros(size(heSlices),'uint8');
heSlicesMask = false(size(heSlices));
clusterMin = zeros(size(heSlices,3),k);
width = size(heSlices,2);
height = size(heSlices,1);
numberOfslices = size(heSlices,3);
if width <= 128
    removeSize = 3;
    sideColumns = 10;
elseif width <= 256
    removeSize = 7;
    sideColumns = 15;
else
    removeSize = 12;
    sideColumns = 25;
end

%% Clustering 3He MRI
for slice = 1:numberOfslices 
    %--------------------->>>  First step of Hierarchical Clustering
    %disp(['------------------------------------------- Slice# ',num2str(slice)])
    I = heSlices(:,:,slice);
    sumOfPixels = sum(I(:));
    if sumOfPixels ~= 0
        % Nov 30 2011 awheat changes
        %I = round(((I - min(I(:)))/(max(I(:))-min(I(:))))*255); % normalized image to [0 255]
        I = double(I);
        I = I - min(I(:));
        I = I / max(I(:));
        I = I * 255;
        I = round(I);
        
        [garbage,clusters] = VentilationFunctions.kmeansModified(I,k);
        %figure; imshow(clusters,[]), title(['Four clusters, result of "',clusteringMethod,'"'])
        mask = true(size(clusters));
        for indx =1:k
            temp = clusters == indx;
            if sum(temp(:)) == 0
                clusterMin(slice,indx) = 0;
            else
                clusterMin(slice,indx) = round(min(min(I(temp))));
            end
            %removing areas of current cluster which their size is less or
            %equal to removeSize pixels in order to remove noise
            temp = VentilationFunctions.miprmdebrisn(temp,removeSize);
            [garbage,col]=find(temp);
            minCol = min(col);
            maxCol = max(col);
            % Delete the current clustre from the ventilated area if the
            % cluster was scattered to the end left and the end right sides
            if ((minCol<sideColumns) & (maxCol>width-sideColumns))
                mask(temp) = 0;
            end
        end
        %primaryMask is to show the difference after finding hypoVV
        %primaryMask = miprmdebrisn(mask,removeSize); %removing the areas smaller than removeSize from ventilated area (assume to be noise)
        %figure; imshow(I,[]), hold on
        %contour(primaryMask,[0 0],'r');
        %title('Primarily mask, result of combining clusters')

        %--------------------->>>  Second step of Hierarchical Clustering
        J = I(~mask); % to find the nonventilated area and the background
        [~,clustersNew]= VentilationFunctions.kmeansModified(J,k2);
        temp = clustersNew == k2-1;
        secondLastClusterMin = round(min(min(J(temp))));
        ct = k2 -2; 
        while isempty(secondLastClusterMin)
            temp = clustersNew == ct;
            secondLastClusterMin = round(min(min(J(temp))));
            ct = ct - 1; 
        end    
        temp = clustersNew == k2;
        lastClusterMax = round(max(max(J(temp))));
        ct = k2 - 1; 
        while isempty(lastClusterMax)
            temp = clustersNew == ct;
            lastClusterMax = round(max(max(J(temp))));
            ct = ct - 1;
        end   
        
        threshold = secondLastClusterMin + round((lastClusterMax-secondLastClusterMin)/7);
        HypoVent1 = (I>threshold)  & (I<=lastClusterMax);
        %figure; imshow(HypoVent1,[]), title(['Slice = ',num2str(slice),'         Cluster = Hypo Ventilated Area'])

        % removing a part of current segmented area (as HypoVV) if it was
        % selected as VV at previous step. Usually is empty.
        HypoVent2 = xor(HypoVent1,mask); 
        %figure; imshow(HypoVent2,[]), title(['Slice = ',num2str(slice),'         Cluster = xor(HypoVent,mask)'])
        HypoVent = HypoVent1 & HypoVent2;
        %figure; imshow(HypoVent,[]), title(['Slice = ',num2str(slice),'         Cluster = Hypo Ventilated Area'])
        if sum(mask(:))>0  % we assume that if there is not any VV in the current He slice, so there is no hypoVV in this slice.
            clusters(HypoVent) = 1; % assining hypoVV to be named cluster 1
            mask(HypoVent) = 1; % Adding the HypoVV to the previously segmented VV
        end
        if ((slice==1) | (slice==numberOfslices) | (slice==numberOfslices-1))
            mask = VentilationFunctions.miprmdebrisn(mask,removeSize*5); %removing the areas smaller than removeSize*5 from ventilated area (assume to be noise)
        else
            mask = VentilationFunctions.miprmdebrisn(mask,removeSize); %removing the areas smaller than removeSize from ventilated area (assume to be noise)
        end
        negativeMask1 = ~mask;
        negativeMask2 = VentilationFunctions.miprmdebrisn(negativeMask1,removeSize); %removing the areas smaller than removeSize from non-ventilated (assume to be noise)
        mask = ~negativeMask2;
        temp3 = xor(negativeMask1,negativeMask2);
        clusters(negativeMask2) = 0; % assining all non-ventilated area to the background cluster, named cluster 0
        clusters(temp3) = 1; % assining small non-ventilated area to the hypoVV cluster, named cluster 1
        %figure; imshow(HypoVent,[]), title(['Slice = ',num2str(slice),'         Cluster = Hypo Ventilated Area'])
        
    else
        clusters = zeros(height,width,'uint8');
        mask = false(height,width);
    end    
    heSlicesClusters(:,:,slice) = clusters;
    heSlicesMask(:,:,slice) = mask;
    %imshow(I,[]), hold on
    %contour(mask,[0 0],'g');
    %%contour(primaryMask,[0 0],'r');
    %title(['Slice# ',num2str(slice),',  VV including hypo-VV.'])
    %%title(['Segmentation after adding hypo-ventilated area from VDV+backgroung mask;  slice# ',num2str(slice)])
end
t = toc; %to measure run time
%disp(['Time for segmenting (clustering) helium images (in second): ',num2str(t)])
% scroll_matrix_contour(heSlices,heSlicesMask)

end %function