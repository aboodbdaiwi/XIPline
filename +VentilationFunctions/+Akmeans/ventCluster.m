function [clusterData, clusterSizes,cluster_centroid] = ventCluster2(imageData, maskData,num_clusters,cluster_option)
% Inputs: imageData        -2D or 3D
%         maskData         -optional
%         num_clusters     -almost always be 4
%         cluster_option   -1) take the first 2 sub-clusters of the first
%                              of the intial 4 clusters as cluster 1;
%                          -2) take the first of the initial 4 clusters as cluster 1;
%                          -3) take the first 3 sub-clusters of the first
%                              of the intial 4 clusters as cluster 1;
% Output: clusterData      -clustered data of the same size as the
%                           imageData;
%         clusterSizes     -a vector records the number of pixels
%                           classified into each cluster;
%
% W. Zha @ 2014
%%
debug =false;
if nargin <2
    num_clusters=4;
    maskData=[];
    cluster_option=2;
elseif nargin<3
    num_clusters=4;
    cluster_option=2;
elseif nargin<4
    cluster_option=2;
else
    if ~isempty(maskData) && numel(imageData) ~= numel(maskData)
        error('Error:  Image data and mask must be the same size');
    end
end
[nx,ny,nz] = size(imageData);

if isempty(maskData)
    index_list_masked_lung=1:numel(imageData);
else
    index_list_masked_lung = find(maskData);
end
%%
if ~isempty(maskData)
    validData = double(imageData(logical(maskData)));
else
    validData = double(imageData(:));
end
if isempty(validData)
    clusterData=zeros(nx,ny,nz);
    clusterSizes=[];
    cluster_centroid=0;
    return;
end
%% initial K-means clustering
% normalize to range of [0 1]
if max(validData)~=min(validData)
    validData_norm = (validData - min(validData)) ./ (max(validData) - min(validData));
    if max(validData)==0 && min(validData)==0
        clusterData=zeros(nx,ny,nz);
        clusterSizes =[];
        return;
    end
else
    clusterData=zeros(nx,ny,nz);
    clusterSizes=[];
    return;
end
KMopts = statset('MaxIter', 250);
nCluster1=num_clusters;
[clust1, cent1] = kmeans(validData_norm, nCluster1, 'start','uniform','distance','sqEuclidean', 'Replicates',5, 'Options', KMopts,'emptyaction','drop');

nClusters=4;
if nCluster1==5
    [clust_first_tmp,cent1_tmp] = kmeans(validData_norm, nClusters, 'start','uniform','distance','sqEuclidean', 'Replicates',5, 'Options', KMopts,'emptyaction','drop');
else
    clust_first_tmp = clust1;
    cent1_tmp = cent1;
end
% sort the centroids
[~,sorted_order] = sort(cent1, 'ascend');
[~,sorted_order_tmp] = sort(cent1_tmp, 'ascend');

if debug
    first_cluster = zeros(nx,ny,nz,'uint8');
    
    for n_clst = 1:num_clusters
        indices = (clust1==sorted_order(n_clst));
        masked_index_list = index_list_masked_lung(indices);
        first_cluster(masked_index_list)=n_clst;
    end
    
    figure,imshow(first_cluster,[]);
    title('first round of kmeans clustering');
    dummy = zeros(size(first_cluster));
    dummy(first_cluster==1)=1;
    figure,imshow(dummy,[]);
end
% Researve bin1 and bin2 for the 2nd round of Kmeans while sort all pixels
% above cent1==1 into 3 bins: bin3,bin4 and bin5.
cluster_centroid=cent1(sorted_order(1));
for n_clst=2:nClusters
%     if n_clst>4 % 4 bins only regardless of nCluster1
%         cluster.(this_bin)=[cluster.(this_bin); find(clust_first4==sorted_order(n_clst))];
%     else
        this_bin=sprintf('bin%d',n_clst+1);
        cluster.(this_bin)=find(clust_first_tmp==sorted_order_tmp(n_clst));
%     end
end
clust_first4_bin1 = find(clust_first_tmp==sorted_order_tmp(1));
%% second round of K-means
% find the smallest centroid from first clustering and recluster
index_list_1st_cluster = find(clust1 == sorted_order(1));
validData2 = double(validData(index_list_1st_cluster));
% renormalize
validData2_norm = (validData2 - min(validData2)) ./ (max(validData2) - min(validData2));

% compute second level clustering (again into 4 clusters, same initial
% condition)

final_num_clusters=4;
if sum(isnan(validData2_norm))<length(validData2_norm) && sum(validData2_norm)~=0 && length(validData2_norm)>nClusters
    %     nClusters=num_clusters;
    [clust2, cent2] = kmeans(validData2_norm, nClusters, 'start','uniform','distance','sqEuclidean',  'Options', KMopts,'emptyaction','drop');
    % group clusters with 2 smallest and 2 largest second level centroids
    [~,refined_sorted_order] = sort(cent2, 'ascend');
    if debug
        refined_C_pixels=zeros(nClusters,1);
        for n_clst=1:nClusters
            refined_C_pixels(n_clst) = sum(clust2==refined_sorted_order(n_clst));
            
        end
        disp('# of pixels in the refined cluster: 1\n');
        disp(refined_C_pixels);
    end
    if cluster_option==1
        %         final_num_clusters=num_clusters+1;
        cluster.bin1=[index_list_1st_cluster(clust2 == refined_sorted_order(1)); index_list_1st_cluster(clust2 == refined_sorted_order(2))];
        cluster.bin2=[index_list_1st_cluster(clust2 == refined_sorted_order(3)); index_list_1st_cluster(clust2 == refined_sorted_order(4))];
        cluster_centroid=cluster_centroid*cent2(refined_sorted_order(2));
    elseif cluster_option==2
        %         final_num_clusters=num_clusters;
        cluster.bin1=index_list_1st_cluster;
        cluster.bin2=[];
%         for n_clst=2:num_clusters
%             this_bin = sprintf('bin%d',n_clst);
%             next_bin = sprintf('bin%d',n_clst+1);
%             cluster.(this_bin)=cluster.(next_bin);
%         end
    elseif cluster_option==3
        %         final_num_clusters=num_clusters+1;
        cluster.bin1=[index_list_1st_cluster(clust2 == refined_sorted_order(1)); index_list_1st_cluster(clust2 == refined_sorted_order(2)); ...
            index_list_1st_cluster(clust2 == refined_sorted_order(3))];
        cluster.bin2=index_list_1st_cluster(clust2 == refined_sorted_order(4));
        cluster_centroid=cluster_centroid*cent2(refined_sorted_order(3));
        
    end
else
    %     final_num_clusters=num_clusters;
    cluster.bin1=[];
    cluster.bin2=[];
%     for n_clst=2:num_clusters
%         this_bin = sprintf('bin%d',n_clst);
%         next_bin = sprintf('bin%d',n_clst+1);
%         cluster.(this_bin)=cluster.(next_bin);
%     end
%     last_bin=sprintf('bin%d',num_clusters+1);
%     cluster=rmfield(cluster,last_bin);
end
%% Reshape the clustered pixels into a matrix of the input data size
% convert the indices in id1-5 into indices in to a 3d volume of the same
% size as the original imageData (and maskData)
% clusterData = zeros(nx,ny,nz,'uint8');
clusterData = zeros(nx,ny,nz);
clusterSizes = zeros(1,final_num_clusters);

%% matching
% bin1- cluster#1
% bin2 & bin3 - cluster#2
% bin4 - cluster#3
% bin5 - cluster#4

index_list.C1 = index_list_masked_lung(cluster.bin1);
if nCluster1==5
    lost_pixels = setdiff(clust_first4_bin1,cluster.bin1);
    index_list.C2 = index_list_masked_lung([lost_pixels;cluster.bin2;cluster.bin3]);
else
    index_list.C2 = index_list_masked_lung([cluster.bin2;cluster.bin3]);
end
index_list.C3 = index_list_masked_lung(cluster.bin4);
index_list.C4 = index_list_masked_lung(cluster.bin5);
for n_clst=1:final_num_clusters
    this_label=sprintf('C%d',n_clst);
    clusterData(index_list.(this_label))=n_clst;
    clusterSizes(n_clst)=length(index_list.(this_label));
end


