function [clustered_stack,ctr] = kmeans_clustering(image_stack,nClusters)
if nargin<2
    nClusters=4;
end
KMopts =statset('MaxIter', 300);

pixel_list = find(image_stack(:)~=0);

[clustered_vec,ctr]=kmeans(image_stack(pixel_list),nClusters,'distance','sqeuclidean', 'Replicates',1, 'Options', KMopts,'Start','uniform');
[~,sorted_order] = sort(ctr, 'ascend');

%% restore the image size
clustered_stack = zeros(size(image_stack));
for n_clst = 1:nClusters
    indices = (clustered_vec==sorted_order(n_clst));
    masked_index_list = pixel_list(indices);
    clustered_stack(masked_index_list)=n_clst;
end

% figure,imshow3D(clustered_stack);