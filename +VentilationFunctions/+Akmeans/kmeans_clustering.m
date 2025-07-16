function [clustered_stack,ctr] = kmeans_clustering(image_stack,nClusters)
% if nargin<2
%     nClusters=4;
% end
% KMopts =statset('MaxIter', 300);
% 
% pixel_list = find(image_stack(:)~=0);
% 
% [clustered_vec,ctr]=kmeans(image_stack(pixel_list),nClusters,'distance','sqeuclidean', 'Replicates',1, 'Options', KMopts,'Start','uniform');
% [~,sorted_order] = sort(ctr, 'ascend');
 if nargin < 2
    nClusters = 4;
end

KMopts = statset('MaxIter', 300);

% Get indices of non-zero pixels
pixel_list = find(image_stack(:) ~= 0);
data = image_stack(pixel_list);

clustered_stack = zeros(size(image_stack));  % default output
ctr = zeros(nClusters, 1);  % default cluster centers

% Return early if not enough points for clustering
if numel(data) < nClusters
    warning('Not enough non-zero pixels (%d) for %d clusters. Returning zeros.', numel(data), nClusters);
    return;
end

% Run kmeans
[clustered_vec, ctr] = kmeans(data, nClusters, ...
    'distance', 'sqeuclidean', ...
    'Replicates', 1, ...
    'Options', KMopts, ...
    'Start', 'uniform');
% Optional: sort clusters by ascending centroid value
[~, sorted_order] = sort(ctr, 'ascend');
sorted_clustered_vec = arrayfun(@(x) find(sorted_order == x), clustered_vec);

% Map back to image space
clustered_stack(pixel_list) = sorted_clustered_vec;
%% restore the image size
clustered_stack = zeros(size(image_stack));
for n_clst = 1:nClusters
    indices = (clustered_vec==sorted_order(n_clst));
    masked_index_list = pixel_list(indices);
    clustered_stack(masked_index_list)=n_clst;
end

% figure,imshow3D(clustered_stack);