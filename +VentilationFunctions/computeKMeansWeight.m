function Weight = computeKMeansWeight(MR, maskarray, nClusters)
% computeKMeansWeight - Generates a normalized weight map using k-means clustering on MR intensities within a mask
%
% Syntax:
%   Weight = computeKMeansWeight(MR, maskarray, nClusters)
%
% Inputs:
%   MR         - Input MR image (2D or 3D)
%   maskarray  - Binary mask of the region of interest (same size as MR)
%   nClusters  - Number of k-means clusters (e.g., 3–5)
%
% Output:
%   Weight     - Normalized weight map (same size as MR)

% Ensure double precision
MR = double(MR);

% Get intensities inside the mask
masked_vals = MR(maskarray);

% Run k-means on masked intensity values
[idx, C] = kmeans(masked_vals, nClusters, 'Replicates', 5);

% Assign weights to clusters based on intensity (higher intensity → higher weight)
[~, sortIdx] = sort(C, 'descend'); % ascend 
clusterWeights = linspace(1, .1, nClusters);      % Adjust scale if needed
clusterWeights = clusterWeights(sortIdx);         % Match cluster index to weights

% Create empty weight map
Weight = zeros(size(MR));

% Fill in weights for masked region
masked_indices = find(maskarray);
for i = 1:nClusters
    cluster_mask = (idx == i);
    Weight(masked_indices(cluster_mask)) = clusterWeights(i);
end

% Normalize weight map to [0, 1]
Weight = Weight ./ max(Weight(:));

end
