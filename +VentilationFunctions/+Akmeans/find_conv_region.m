function [convex_list,concave_list]=find_conv_region(lobe_contours,conv_points,conv_indices)

% Input: lobe_contours: Nx2 matrix of contour points
%        
%
% Output:-convex_list: index list corresponds to the convex portion of the
%                      contour
%        -concave_list: index list corresponds to the concave portion of
%                      the contour
%
% Copyright W. Zha @ 3/30/2015

nPoints=size(lobe_contours,1);
if nargin < 2
    conv_indices = convhull(lobe_contours(:,1),lobe_contours(:,2));
    conv_points = lobe_contours(conv_indices,:);
end
adjacent_dist=sqrt(sum((conv_points(1:(end-1),:)-conv_points(2:end,:)).^2,2));
[~,max_index]=max(adjacent_dist);
vertex_1=conv_indices(max_index);
vertex_2=conv_indices(max_index+1);
partition_list = min(vertex_1,vertex_2):max(vertex_1,vertex_2);
complement_list = [max(vertex_1,vertex_2):nPoints,1:1:min(vertex_1,vertex_2)];

if ~isempty(intersect(conv_indices,partition_list(2:end-1)))
    convex_list = partition_list;
    concave_list=complement_list;
else
    convex_list = complement_list;
    concave_list= partition_list;
end

