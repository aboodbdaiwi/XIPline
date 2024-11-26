function overlap_vec = compare_mask_overlap(ventmap,mu_map)
% calculate Dice coefficient of each cluster from the two clustered masks;
% W. Zha 10/21/2016
%%
[nRows,nCols,nSlices] = size(ventmap);
cluster_list = union(unique(ventmap),unique(mu_map));
cluster_list(cluster_list == 0) =[];
overlap_vec = zeros(length(cluster_list),1);

for nc = 1:length(cluster_list)
    vent_nc = zeros(nRows,nCols,nSlices);
    mu_nc = zeros(nRows,nCols,nSlices);
    vent_nc(ventmap==cluster_list(nc)) = 1;
    mu_nc(mu_map==cluster_list(nc))=1;
    if sum(vent_nc(:))==0 && sum(mu_nc(:))==0
        overlap_vec(nc)=1;
    elseif (sum(vent_nc(:))==0 && sum(mu_nc(:))~=0) || (sum(vent_nc(:))~=0 && sum(mu_nc(:))==0)
        overlap_vec(nc)=0;
    else
        overlap_mask = vent_nc & mu_nc;
        union_mask = vent_nc + mu_nc;
        overlap_vec(nc) = 2*sum(overlap_mask(:))/sum(union_mask(:));
    end
end