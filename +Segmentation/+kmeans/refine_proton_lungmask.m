function [mask3D,contour_struct] = refine_proton_lungmask(proton,proton_mask)

[nRows,nCols,nSlices]=size(proton_mask);
lungPixelList=squeeze(sum(sum(proton_mask,1),2));
lungRight_sl=find(lungPixelList).';
mask3D=zeros(nRows,nCols,nSlices);
contour_struct=struct;
nClusters=4;
% se2=strel('disk',2);
spacing=10;
for nsl=lungRight_sl
%     clst=kmeans_clustering(proton(:,:,nsl),nClusters);
%     clst(clst>1)=0;
    slice_proton=proton_mask(:,:,nsl);
    [labeled_mask,nObjs]=bwlabel(slice_proton);
    for nob=1:nObjs
        this_piece=zeros(nRows,nCols);
        this_piece(labeled_mask==nob)=1;
        
        this_piece=imfill(this_piece,'holes');
        contours=mask2poly(this_piece,'CW');
        inValids=contours(:,1)<1 | contours(:,2)>nRows | contours(:,2)<1 | contours(:,2)>nCols;
        contours(inValids,:)=[];
        sampled_points=contours(1:spacing:end,:);
        [mask2d,livewire_contours]=optimize_contours_using_livewire(proton(:,:,nsl),sampled_points);
        
        mask3D(:,:,nsl)=mask3D(:,:,nsl)+mask2d*nob;
        slice_name=sprintf('sl%d',nsl);
        contour_struct.(slice_name)=livewire_contours;
    end
end
