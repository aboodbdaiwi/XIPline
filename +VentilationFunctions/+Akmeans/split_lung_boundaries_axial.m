function [convex_boundary,concave_boundary] = split_lung_boundaries_axial(temp_handles,slice_range)
% Input: -temp_handles.lungsmask: binary lung masks
%        -slice_range: the list of slices for boundary split.
%
% Output:-convex_boundary: 3D matrix that contains the edge masks of "outter
%                          boundary" of the lungs
%        -concave_boundary: 3D matrix that contains the edge masks of "inner
%                          boundary" of the lungs where the vessels may be
%                          attached to.
%
% Copyright W. Zha @ 3/30/2015
%%
if nargin<2
    slice_range=1:size(temp_handles.lungsmask,3);
end
lungsmask=temp_handles.lungsmask;
[nRows,nCols,nSlices]=size(lungsmask);
convex_boundary=zeros(nRows,nCols,nSlices);
concave_boundary=zeros(nRows,nCols,nSlices);
[temp_handles,lobe_label_values]=identify_lobes(temp_handles);
%% Save the binary right and left lung masks into lobemask.RUL/LUL
% The possible gaps among RUL, RML and RLL labeled in SARP data may 
% cause fake lung boundaries, so RUL and LUL are automatically 
% relabeled from lung masks here regardless.

if length(lobe_label_values)==5
    RLL=zeros(nRows,nCols,nSlices);
    RML=zeros(nRows,nCols,nSlices);
    LLL=zeros(nRows,nCols,nSlices);
    se=strel('disk',2);
    for nsl=1:nSlices
        RLL(:,:,nsl)=imdilate(temp_handles.lobemaskRLL(:,:,nsl),se);
        RML(:,:,nsl)=imdilate(temp_handles.lobemaskRML(:,:,nsl),se);
        LLL(:,:,nsl)=imdilate(temp_handles.lobemaskLLL(:,:,nsl),se);
    end
    lobemask.RUL=(temp_handles.lobemaskRUL+RML+RLL).*lungsmask;
    lobemask.RUL(lobemask.RUL>0)=1;
    lobemask.LUL=(temp_handles.lobemaskLUL+LLL).*lungsmask;
    lobemask.LUL(lobemask.LUL>0)=1;
elseif length(lobe_label_values)==2
    lobemask.RUL=temp_handles.lobemaskRUL;
    lobemask.LUL=temp_handles.lobemaskLUL;
end
%% separate the convex and concave portion of the lung boundaries.
lobe_list={'RUL','LUL'};
for n_lobe=1:numel(lobe_list)
    lobe_name=lobe_list{n_lobe};
    this_lobe=lobemask.(lobe_name);
    for n_ind=1:length(slice_range)
        n_sl=slice_range(n_ind);
        slice_mask=this_lobe(:,:,n_sl);
        if sum(slice_mask(:))>0
            % the right/left lung masks may be splited toward the base due to the appearance of diaphragm.
            [labeled_mask,nObj]=bwlabel(slice_mask);
            for n=1:nObj
                object_mask=zeros(nRows,nCols);
                object_mask(labeled_mask==n)=1;
                lobe_contours = mask2poly(object_mask,'Exact','CW');
                invalids = (lobe_contours(:,1)<0 | lobe_contours(:,2)<0);
                lobe_contours(invalids,:)=[];
                [convex_list,concave_list]=find_conv_region(lobe_contours);
                
                convexBoundary=lobe_contours(convex_list,:);
                convexboundary_indices = sub2ind([nRows,nCols],convexBoundary(:,2),convexBoundary(:,1));
                convex_edge = zeros(nRows,nCols);
                convex_edge(convexboundary_indices)=1;
                convex_boundary(:,:,n_sl)=convex_boundary(:,:,n_sl)+convex_edge;
                
                concaveBoundary=lobe_contours(concave_list,:);
                concaveboundary_indices = sub2ind([nRows,nCols],concaveBoundary(:,2),concaveBoundary(:,1));
                concave_edge = zeros(nRows,nCols);
                concave_edge(concaveboundary_indices)=1;
                concave_boundary(:,:,n_sl)=concave_boundary(:,:,n_sl)+concave_edge;
            end
        end
    end
end