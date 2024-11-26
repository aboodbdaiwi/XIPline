function [convex_boundary,concave_boundary,between_lungs] = split_lung_boundaries_coronal(lungsmask,slice_range)
% Input: -lungsmask: binary lung masks
%        -slice_range: the list of slices for boundary split.
%
% Output:-convex_boundary: 3D matrix that contains the edge masks of "outter
%                          boundary" of the lungs
%        -concave_boundary: 3D matrix that contains the edge masks of "inner
%                          boundary" of the lungs where the vessels may be
%                          attached to.
%
% W. Zha @ 3/30/2015
%%
if nargin<2
    slice_range=1:size(lungsmask,3);
end
[nRows,nCols,nSlices]=size(lungsmask);
convex_boundary=zeros(nRows,nCols,nSlices);
concave_boundary=zeros(nRows,nCols,nSlices);
%% Save the binary right and left lung masks into lobemask.RUL/LUL
% The possible gaps among RUL, RML and RLL labeled in SARP data may
% cause fake lung boundaries, so RUL and LUL are automatically
% relabeled from lung masks here regardless.
temp_handles.lungsmask = lungsmask;
temp_handles.lobemask = lungsmask;
[temp_handles,lobe_label_values] = identify_lobes(temp_handles);

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
    lobemask.RUL=temp_handles.lobemaskRUL.*lungsmask;
    lobemask.LUL=temp_handles.lobemaskLUL.*lungsmask;
end

se2=strel('disk',2);
%% separate the convex and concave portion of the lung boundaries.
lobe_list={'RUL','LUL'};
between_lungs = zeros(nRows, nCols, nSlices);
for n_ind=1:length(slice_range)
    n_sl=slice_range(n_ind);
            all_concave = [];
    for n_lobe=1:numel(lobe_list)
        lobe_name=lobe_list{n_lobe};
        this_lobe=lobemask.(lobe_name);

        slice_mask=this_lobe(:,:,n_sl);
        slice_mask=imerode(imdilate(imfill(slice_mask,'holes'),se2),se2);
        [labeled_slice,num]=bwlabel(slice_mask);
        

        if num>1
            slicepixels = zeros(num,1);
                    for np=1:num
                        slicepixels(np)=sum(labeled_slice(:)==np);
                    end
        [~,piece_index] = max(slicepixels);
        slice_mask(labeled_slice~=piece_index)=0;
        end
        if sum(slice_mask(:))>0
            % the right/left lung masks may be splited toward the base due to the appearance of diaphragm.
            
            lobe_contours = mask2poly(slice_mask,'Exact','CW');
            invalids = (lobe_contours(:,1)<=0 | lobe_contours(:,2)<=0 | lobe_contours(:,1)>nCols | lobe_contours(:,2)>nRows);
            lobe_contours(invalids,:)=[];
%             nPoints=size(lobe_contours,1);
            this_lobe_center=mean(lobe_contours);
            
            switch lobe_name
                case 'RUL'
                    %%
                    quarter4_ind=find(lobe_contours(:,1)>this_lobe_center(1) & lobe_contours(:,2)>this_lobe_center(2));
                    quarter1_ind=find(lobe_contours(:,1)>this_lobe_center(1) & lobe_contours(:,2)<this_lobe_center(2));
                    %%
                case 'LUL'
                    %%
                    quarter4_ind=find(lobe_contours(:,1)<this_lobe_center(1) & lobe_contours(:,2)>this_lobe_center(2));
                    quarter1_ind=find(lobe_contours(:,1)<this_lobe_center(1) & lobe_contours(:,2)<this_lobe_center(2));
                    if length(quarter4_ind)<length(quarter1_ind)/4
                        end_ind=quarter1_ind(1)-length(quarter1_ind);
                        quarter4_ind=(quarter1_ind(1)-1):-1:max(end_ind,20);
                    end
                    %%
            end
            dist_list=(lobe_contours(quarter1_ind,1)-this_lobe_center(1)).^2+...
                (lobe_contours(quarter1_ind,2)-this_lobe_center(2)).^2;
            [~,max_ind1]=max(dist_list) ;
            dist_list2=(lobe_contours(quarter4_ind,1)-this_lobe_center(1)).^2+...
                (lobe_contours(quarter4_ind,2)-this_lobe_center(2)).^2;
            [~,max_ind2]=max(dist_list2);
            apex_ind=quarter1_ind(max_ind1);
            base_ind=quarter4_ind(max_ind2);
            %%
            %             figure,imshow(slice_mask,[])
            %            hold on; plot(lobe_contours(:,1),lobe_contours(:,2));
            %             hold on;plot(this_lobe_center(1),this_lobe_center(2),'m*')
            %             hold on,plot(lobe_contours(quarter1_ind,1),lobe_contours(quarter1_ind,2),'g.-');
            %             hold on,plot(lobe_contours(apex_ind,1),lobe_contours(apex_ind,2),'go');
            %             hold on,plot(lobe_contours(quarter4_ind,1),lobe_contours(quarter4_ind,2),'r.-');
            %             hold on,plot(lobe_contours(base_ind,1),lobe_contours(base_ind,2),'ro');
            %%
            concave_list= base_ind:apex_ind;
            convex_list=setdiff(1:size(lobe_contours,1),concave_list);
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
            %                 figure,imagesc(convex_edge*10+concave_edge)
            %%
            if n_lobe ==1
            all_concave = [all_concave;concaveBoundary];
            else
                all_concave = [all_concave;concaveBoundary(end:-1:1,:)];
            end
        end
    end
    if sum(concave_boundary(:))>0 && ~isempty(all_concave)
        between_lungs(:,:,n_sl) = poly2mask(all_concave(:,1),all_concave(:,2),nRows,nCols);
    end
end