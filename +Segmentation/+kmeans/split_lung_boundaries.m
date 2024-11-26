function [convex_boundary,concave_boundary,...
    region_inbetween_lungs,mask_spline] = split_lung_boundaries(mask_struct,slice_range,binary_flag)
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
lungsmask = mask_struct.lungsmask;
if nargin<2
    slice_range = 1:size(lungsmask,3);
    binary_flag = true;
elseif nargin < 3
    binary_flag = true;
    
end


[nRows,nCols,nSlices] = size(lungsmask);
convex_boundary = zeros(nRows,nCols,nSlices);
concave_boundary = zeros(nRows,nCols,nSlices);

[mask_struct,lobe_label_values] = identify_lobes(mask_struct);
mask_spline = zeros(nRows,nCols,nSlices);

%% Save the binary right and left lung masks into lobemask.RUL/LUL
% The possible gaps among RUL, RML and RLL labeled in SARP data may
% cause fake lung boundaries, so RUL and LUL are automatically
% relabeled from lung masks here regardless.

if length(lobe_label_values) == 5
    RLL = zeros(nRows,nCols,nSlices);
    RML = zeros(nRows,nCols,nSlices);
    LLL = zeros(nRows,nCols,nSlices);
    se = strel('disk',2);
    for nsl = 1:nSlices
        RLL(:,:,nsl) = imdilate(mask_struct.lobemaskRLL(:,:,nsl),se);
        RML(:,:,nsl) = imdilate(mask_struct.lobemaskRML(:,:,nsl),se);
        LLL(:,:,nsl) = imdilate(mask_struct.lobemaskLLL(:,:,nsl),se);
    end
    lobemask.RUL = (mask_struct.lobemaskRUL+RML+RLL).*lungsmask;
    lobemask.RUL(lobemask.RUL>0) = 1;
    lobemask.LUL = (mask_struct.lobemaskLUL+LLL).*lungsmask;
    lobemask.LUL(lobemask.LUL>0) = 1;
elseif length(lobe_label_values)  == 2
    lobemask.RUL = mask_struct.lobemaskRUL;
    lobemask.LUL = mask_struct.lobemaskLUL;
end
%% separate the convex and concave portion of the lung boundaries.
lobe_list = {'RUL','LUL'};
lobe_values = [2, 5];
region_inbetween_lungs = zeros(nRows,nCols,nSlices);
% [~,maxLungSL] = max(squeeze(sum(sum(lungsmask,1),2)));

% slice_reordered = [maxLungSL:slice_range(end), (maxLungSL-1):-1:slice_range(1)];
slice_reordered = slice_range;
for n_ind = 1:length(slice_range)
    n_sl = slice_reordered(n_ind);
    slname = sprintf('sl%d',n_sl);
    for n_lobe = 1:numel(lobe_list)
        lobe_name = lobe_list{n_lobe};
        this_lobe = lobemask.(lobe_name);
        slice_mask = this_lobe(:,:,n_sl);
        if sum(slice_mask(:))>0
                [labeled_mask,nObj] = bwlabel(slice_mask);
                for n = 1:nObj
                    object_mask = zeros(nRows,nCols);
                    object_mask(labeled_mask == n) = 1;
                    lobe_contours = mask2poly(object_mask,'Exact','CW');
%                     invalids = (lobe_contours(:,1)<0 | lobe_contours(:,2)<0 | lobe_contours(:,1)>nRows | lobe_contours(:,2)>nCols);
                    invalids = (lobe_contours(:,1)<0 | lobe_contours(:,2)<0 | lobe_contours(:,1)>nCols | lobe_contours(:,2)>nRows);
                    lobe_contours(invalids,:) = [];
                    
                    
                    [convex_list,concave_list] = find_conv_region(lobe_contours);
                    
                    convexBoundary = lobe_contours(convex_list,:);
                    
                    
                    % fit b_spline
                    %                 nPoints = size(convexBoundary,1);
                    %                     pp1 = splinefit(1:nPoints,[convexBoundary(:,1).';convexBoundary(:,2).'],8,spline_order,'r');
                    %                  fitted_convex = ppval(pp1,1:nPoints);
                    %                 convexBoundary = round(fitted_convex.');
                    
                    convexboundary_indices = sub2ind([nRows,nCols],convexBoundary(:,2),convexBoundary(:,1));
                    %                 if n_sl < slice_range(1)+2 % check consistency towards apex
                    %                     %%
                    %                     prev_convex_mask = imdilate(convex_boundary(:,:,n_sl+1),se1);
                    %
                    %                      on_prev_convex = prev_convex_mask(convexboundary_indices);
                    %                     nOverlaps = zeros(1000,1);
                    %
                    %                     nOverlaps(2) = sum(on_prev_convex);
                    %                     nCount = 2;
                    %                     while nOverlaps(nCount) >= nOverlaps(nCount-1) && nOverlaps(nCount)< length(convexboundary_indices)
                    %                         nCount = nCount + 1;
                    %
                    %                          prev_convex_mask = imdilate(prev_convex_mask,se1);
                    %                          on_prev_convex = prev_convex_mask(convexboundary_indices);
                    %                          nOverlaps(nCount) = sum(on_prev_convex);
                    %                     end
                    %%
                    %                     figure,imshow(prev_convex_mask,[]);
                    %                     hold on;plot(lobe_contours(:,1),lobe_contours(:,2),'.')
                    %                     if sum(on_prev_convex)< .5*length(convexboundary_indices)
                    %                         all_contour_indices = sub2ind([nRows,nCols],lobe_contours(:,2),lobe_contours(:,1));
                    %                         on_prev_convex = prev_convex_mask(all_contour_indices);
                    %                         convex_indices = find(on_prev_convex == 1);
                    %                         convex_points = lobe_contours(convex_indices,:);
                    %                         [convex_list,concave_list] = find_conv_region(lobe_contours,convex_points,convex_indices);
                    %                         convexBoundary = lobe_contours(convex_list,:);
                    %                         convexboundary_indices = sub2ind([nRows,nCols],convexBoundary(:,2),convexBoundary(:,1));
                    %                     end
                end
                
       
            convex_edge = zeros(nRows,nCols);
            convex_edge(convexboundary_indices) = 1;
            
            concaveBoundary = lobe_contours(concave_list,:);
            
            % fit b_spline
            %                 nPoints = size(concaveBoundary,1);
            %                     pp1 = splinefit(1:nPoints,[concaveBoundary(:,1).';concaveBoundary(:,2).'],8,spline_order,'r');
            %                  fitted_concave = ppval(pp1,1:nPoints);
            %                 concaveBoundary  =  round(fitted_concave.');
            
            concaveboundary_indices  =  sub2ind([nRows,nCols],concaveBoundary(:,2),concaveBoundary(:,1));
            concave_pts.(lobe_name).(slname) = concaveBoundary;
            concave_edge = zeros(nRows,nCols);
            concave_edge(concaveboundary_indices) = 1;
            
            if binary_flag
                concave_boundary(:,:,n_sl) = concave_boundary(:,:,n_sl) + concave_edge;
                convex_boundary(:,:,n_sl) = convex_boundary(:,:,n_sl) + convex_edge;
            else
                concave_boundary(:,:,n_sl) = concave_boundary(:,:,n_sl) + concave_edge*lobe_values(n_lobe);
                convex_boundary(:,:,n_sl) = convex_boundary(:,:,n_sl) + convex_edge*lobe_values(n_lobe);
                
            end
            
            
            mask_spline(:,:,n_sl) = mask_spline(:,:,n_sl) + poly2mask([convexBoundary(:,1);concaveBoundary(:,1)],...
                [convexBoundary(:,2);concaveBoundary(:,2)],nRows,nCols);
        else
            concave_pts.RUL = [];
            concave_pts.LUL = [];
        end
    end
    if isfield(concave_pts.RUL,slname) && isfield(concave_pts.LUL,slname)
%         middle_region = poly2mask([concave_pts.RUL.(slname)(:,1);concave_pts.LUL.(slname)(end:-1:1,1)],...
%             [concave_pts.RUL.(slname)(:,2);concave_pts.LUL.(slname)(end:-1:1,2)],nRows,nCols);
                middle_region = poly2mask([concave_pts.RUL.(slname)(:,1);concave_pts.LUL.(slname)(end:-1:1,1)],...
            [concave_pts.RUL.(slname)(:,2);concave_pts.LUL.(slname)(end:-1:1,2)],nRows,nCols);
        region_inbetween_lungs(:,:,n_sl) = middle_region;
    end
 
end

%%
% function ctr3D = locate_3D_center(mask3D)
% masksum = squeeze(sum(mask3D,3));
% % figure,imagesc(lungsum)
% row_list = find(sum(masksum,1)>0);
% col_list= find(sum(masksum,2)>0);
% ctr3D = [mean(row_list), mean(col_list)];
