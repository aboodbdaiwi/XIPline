function [temp_handles,lobe_label_values,lobe_label_list,valid_lobe_labels] = identify_lobes(temp_handles)
%
% Pair label label values with lobe names;
% if lobemask does not exist, it will try to separate left and right lungs
% from lungsmask with lobe values 5 and 2 respectively. However, it
% currently can NOT separate the left and right lungs when connected.
%
%                  temp_handles.lungmask (binary)
%                              .lobemask               values
%                              .lobemaskRUL                 2
%                              .lobemaskRML                 3
%                              .lobemaskRLL                 4
%                              .lobemaskLUL                 5
%                              .lobemaskLLL                 6
% W. Zha @2014
%%

if ~isfield(temp_handles,'lobemask')
    temp_handles.lobemask = temp_handles.lungsmask;
elseif ~isfield(temp_handles,'lungsmask')
    error('detect_defects_by_kmeans_sarp:identify_lobes','lungsmask does not existed. \n');
end
[nRows,nCols,nSlices]=size(temp_handles.lungsmask);
contoured_slices = find(squeeze(sum(sum(temp_handles.lungsmask,1),2))>0);
if isempty(contoured_slices)
   lobe_label_values = [];
   lobe_label_list = [];
   valid_lobe_labels = [];
   return;
end
label_list = double(unique(temp_handles.lobemask(:)));
valid_lobe_labels = label_list(label_list>1 & label_list<7);
if isempty(valid_lobe_labels) || length(valid_lobe_labels)<5
    
    lobe_label_list = {'rul','lul'};
    lobe_label_values =[2,5];
    valid_lobe_labels = lobe_label_values;
        if isfield(temp_handles,'lobemaskRUL')
            temp_handles.lobemask=temp_handles.lobemaskRUL*lobe_label_values(1)+temp_handles.lobemaskLUL*lobe_label_values(2);
        else
            [nRows,nCols,nSlices]=size(temp_handles.lobemask);
            %% use mid-slice lung to label left and right lungs
            n_sl = round(nSlices/2);
            mid_lung = temp_handles.lungsmask(:,:,n_sl);
            
            [labeled_mask,nLungPieces] = bwlabel(mid_lung);
            if nLungPieces > 0
                pixels = zeros(nLungPieces,1);
                
                for np = 1:nLungPieces
                    pixels(np) = sum(labeled_mask(:) == np);
                end
                [sorted_pixels, sorted_ind] = sort(pixels, 'descend');
                
                if nLungPieces>1 && (sorted_pixels(1) - sorted_pixels(2) < 1000 || (sorted_pixels(2)>1000)) % two separated lungs

                     mid_lung(labeled_mask ~= sorted_ind(1) & labeled_mask ~=sorted_ind(2)) = 0;
                     [labeled_mask,nLungPieces] = bwlabel(mid_lung);
                     
                else % connected lungs
                    mid_lung(labeled_mask ~= sorted_ind(1)) = 0;
                    separated_lung = VentilationFunctions.Akmeans.separate_connected_lung(mid_lung);
                    [labeled_mask,nLungPieces] = bwlabel(separated_lung);

                    
                end
            else
                   lobe_label_values = [];
                   lobe_label_list = [];
                   valid_lobe_labels = [];
                   fprintf(' Unable to identify lobes in the mid lung slice #%d.\n',n_sl);
                   return;
            end
            %% Calculate the centers of the left and right lungs
            lobe_centers = zeros(nLungPieces, 2);
            RUL_center = [NaN, NaN];  % default if not assigned
            LUL_center = [NaN, NaN];
            
            for n_ob = 1:nLungPieces
                this_lobe = zeros(nRows, nCols);
                this_lobe(labeled_mask(:) == n_ob) = 1;
            
                lobe_contours = VentilationFunctions.Akmeans.mask2poly(this_lobe, 'Exact', 'CW');
                invalids = (lobe_contours(:,1) < 0 | lobe_contours(:,2) < 0);
                lobe_contours(invalids,:) = [];
            
                if ~isempty(lobe_contours)
                    lobe_centers(n_ob,:) = mean(lobe_contours);
                else
                    lobe_centers(n_ob,:) = [NaN, NaN];  % fallback if no valid contour
                end
            
                % Assign RUL and LUL if there are at least two lung pieces
                if nLungPieces >= 2 && n_ob == 2
                    if lobe_centers(1,1) < lobe_centers(2,1)
                        RUL_center = lobe_centers(1,:);
                        LUL_center = lobe_centers(2,:);
                    else
                        RUL_center = lobe_centers(2,:);
                        LUL_center = lobe_centers(1,:);
                    end
                end
            end

            % Initialize masks
            temp_handles.lobemaskRUL = zeros(nRows, nCols, nSlices);
            temp_handles.lobemaskLUL = zeros(nRows, nCols, nSlices);
            
            for sl = contoured_slices.'
                slice_lobemask = temp_handles.lungsmask(:,:,sl);
                if sum(slice_lobemask(:)) > 0
                    slice_lobemask = bwareaopen(slice_lobemask, 30);  % Remove small areas
                    [labeled_mask, nObj] = bwlabel(slice_lobemask);
            
                    % Attempt to separate connected lungs in mid-slices
                    if (sl > contoured_slices(1) + 1 && sl < contoured_slices(end) - 2) && nObj < 2
                        separated_lung = VentilationFunctions.Akmeans.separate_connected_lung(slice_lobemask);
                        [labeled_mask, nObj] = bwlabel(separated_lung);
                    end
            
                    nob_pixels = zeros(nObj, 1);
                    for n_ob = 1:nObj
                        nob_pixels(n_ob) = sum(labeled_mask(:) == n_ob);
                    end
            
                    [~, sorted_ob_index] = sort(nob_pixels, 'descend');
            
                    % Initialize for current slice
                    lobe_centers = NaN(nObj, 2);
            
                    for n_ob = 1:nObj
                        this_lobe = zeros(nRows, nCols);
                        this_lobe(labeled_mask(:) == sorted_ob_index(n_ob)) = 1;
            
                        % Extract contours and compute centroid
                        lobe_contours = VentilationFunctions.Akmeans.mask2poly(this_lobe, 'Exact', 'CW');
                        invalids = (lobe_contours(:,1) < 0 | lobe_contours(:,2) < 0);
                        lobe_contours(invalids, :) = [];
            
                        if ~isempty(lobe_contours)
                            lobe_centers(n_ob, :) = mean(lobe_contours);
            
                            % Decide assignment using known RUL and LUL centers
                            if all(~isnan(RUL_center)) && all(~isnan(LUL_center))
                                if abs(lobe_centers(n_ob, 1) - RUL_center(1)) < abs(lobe_centers(n_ob, 1) - LUL_center(1))
                                    temp_handles.lobemaskRUL(:,:,sl) = temp_handles.lobemaskRUL(:,:,sl) + this_lobe;
                                else
                                    temp_handles.lobemaskLUL(:,:,sl) = temp_handles.lobemaskLUL(:,:,sl) + this_lobe;
                                end
                            else
                                % If only one known center or fallback case, assign to RUL
                                temp_handles.lobemaskRUL(:,:,sl) = temp_handles.lobemaskRUL(:,:,sl) + this_lobe;
                            end
                        end
                    end
                end
            end
            
            % Combine masks
            temp_handles.lobemask = temp_handles.lobemaskRUL * 2 + temp_handles.lobemaskLUL * 5;

        end
        
else
        lobe_label_values = 2:6;
    lobe_label_list = {'rul','rml','rll','lul','lll'};
    for n_lobe = 1:numel(lobe_label_list)
        lobename = sprintf('lobemask%s',upper(lobe_label_list{n_lobe}));
        if ~isfield(temp_handles,lobename)
            temp_handles.(lobename) = zeros(nRows,nCols,nSlices);
            temp_handles.(lobename)(temp_handles.lobemask == lobe_label_values(n_lobe)) = 1;
            
        end
    end

end


