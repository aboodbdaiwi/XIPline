function [combined_mask,combined_lobemask, gasmask,trachea,bodymask,seeds] = segment_proton_gas(proton,gas,seeds,RGthresh,acquisition)
% Use a seeded region-growing segmentation on proton and then resue these
% seeds on gas for a coarse gas mask to pick up apparent gas signals
% outside the proton lung mask. Trachea was estimated using seeded
% region-growing on gas image and excluded from both the proton and 
% gas mask. All seeds are automated generated. However, the manual
% selection of a trachea seed may be required.
%
%
% W. Zha@ Spring, 2016
%%
if nargin < 3 || isempty(seeds)
    seeds.right = [];
    seeds.left = [];
    seeds.trachea = [];
end
if nargin<2
    gas = zeros(size(proton));
    
end
if nargin<1
    error('proton_segmentation:segment_proton_gas','Proton data was not detected.\n');
end

 gasmask = zeros(size(gas));
    trachea = zeros(size(proton));
    
    proton = double(proton);
    proton = 255*(proton-min(proton(:)))/(max(proton(:))-min(proton(:)));
reg_tran_proton = hist_transformation(proton,385);
nSlices = size(reg_tran_proton,3);
%% BODY MASK from K-means
bodymask = kmeans_clustering(reg_tran_proton,4);
bodymask(bodymask<=1)=0;
bodymask(bodymask>0)=1;

se2 = strel('disk',2);
for nsl = 1:nSlices
    bodymask(:,:,nsl) = imerode(imfill(bodymask(:,:,nsl),'holes'),se2);
end

%% leakage
for nsl = 2:(nSlices-1)
    slice_body = bodymask(:,:,nsl);
    mutual = bodymask(:,:,nsl-1) & bodymask(:,:,nsl+1);
    diff = mutual - slice_body;
    diff(diff<0) = 0;
    [labeled_diff,nDiff] = bwlabel(diff);
    if nDiff>0
        diffPixels = zeros(nDiff,1);
        for np = 1: nDiff
            diffPixels(np) = sum(labeled_diff(:) == np);
        end
        [maxPixel,max_diff] = max(diffPixels);
    if maxPixel>1000
        slice_body(labeled_diff ==max_diff) = 1;
        bodymask(:,:,nsl) = imfill(slice_body,'holes');
    end
    end
end

proton = 255*(proton-min(proton(:)))/(max(proton(:))-min(proton(:)));

%% coarse proton lung mask from seeded region growing
[RGlungmask,boundingbox,seeds] = seeded_region_growing_on_lung(proton, bodymask, seeds,RGthresh);

lungVolList = squeeze(sum(sum(RGlungmask,1),2));
[~,maxSL] = max(lungVolList);


diaphragm_sl = locate_diaphragm_slice(RGlungmask,maxSL);
[~,~,closed_reserved_mask] = parse_mask(RGlungmask);
 figure,data_mask_overlay(gas,RGlungmask);
% figure,data_mask_overlay(gas,reserved_mask);
%% smooth mask
se1 = strel('disk',1);
se2 = strel('disk',2);

for nsl = 1:nSlices
    closed_reserved_mask(:,:,nsl) = imdilate(imerode(imclose(closed_reserved_mask(:,:,nsl),se2),se1),se1);
end
% figure,data_mask_overlay(gas,closed_reserved_mask);
%% remove small pieces
closed_reserved_mask = remove_small_pieces(closed_reserved_mask,1:maxSL);
% figure,data_mask_overlay(gas,closed_reserved_mask);
% figure,data_mask_overlay(reg_tran_proton,closed_reserved_mask);

%% KMEANS CLUSTERING on masked proton
nClusters = 4;
clst_proton = kmeans_clustering(reg_tran_proton.*closed_reserved_mask,nClusters);
mask_tmp = closed_reserved_mask;
mask_tmp(clst_proton == nClusters)=0;
% figure,data_mask_overlay(reg_tran_proton,mask_tmp);

%%
mask_tmp2 = remove_small_pieces(mask_tmp,1:diaphragm_sl);
clean_mask = mask_tmp-mask_tmp2;
for nsl = 1:nSlices
    clean_mask(:,:,nsl) = imdilate(clean_mask(:,:,nsl),se1);
end

corrected_mask = closed_reserved_mask - clean_mask;
corrected_mask(corrected_mask<0) = 0;
% figure,data_mask_overlay(gas,corrected_mask);

% figure,imshow3Dfull([clean_mask,removed_region])
%% CHECK NEIGHBOR SLICES
for nsl = 2: (nSlices-1)
    combined_neighbor = corrected_mask(:,:,nsl-1) | corrected_mask(:,:,nsl+1);
    mutual_neighbor =  corrected_mask(:,:,nsl-1) & corrected_mask(:,:,nsl+1);
    unique_region = corrected_mask(:,:,nsl) - combined_neighbor;
    unique_region(unique_region<0) = 0;
    corrected_mask(:,:,nsl) = corrected_mask(:,:,nsl) - unique_region;
    
    missing_region = mutual_neighbor - corrected_mask(:,:,nsl);
    missing_region(missing_region<0) = 0;
    if nsl > 3 && nsl<diaphragm_sl-1
        corrected_mask(:,:,nsl) = corrected_mask(:,:,nsl) + missing_region;
    end
end
for ns = 1:maxSL
    corrected_mask(:,:,ns) = imfill(corrected_mask(:,:,ns),'holes');
end
if sum(corrected_mask(:))==0
    combined_lobemask = [];
    combined_mask = [];
    return;
end

%% DETECT BETWEE LUNG BINARY MASK FROM THE CORRECTED PROTON LUNG MASK
switch acquisition
    case 'Coronal'
        
        [~,~,between_lungs] = split_lung_boundaries_coronal(corrected_mask,1:nSlices);
    case 'Axial'

        between_lungs = ones(size(proton));
end

%% REGION GROWING ON GAS -optional
if sum(gas(:))>0 && sum(gas(:))~=Inf
    gasmaskR = regiongrowing(gas/max(gas(:)),0.08,seeds.right);
    gasmaskL = regiongrowing(gas/max(gas(:)),0.08,seeds.left);
    gasmask = gasmaskR | gasmaskL;
    %% figure,imshow3Dfull(hemask)
    for nsl = 1:nSlices
        gasmask(:,:,nsl) = imfill(imclose(gasmask(:,:,nsl),se2),'holes');
        gasmask(:,:,nsl) = bwareaopen(gasmask(:,:,nsl),30);
    end
    % figure,data_mask_overlay(gas,hemask)
    
    [trachea, seeds] = estimate_trachea(gas, bodymask, seeds,acquisition);
    % if sum(trachea(:))>5000
    %     trachea = zeros(nRows,nCols,nCols);
    % end
    
    trachea = trachea.*between_lungs;
    gas_slices = find(squeeze(sum(sum(gasmask,1),2))>0);
    trachea_sl = gas_slices(1) + floor(length(gas_slices)/4):floor(3*length(gas_slices)/4);
    non_trachea = setdiff(1:nSlices,trachea_sl);
    between_lungs(:,:,non_trachea) = 0;
    %% REMOVE DETECTED TRACHEA AND MERGE GAS SIGNAL INTO PROTON LUNG MASK
    
    corrected_mask(trachea == 1) = 0;
    gasmask = gasmask.*boundingbox;
    gasmask(trachea == 1) = 0;
    gasmask = gasmask.*~between_lungs();
    % c = combine_proton_he_mask(corrected_mask,hemask);
    if sum(gasmask(:))> sum(corrected_mask(:))+1e4
        combined_mask = corrected_mask;
    else
        %     switch acquisition
        %         case 'Coronal'
        combined_mask = merge_he_into_proton(corrected_mask,gasmask);
        %         case 'Axial'
        %             combined_mask = combine_proton_he_mask(corrected_mask,gasmask);
        %             combined_mask = merge_he_into_proton(combined_mask,gasmask);
        %     end
        
    end
else
%     gasmask = zeros(size(gas));
    combined_mask = corrected_mask;
%     trachea = zeros(size(proton));
end
% combined_lobemask = merge_ctlobemask_into_lungmask(lobemask,combined_mask);
combined_mask = morphometric_closing(combined_mask);
combined_lobemask = [];
%% -----------------------------------------------------------------------------
%% SEEDED REGION GROWING

function [protonmask, boundingbox,seeds,lobemask] = seeded_region_growing_on_lung(proton,bodymask, seeds,RGthresh)
[nRows, nCols, nSlices] = size(proton);
% proton_uint8 = uint8(proton);
% cutoff_intn=locate_histogram_landmark_intensity(proton,.999);
% if cutoff_intn<100
% proton_uint8(proton_uint8>cutoff_intn) = 0;%cutoff_intn;
% proton = double(proton_uint8);
% proton = 255*(proton-min(proton(:)))/(max(proton(:))-min(proton(:)));
% else
reg_tran_proton = hist_transformation(proton,385);
reg_tran_proton = reg_tran_proton.*bodymask;


% end
% [~,cropmask]=get_plotting_range(bodymask,-30);
% bodymask = cropmask.*bodymask;
if isempty(seeds.right)
[seedsLung, midLungmask, ~, boundingbox] = generate_seed_point(proton,bodymask);

if isempty(seeds.right)
    seeds.right = seedsLung.right;
end

if isempty(seeds.left)
    seeds.left = seedsLung.left;
end
else
    plotting_range = get_plotting_range(bodymask,2);
boundingbox = zeros(nRows, nCols, nSlices);
boundingbox(plotting_range.y,plotting_range.x,:) = 1;
midLungmask=zeros(nRows,nCols);
end
%figure,data_mask_overlay(reg_tran_proton,bodymask)
if ~isempty(seeds.right)

    [maskRight,seeds.right] = regiongrowing(reg_tran_proton/max(reg_tran_proton(:)),RGthresh,seeds.right);
else
        [maskRight,seeds.right] = regiongrowing(reg_tran_proton/max(reg_tran_proton(:)),RGthresh);
end

if ~isempty(seeds.left)
    [maskLeft,seeds.left] = regiongrowing(reg_tran_proton/max(reg_tran_proton(:)),RGthresh,seeds.left);
else
        [maskLeft,seeds.left] = regiongrowing(reg_tran_proton/max(reg_tran_proton(:)),RGthresh);
end

RGlungmask = (maskRight | maskLeft).*boundingbox;
bbox_boundary = zeros(nRows,nCols,nSlices);
se1 = strel('disk',1);
for nsl = 1:nSlices
    RGlungmask(:,:,nsl) = imfill(RGlungmask(:,:,nsl),'holes');
    bbox_boundary(:,:,nsl) = boundingbox(:,:,nsl) - imerode(boundingbox(:,:,nsl),se1);
end

lung_slices = find(squeeze(sum(sum(RGlungmask,1),2))>0);
midSlice = round(nSlices/2);
reordered_slices = [(midSlice+1): lung_slices(end), (midSlice-1):-1:lung_slices(1)];
lobevalues = unique(midLungmask);
lobevalues(lobevalues == 0) = [];
%% Use right- and left-lung labeled midLungmask to refine RGlungmask
if ~isempty(lobevalues)
    %%
    
lobemask = zeros(nRows,nCols,nSlices);
lobemask(:,:,midSlice) = midLungmask;
    for nlobe = 1:length(lobevalues)
        template_lobe = zeros(nRows, nCols);
        template_lobe(midLungmask == lobevalues(nlobe)) = 1;
        for nsl = reordered_slices
            slice_lung = RGlungmask(:,:,nsl).*boundingbox(:,:,nsl);
            slice_lung(slice_lung<0)=0;
            
            [labeled_lung, nLungPieces] = bwlabel(slice_lung);
            
            for np = 1: nLungPieces
                piecemask = zeros(nRows,nCols);
                piecemask(labeled_lung == np) = 1;
                overlap_lung = piecemask.*template_lobe;
                bound = bbox_boundary(:,:,nsl);
                overlap_bbox = piecemask.*bound;
                if sum(overlap_lung(:)>0) && sum(overlap_bbox(:))<sum(bound(:))/8
                    lobemask(:,:,nsl) =  lobemask(:,:,nsl) + piecemask.*lobevalues(nlobe);
                end
            end
            
            if nsl>midSlice && nsl<max(reordered_slices)
            template_lobe = lobemask(:,:,nsl-1);
            elseif nsl==max(reordered_slices)
                        template_lobe = zeros(nRows, nCols);
        template_lobe(midLungmask == lobevalues(nlobe)) = 1;

            else
                template_lobe = lobemask(:,:,nsl+1);
            end
        end
    end
    %%
 %   figure,imshow3Dfull(lobemask);
    protonmask = lobemask;
protonmask(protonmask>0) = 1;
else
    protonmask = RGlungmask;
    lobemask = [];
end
