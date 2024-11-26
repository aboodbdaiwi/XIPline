function localScale = get_local_scale(lungsmask)
%% Reference: Tzeng et al, "The difference in ventilation heterogeneity between
% asthmatic and healthy subjects quantified using hyperpolarized 3He MRI."
% Journal of applied physiology 2009.
% W. Zha @ 08/30/2016
%%
contoured_slices =  find(squeeze(sum(sum(lungsmask,1),2)));
nSlices = length(contoured_slices);
[nRows, nCols, ~] = size(lungsmask);
globalLengthScale = zeros(nSlices,2);

mask_struct.lungsmask = lungsmask;
mask_struct.lobemask = lungsmask;
mask_struct = identify_lobes(mask_struct);

lobes = {'lobemaskRUL','lobemaskLUL'};
for sl = 1:nSlices
    for n_lobe = 1:numel(lobes)
        lobename = lobes{n_lobe};
    sliceLung = mask_struct.(lobename)(:,:,contoured_slices(sl));
    [labeled_lung,nPieces] = bwlabel(sliceLung);
     for np = 1: nPieces
         pieceMask = zeros(nRows, nCols);
         pieceMask(labeled_lung == np) = 1;
         lobe_contours = mask2poly(pieceMask,'Exact','CW');
         invalids = (lobe_contours(:,1)<0 | lobe_contours(:,2)<0 | lobe_contours(:,1)>nRows | lobe_contours(:,2)>nCols);
         lobe_contours(invalids,:) = [];
         nPoints = size(lobe_contours,1);
         boundary_indices = sub2ind([nRows,nCols],lobe_contours(:,2)-1,lobe_contours(:,1));
         maxWidthCal = zeros(nPoints,1);
         if nPoints>0
         parfor n = 1:nPoints
%              line_seg = [1:nCols; ones(1,nCols)*lobe_contours(n,2)];
                          line_img = zeros(nRows,nCols);
             line_img(ones(1,nCols)*lobe_contours(n,2)-1,:)=1;
             intersects = find(line_img+pieceMask ==2);
             if ~isempty(intersects)
             maxWidthCal(n) = max(abs(boundary_indices(n)-intersects+1))/nCols;
             end
         end
              globalLengthScale(sl,n_lobe) = max(maxWidthCal);
         end
%          row_indices = find(sum(pieceMask,1));
     end

    end
end
sorted_globalLengthScale = sort(globalLengthScale,2,'descend');
global_max_length = .5*max(sorted_globalLengthScale(:,1))+.5*max(sorted_globalLengthScale(:,2));
localScale = round(.1*global_max_length);
