function clean_mask = remove_small_pieces(ori_mask,slice_range)
if nargin<2
    slice_range = 1:size(ori_mask,3);
end
if size(slice_range,1)~=1
    slice_range = slice_range.';
end
% [L,num] = bwlabeln(closed_reserved_mask);
clean_mask = ori_mask;
for nsl = slice_range
    slicemask = ori_mask(:,:,nsl);
    [labeled_slice,num] = bwlabel(slicemask);
    if num>2
        pixels = zeros(num,1);
        for nob = 1:num
            pixels(nob) = sum(labeled_slice(:)==nob);
        end
        [~, sorted_ind ] = sort(pixels,'descend');
        slicemask(labeled_slice~=sorted_ind(1) & labeled_slice~=sorted_ind(2)) = 0;
    end
    clean_mask(:,:,nsl) = slicemask;
end