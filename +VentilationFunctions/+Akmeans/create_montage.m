function montage_img = create_montage(data,mask,SlicePerRow)
% save multi slices into a single matrix for display using imshow so that
% the window/level is the same for the entire 3D data set.
% W. Zha @ Fall 2016
if nargin<3
    SlicePerRow = 7;
end
range = get_plotting_range(mask);

snap.nRows = length(range.y);
snap.nCols = length(range.x);
snap.nSlices = length(range.z);

numRows = ceil(length(range.z)/SlicePerRow);
montage_img = zeros(snap.nRows*numRows, snap.nCols*SlicePerRow);

for nsl = 1:length(range.z)
    row_index = ceil(nsl/SlicePerRow);
    col_index = rem(nsl,SlicePerRow)+(rem(nsl,SlicePerRow)==0)*SlicePerRow;
    montage_img((1:snap.nRows)+snap.nRows*(row_index-1), (1:snap.nCols)+snap.nCols*(col_index-1)) = data(range.y,range.x,range.z(nsl));
end