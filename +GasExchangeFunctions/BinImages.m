function [BinnedImage] = BinImages(Image, Thresholds)
%Allocate binned image
BinnedImage = zeros(size(Image));
%Bin according to thresholds
for ii = 1:size(Image,3)
    for jj = 1:size(Image,2)
        for kk = 1:size(Image,1)
            value = Image(kk,jj,ii); %get value
            bin = find((value>Thresholds),1,'last')+1;%determine bin
            if isempty(bin)
                bin = 1;
            end
            BinnedImage(kk,jj,ii) = bin;%set bin
        end
    end
end

end