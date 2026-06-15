% Median filter applied only in the first two dimensions.
% kernelSize can be 2 or 3.
%
% Example:
%   newarray = medFilter(oldarray);    % default 3x3
%   newarray = medFilter(oldarray,2);  % 2x2
%   newarray = medFilter(oldarray,3);  % 3x3

function newarray = medFilter_integer(oldarray,kernelSize)

if nargin < 2 || isempty(kernelSize)
    kernelSize = 3;
end

if ~ismember(kernelSize,[2 3])
    error('kernelSize must be either 2 or 3.');
end

newarray = oldarray;

if ndims(oldarray) == 2

    if kernelSize == 3
        for i = 2:size(oldarray,1)-1
            for j = 2:size(oldarray,2)-1
                list = oldarray((i-1):(i+1),(j-1):(j+1));
                newarray(i,j) = round(median(list(:),'omitnan'));
            end
        end
    else
        for i = 1:size(oldarray,1)-1
            for j = 1:size(oldarray,2)-1
                list = oldarray(i:i+1,j:j+1);
                newarray(i,j) = round(median(list(:),'omitnan'));
            end
        end
    end

elseif ndims(oldarray) == 3

    if kernelSize == 3
        for k = 1:size(oldarray,3)
            for i = 2:size(oldarray,1)-1
                for j = 2:size(oldarray,2)-1
                    list = oldarray((i-1):(i+1),(j-1):(j+1),k);
                    newarray(i,j,k) = round(median(list(:),'omitnan'));
                end
            end
        end
    else
        for k = 1:size(oldarray,3)
            for i = 1:size(oldarray,1)-1
                for j = 1:size(oldarray,2)-1
                    list = oldarray(i:i+1,j:j+1,k);
                    newarray(i,j,k) = round(median(list(:),'omitnan'));
                end
            end
        end
    end

elseif ndims(oldarray) == 4

    if kernelSize == 3
        for b = 1:size(oldarray,4)
            for k = 1:size(oldarray,3)
                for i = 2:size(oldarray,1)-1
                    for j = 2:size(oldarray,2)-1
                        list = oldarray((i-1):(i+1),(j-1):(j+1),k,b);
                        newarray(i,j,k,b) = round(median(list(:),'omitnan'));
                    end
                end
            end
        end
    else
        for b = 1:size(oldarray,4)
            for k = 1:size(oldarray,3)
                for i = 1:size(oldarray,1)-1
                    for j = 1:size(oldarray,2)-1
                        list = oldarray(i:i+1,j:j+1,k,b);
                        newarray(i,j,k,b) = round(median(list(:),'omitnan'));
                    end
                end
            end
        end
    end

else
    error('Input array must be 2D, 3D, or 4D.');
end

newarray(isnan(newarray)) = 0;
newarray = round(newarray);

end