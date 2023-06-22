
%median Filter (3x3) kernel, applied only first two dimensions, 5-23-2016, 
function newarray=medFilter(oldarray)


if ndims(oldarray)==2;
    newarray = zeros(size(oldarray,1),size(oldarray,2));
    for i=2:size(oldarray,1)-1;
        for j=2:size(oldarray,2)-1;
            list = oldarray((i-1):(i+1),(j-1):(j+1));
            newarray(i,j) = median(median(list));
        end
    end
    newarray(isnan(newarray)) = 0;
end
 

if ndims(oldarray)==3;
    newarray = zeros(size(oldarray,1),size(oldarray,2),size(oldarray,3));
    for k=1:size(oldarray,3);
        for i=2:size(oldarray,1)-1;
            for j=2:size(oldarray,2)-1;
                list = oldarray((i-1):(i+1),(j-1):(j+1),k);
                newarray(i,j,k) = median(median(list));
            end
        end
    end
    newarray(isnan(newarray)) = 0;
end
 

if ndims(oldarray)==4;
    newarray = zeros(size(oldarray,1),size(oldarray,2),size(oldarray,3),size(oldarray,4));
    for b=1:size(oldarray,4);
        for k=1:size(oldarray,3);
            for i=2:size(oldarray,1)-1;
                for j=2:size(oldarray,2)-1;
                    list = oldarray((i-1):(i+1),(j-1):(j+1),k,b);
                    newarray(i,j,k,b) = median(median(list));
                end
            end
        end
    end
    newarray(isnan(newarray)) = 0;
end
 


