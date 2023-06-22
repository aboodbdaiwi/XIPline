function zc = zero_cross(data)
%Function to find approximately where data crosses zero

zc = [];
for i = 2:length(data)
    if data(i) >= 0 && data(i-1) < 0
        if abs(data(i)) > abs(data(i-1))
            zc = [zc i-1];
        else
            zc = [zc i];
        end
    elseif data(i) <= 0 && data(i-1) > 0
        if abs(data(i)) > abs(data(i-1))
            zc = [zc i-1];
        else
            zc = [zc i];
        end
    end
end