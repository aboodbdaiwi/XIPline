function [first_max,first_max_value] = find_local_max(curve)
nLength = length(curve);
ref_value = curve(1);
for np = 2: nLength
    this_value = curve(np);
    if this_value > ref_value
        ref_value = this_value;
    else
        first_max = np-1;
        first_max_value = ref_value;
        break;
    end
end