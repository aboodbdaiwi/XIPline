function x = coerceToDouble(x)
%COERCETODOUBLE Convert mixed spreadsheet column to double.

    if isnumeric(x)
        x = double(x);
        return;
    end

    x = string(x);
    x = strtrim(x);
    x(x == "") = missing;
    x = double(str2double(x));
end