function out = normalizeScanDate(x)
%NORMALIZESCANDATE Convert mixed date representations to yyyy-MM-dd string.
%
% Handles inputs like:
%   "3/5/2026"
%   "20240626"
%   numeric Excel-ish or imported numeric-looking values as strings

    x = string(x);
    x = strtrim(x);

    out = strings(size(x));

    for i = 1:numel(x)
        xi = x(i);

        if ismissing(xi) || xi == ""
            out(i) = "";
            continue;
        end

        % Case 1: yyyyMMdd
        if strlength(xi) == 8 && all(isstrprop(char(xi), 'digit'))
            try
                dt = datetime(xi, "InputFormat", "yyyyMMdd");
                out(i) = string(dt, "yyyyMMdd");
                continue;
            catch
            end
        end

        % Case 2: M/d/yyyy
        try
            dt = datetime(xi, "InputFormat", "M/d/yyyy");
            out(i) = string(dt, "yyyyMMdd");
            continue;
        catch
        end

        % Case 3: already parseable by datetime
        try
            dt = datetime(xi);
            out(i) = string(dt, "yyyyMMdd");
            continue;
        catch
        end

        % Fallback: preserve original
        out(i) = xi;
    end
end