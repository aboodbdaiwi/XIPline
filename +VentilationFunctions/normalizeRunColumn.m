function T = normalizeRunColumn(T)

runMode = T.RunMode;

for i = 1:numel(runMode)
    x = runMode{i};   % because your downstream uses cell-of-numeric

    if isempty(x) || (isnumeric(x) && isscalar(x) && isnan(x))
        continue
    end

    if ~isnumeric(x) || ~isscalar(x)
        error("Unknown RunMode type at row %d.", i)
    end

    switch x
        case 0
            runMode{i} = 2;
        case 1
            runMode{i} = 0;
        case 2
            runMode{i} = 1;
        case 3
            runMode{i} = 3;
        otherwise
            error("Unknown RunMode value at row %d: %g", i, x)
    end
end

T.RunMode = runMode;

end