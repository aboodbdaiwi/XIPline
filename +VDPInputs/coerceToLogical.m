function x = coerceToLogical(x)
%COERCETOLOGICAL Convert spreadsheet-style logical representations to logical.

    if islogical(x)
        return;
    end

    if isnumeric(x)
        x = x ~= 0;
        return;
    end

    x = string(x);
    x = upper(strtrim(x));

    trueMask  = ismember(x, ["TRUE","T","YES","Y","1"]);
    falseMask = ismember(x, ["FALSE","F","NO","N","0",""]);

    out = false(size(x));
    out(trueMask) = true;
    out(falseMask) = false;

    unknownMask = ~(trueMask | falseMask);
    if any(unknownMask)
        warning("coerceToLogical:UnknownValues", ...
            "Some Selected values were not recognized; setting them false.");
    end

    x = out;
end