function T = addMissingVar(T, varName, defaultValue)
%ADDMISSINGVAR Add variable only if missing.

    if ~ismember(varName, T.Properties.VariableNames)
        T.(varName) = defaultValue;
    end
end