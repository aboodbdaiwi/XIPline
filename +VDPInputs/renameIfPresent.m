function T = renameIfPresent(T, oldName, newName)
%RENAMEIFPRESENT Rename table variable if old exists and new does not.

    vars = T.Properties.VariableNames;

    if ismember(oldName, vars) && ~ismember(newName, vars)
        T = renamevars(T, oldName, newName);
    end
end