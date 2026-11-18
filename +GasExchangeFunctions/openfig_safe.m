function h = openfig_safe(figName)

    if isfile(figName)
        figPath = figName;
    else
        d = dir(fullfile(pwd,'**',figName));
        if isempty(d)
            error('Could not locate figure: %s',figName);
        end
        figPath = fullfile(d(1).folder,d(1).name);
    end

    h = openfig(figPath);
end