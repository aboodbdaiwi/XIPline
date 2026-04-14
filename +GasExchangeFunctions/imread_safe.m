function h = imread_safe(figName)

    % Look for the figure anywhere under the current analysis folder
    d = dir(fullfile(pwd,'**',figName));

    if isempty(d)
        error('Could not locate %s anywhere under %s', figName, pwd);
    end

    % If multiple matches exist, take the first one
    figPath = fullfile(d(1).folder, d(1).name);

    fprintf('Reading image: %s\n', figPath);

    h = imread(figPath);

end