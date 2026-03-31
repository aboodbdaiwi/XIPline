function T = normalizeVentPath(T)

paths = T.Vent_path;
H_paths = T.PROT_FILEPATH;

for i = 1:numel(paths)

    p = paths{i};          % char vector inside cell
    p = strtrim(p);
    h = H_paths{i};
    h = strtrim(h);


    % If already a UNC path, leave it
    if startsWith(p, '\\')
        continue
    end
    if startsWith(h, '\\')
        continue
    end

    if startsWith(p,'/\\')
        p = p(2:end);
    end
    if startsWith(h,'/\\')
        h = h(2:end);
    end

    % normalize separators
    p = strrep(p,'/','\');
    h = strrep(h,'/','\');

    

    % If it looks like a network path but missing the prefix
    if ~startsWith(p,'\\') && startsWith(p,'rds','IgnoreCase',true)
        p = ['\\' p];
    end
    if ~startsWith(h,'\\') && startsWith(h,'rds','IgnoreCase',true)
        h = ['\\' h];
    end

    paths{i} = p;
    H_paths{i} = h;

end

T.Vent_path = paths;
T.PROT_FILEPATH = H_paths;
end