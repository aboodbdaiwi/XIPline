function T = normalizeSexColumn(T)


map = containers.Map({'M' 'F'}, {'Male' 'Female'});

sex = string(T.Sex);
sex = strtrim(sex);

for i = 1:numel(sex)

    s = sex(i);

    if ismember(s, ["Male","Female"])
        continue
    end

    if ~isempty(regexp(s, '^[MF]$', 'once'))
        sex(i) = map(char(s));
    else
        error("Unknown sex entry: %s", s)
    end

end

T.Sex = cellstr(sex);

end