function T = normalizeStudyColumn(T)
import VentilationFunctions.*

map = studyDescriptionMap();

study = string(T.Study);
study = strtrim(study);

for i = 1:numel(study)

    s = study(i);

    % If already a code like IRC1046 leave it alone
    if ~isempty(regexp(s,'^[A-Z]{3,4}\d{3,4}$','once'))
        continue
    end

    if isKey(map,s)
        study(i) = map(s);
    else
        error("Unknown study description: %s", s)
    end

end

T.Study = cellstr(study);

end