function exportStructToJSON(dataStruct, outputJSONFile)
    % Check that the input is a structure
    if ~isstruct(dataStruct)
        error('Input must be a structure.');
    end

    % Use containers.Map for storing flattened key-value pairs
    flatMap = containers.Map();

    % Recursive helper to flatten structure fields
    function flattenFields(s, prefix)
        fields = fieldnames(s);
        for i = 1:length(fields)
            key = fields{i};
            val = s.(key);
            fullKey = key;
            if ~isempty(prefix)
                fullKey = [prefix '.' key];
            end

            % Skip figures and specific fields
            if isa(val, 'matlab.ui.Figure') || (strcmp(key, 'DicomInfo') && isstruct(val))
                continue;
            end

            % Handle value types
            if ischar(val) || isstring(val)
                flatMap(fullKey) = char(val);

            elseif isnumeric(val) || islogical(val)
                if isscalar(val)
                    flatMap(fullKey) = val;
                elseif numel(val) <= 100
                    flatMap(fullKey) = val(:)';  % row vector
                else
                    flatMap(fullKey) = sprintf('[%s]', num2str(size(val)));
                end

            elseif isstruct(val)
                flattenFields(val, fullKey);

            else
                flatMap(fullKey) = 'Unsupported Type';
            end
        end
    end

    % Flatten all fields from top-level
    flattenFields(dataStruct, '');

    % Convert map to struct for JSON encoding
    jsonStruct = struct();
    keys = flatMap.keys;
    for i = 1:numel(keys)
        jsonStruct.(matlab.lang.makeValidName(keys{i})) = flatMap(keys{i});
    end

    % Encode and write
    jsonStr = jsonencode(jsonStruct, 'PrettyPrint', true);
    fid = fopen(outputJSONFile, 'w');
    if fid == -1
        error('Cannot open file for writing: %s', outputJSONFile);
    end
    fwrite(fid, jsonStr, 'char');
    fclose(fid);

    fprintf('Exported structure to %s\n', outputJSONFile);
end
