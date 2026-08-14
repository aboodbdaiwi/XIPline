function exportStructToJSON(dataStruct, outputJSONFile)
%EXPORTSTRUCTTOJSON Flatten a MATLAB structure and export to JSON.
%
% Handles:
%   char
%   string
%   numeric
%   logical
%   datetime
%   duration
%   categorical
%   cell
%   table
%   struct
%
% Large numeric arrays are represented by their dimensions instead of
% being written completely to JSON.

    if ~isstruct(dataStruct)
        error('Input must be a structure.');
    end

    flatMap = containers.Map();

    % ================================================================
    % Recursive structure flattening
    % ================================================================
    function flattenFields(s, prefix)

        fields = fieldnames(s);

        for i = 1:numel(fields)

            key = fields{i};
            val = s.(key);

            if isempty(prefix)
                fullKey = key;
            else
                fullKey = [prefix '.' key];
            end

            % --------------------------------------------------------
            % Skip fields that should not be exported
            % --------------------------------------------------------
            if isa(val, 'matlab.ui.Figure')
                continue;
            end

            if strcmp(key, 'DicomInfo') && isstruct(val)
                continue;
            end

            % --------------------------------------------------------
            % STRUCT
            % --------------------------------------------------------
            if isstruct(val)

                if isscalar(val)
                    flattenFields(val, fullKey);
                else
                    flatMap(fullKey) = sprintf( ...
                        '[struct array %s]', ...
                        mat2str(size(val)));
                end

            % --------------------------------------------------------
            % CHAR
            % --------------------------------------------------------
            elseif ischar(val)

                flatMap(fullKey) = val;

            % --------------------------------------------------------
            % STRING
            % --------------------------------------------------------
            elseif isstring(val)

                if isempty(val)
                    flatMap(fullKey) = '';

                elseif isscalar(val)

                    if ismissing(val)
                        flatMap(fullKey) = '';
                    else
                        flatMap(fullKey) = char(val);
                    end

                else
                    temp = val;
                    temp(ismissing(temp)) = "";

                    flatMap(fullKey) = cellstr(temp);
                end

            % --------------------------------------------------------
            % CATEGORICAL
            % --------------------------------------------------------
            elseif iscategorical(val)

                if isempty(val)
                    flatMap(fullKey) = '';

                elseif isscalar(val)

                    if isundefined(val)
                        flatMap(fullKey) = '';
                    else
                        flatMap(fullKey) = char(string(val));
                    end

                else
                    temp = string(val);
                    temp(ismissing(temp)) = "";

                    flatMap(fullKey) = cellstr(temp);
                end

            % --------------------------------------------------------
            % DATETIME
            % --------------------------------------------------------
            elseif isdatetime(val)

                if isempty(val)
                    flatMap(fullKey) = '';

                elseif isscalar(val)

                    if isnat(val)
                        flatMap(fullKey) = '';
                    else
                        flatMap(fullKey) = char(string(val));
                    end

                else
                    temp = string(val);
                    temp(ismissing(temp)) = "";

                    flatMap(fullKey) = cellstr(temp);
                end

            % --------------------------------------------------------
            % DURATION
            % --------------------------------------------------------
            elseif isduration(val)

                if isempty(val)
                    flatMap(fullKey) = '';

                elseif isscalar(val)
                    flatMap(fullKey) = char(string(val));

                else
                    flatMap(fullKey) = cellstr(string(val));
                end

            % --------------------------------------------------------
            % NUMERIC / LOGICAL
            % --------------------------------------------------------
            elseif isnumeric(val) || islogical(val)

                if isempty(val)

                    flatMap(fullKey) = [];

                elseif isscalar(val)

                    if isnumeric(val) && isnan(val)
                        flatMap(fullKey) = [];
                    else
                        flatMap(fullKey) = val;
                    end

                elseif numel(val) <= 100

                    % Keep small arrays
                    flatMap(fullKey) = val;

                else

                    % Do not write large imaging arrays
                    flatMap(fullKey) = sprintf( ...
                        '[%s]', ...
                        strtrim(sprintf('%d ', size(val))));
                end

            % --------------------------------------------------------
            % CELL
            % --------------------------------------------------------
            elseif iscell(val)

                if isempty(val)

                    flatMap(fullKey) = [];

                elseif isscalar(val)

                    flatMap(fullKey) = convertValue(val{1});

                elseif numel(val) <= 100

                    temp = cell(size(val));

                    for j = 1:numel(val)
                        temp{j} = convertValue(val{j});
                    end

                    flatMap(fullKey) = temp;

                else

                    flatMap(fullKey) = sprintf( ...
                        '[cell %s]', ...
                        mat2str(size(val)));
                end

            % --------------------------------------------------------
            % TABLE
            % --------------------------------------------------------
            elseif istable(val)

                if height(val) <= 20 && width(val) <= 20

                    tempStruct = table2struct(val);

                    try
                        flatMap(fullKey) = tempStruct;
                    catch
                        flatMap(fullKey) = sprintf( ...
                            '[table %d x %d]', ...
                            height(val), width(val));
                    end

                else

                    flatMap(fullKey) = sprintf( ...
                        '[table %d x %d]', ...
                        height(val), width(val));
                end

            % --------------------------------------------------------
            % FUNCTION HANDLE
            % --------------------------------------------------------
            elseif isa(val, 'function_handle')

                flatMap(fullKey) = func2str(val);

            % --------------------------------------------------------
            % OTHER OBJECTS
            % --------------------------------------------------------
            elseif isobject(val)

                try
                    flatMap(fullKey) = char(string(val));
                catch
                    flatMap(fullKey) = sprintf( ...
                        '[%s]', ...
                        class(val));
                end

            % --------------------------------------------------------
            % FALLBACK
            % --------------------------------------------------------
            else

                try
                    flatMap(fullKey) = char(string(val));
                catch
                    flatMap(fullKey) = sprintf( ...
                        '[Unsupported type: %s]', ...
                        class(val));
                end

            end
        end
    end

    % ================================================================
    % Convert individual values, especially values inside cells
    % ================================================================
    function out = convertValue(val)

        if isempty(val)

            out = [];

        elseif ischar(val)

            out = val;

        elseif isstring(val)

            if isscalar(val)
                if ismissing(val)
                    out = '';
                else
                    out = char(val);
                end
            else
                out = cellstr(val);
            end

        elseif iscategorical(val)

            if isscalar(val)
                if isundefined(val)
                    out = '';
                else
                    out = char(string(val));
                end
            else
                out = cellstr(string(val));
            end

        elseif isdatetime(val)

            if isscalar(val)
                if isnat(val)
                    out = '';
                else
                    out = char(string(val));
                end
            else
                out = cellstr(string(val));
            end

        elseif isnumeric(val) || islogical(val)

            if isscalar(val) && isnumeric(val) && isnan(val)
                out = [];
            elseif numel(val) <= 100
                out = val;
            else
                out = sprintf( ...
                    '[%s]', ...
                    strtrim(sprintf('%d ', size(val))));
            end

        elseif iscell(val)

            out = cell(size(val));

            for k = 1:numel(val)
                out{k} = convertValue(val{k});
            end

        else

            try
                out = char(string(val));
            catch
                out = sprintf('[%s]', class(val));
            end
        end
    end

    % ================================================================
    % Flatten top-level structure
    % ================================================================
    flattenFields(dataStruct, '');

    % ================================================================
    % Convert map to structure
    % ================================================================
    jsonStruct = struct();

    keys = flatMap.keys;

    for i = 1:numel(keys)

        validKey = matlab.lang.makeValidName(keys{i});

        jsonStruct.(validKey) = flatMap(keys{i});
    end

    % ================================================================
    % Encode JSON
    % ================================================================
    jsonStr = jsonencode( ...
        jsonStruct, ...
        'PrettyPrint', true);

    % ================================================================
    % Write JSON
    % ================================================================
    fid = fopen(outputJSONFile, 'w');

    if fid == -1
        error( ...
            'Cannot open file for writing: %s', ...
            outputJSONFile);
    end

    cleaner = onCleanup(@() fclose(fid));

    fwrite(fid, jsonStr, 'char');

    fprintf('Exported structure to %s\n', outputJSONFile);

end