function exportStructToExcel(dataStruct, outputExcelFile)
    % Check that the input is a structure
    if ~isstruct(dataStruct)
        error('Input must be a structure.');
    end

    % Initialize output with header
    outputData = {};
    outputData{1,1} = 'VariableName';
    outputData{1,2} = 'Value(s)';
    row = 2;

    % Recursive helper function to process structure
    function row = flattenStruct(s, prefix, row)
        varNames = fieldnames(s);

        for i = 1:length(varNames)
            varName = varNames{i};

            % Skip DicomInfo structures by name
            if strcmp(varName, 'DicomInfo') && isstruct(s.(varName))
                continue;
            end

            fullVarName = [prefix, varName];
            varValue = s.(varName);

            % Skip figures
            if isa(varValue, 'matlab.ui.Figure')
                continue;
            end

            if ischar(varValue)
                outputData{row,1} = fullVarName;
                outputData{row,2} = varValue;
                row = row + 1;

            elseif isstring(varValue)
                outputData{row,1} = fullVarName;
                outputData{row,2} = char(varValue);
                row = row + 1;

            elseif islogical(varValue)
                if isscalar(varValue)
                    outputData{row,1} = fullVarName;
                    outputData{row,2} = mat2str(varValue);
                else
                    outputData{row,1} = fullVarName;
                    outputData{row,2} = sprintf('Size: [%s]', num2str(size(varValue)));
                end
                row = row + 1;

            elseif isnumeric(varValue)
                sz = size(varValue);
                numelVal = numel(varValue);

                if isscalar(varValue)
                    outputData{row,1} = fullVarName;
                    outputData{row,2} = varValue;
                    row = row + 1;

                elseif isvector(varValue) && numelVal <= 100
                    outputData{row,1} = fullVarName;
                    outputData(row, 2:(1+numelVal)) = num2cell(varValue);
                    row = row + 1;

                else
                    outputData{row,1} = fullVarName;
                    outputData{row,2} = sprintf('Size: [%s]', num2str(sz));
                    row = row + 1;
                end

            elseif isstruct(varValue)
                row = flattenStruct(varValue, [fullVarName '.'], row);

            else
                % Skip unsupported types
                continue;
            end
        end
    end

    % Start flattening from top-level
    row = flattenStruct(dataStruct, '', row);

    % Write to Excel
    writecell(outputData, outputExcelFile);
    fprintf('Exported structure to %s\n', outputExcelFile);
end
