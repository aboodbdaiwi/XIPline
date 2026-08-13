function outputFile = download_filtered_list(filtered_inputs)
%DOWNLOAD_FILTERED_LIST Save filtered batch-processing data to Excel.

    outputFile = '';

    % Check input
    if nargin < 1 || isempty(filtered_inputs)
        error('No filtered data available to save.');
    end

    if ~istable(filtered_inputs)
        error('Filtered data must be a MATLAB table.');
    end

    % Default filename
    defaultName = sprintf( ...
        'Filtered_Analysis_List_%s.xlsx', ...
        datestr(datetime('today'), 'yyyymmdd'));

    % Ask user where to save
    [fileName, filePath] = uiputfile( ...
        {'*.xlsx', 'Excel Files (*.xlsx)'}, ...
        'Save Filtered Analysis List', ...
        defaultName);

    % User cancelled
    if isequal(fileName, 0)
        return;
    end

    outputFile = fullfile(filePath, fileName);

    % Save
    writetable(filtered_inputs, outputFile);

    fprintf('\nFiltered data saved successfully:\n');
    fprintf('%s\n', outputFile);
    fprintf('Datasets saved: %d\n\n', height(filtered_inputs));

end