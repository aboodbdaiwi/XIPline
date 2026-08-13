function filtered_inputs = filter_initial_inputs( ...
    initial_inputs, startDate, endDate, polarizer, selectSpecific, rowNumber)

    % ---------------------------------------------------------------
    % Select specific row
    % ---------------------------------------------------------------
    if selectSpecific

        % Convert text to number if RowEditField is a text field
        if ischar(rowNumber) || isstring(rowNumber)
            rowNumber = str2double(rowNumber);
        end

        if ~isnumeric(rowNumber) || isempty(rowNumber) || isnan(rowNumber)
            error('Please enter a valid numeric Row#.');
        end

        rowNumber = round(rowNumber);

        if rowNumber < 1 || rowNumber > height(initial_inputs)
            error('Row# must be between 1 and %d.', height(initial_inputs));
        end

        % Directly select requested row
        filtered_inputs = initial_inputs(rowNumber, :);

        return;
    end

    % ---------------------------------------------------------------
    % Start with all data
    % ---------------------------------------------------------------
    filtered_inputs = initial_inputs;

    % ---------------------------------------------------------------
    % Dates
    % ---------------------------------------------------------------
    scanDates = filtered_inputs.Date;

    if ~isdatetime(scanDates)
        scanDates = datetime(scanDates);
    end

    if ~isdatetime(startDate)
        startDate = datetime(startDate);
    end

    if ~isdatetime(endDate)
        endDate = datetime(endDate);
    end

    scanDates = dateshift(scanDates, 'start', 'day');
    startDate = dateshift(startDate, 'start', 'day');
    endDate   = dateshift(endDate, 'start', 'day');

    if startDate > endDate
        error('Start Date cannot be later than End Date.');
    end

    % Date filtering
    keepRows = scanDates >= startDate & scanDates <= endDate;

    % ---------------------------------------------------------------
    % Polarizer filtering
    % ---------------------------------------------------------------
    if ~strcmpi(string(polarizer), "All")

        dataPolarizer = strtrim(string(filtered_inputs.Polarizer));
        selectedPolarizer = strtrim(string(polarizer));

        keepRows = keepRows & ...
            strcmpi(dataPolarizer, selectedPolarizer);
    end

    % ---------------------------------------------------------------
    % Apply filter
    % ---------------------------------------------------------------
    filtered_inputs = filtered_inputs(keepRows, :);

end