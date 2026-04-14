function excelFile = resolve_excel_file()
%RESOLVE_EXCEL_FILE Find newest VDP_Inputs*.xlsx in Downloads.
% Falls back to file picker if none are found.
% If multiple matches exist, keeps the newest and deletes the older ones.

    downloadsDir = fullfile(getenv('USERPROFILE'), 'Downloads');

    if ~isfolder(downloadsDir)
        excelFile = prompt_for_excel_file();
        return
    end

    listing = dir(fullfile(downloadsDir, 'VDP_Inputs*.xlsx'));

    % Keep only:
    %   VDP_Inputs.xlsx
    %   VDP_Inputs (N).xlsx
    keep = false(numel(listing), 1);
    for k = 1:numel(listing)
        keep(k) = ~isempty(regexp(listing(k).name, '^VDP_Inputs( \(\d+\))?\.xlsx$', 'once'));
    end
    listing = listing(keep);

    if isempty(listing)
        excelFile = prompt_for_excel_file();
        return
    end

    if numel(listing) == 1
        excelFile = fullfile(listing(1).folder, listing(1).name);
        fprintf('Using Excel file:\n%s\n\n', excelFile);
        return
    end

    % Multiple matches: choose newest
    [~, idx] = max([listing.datenum]);
    latest = listing(idx);
    excelFile = fullfile(latest.folder, latest.name);

    fprintf('Multiple matching files found. Using most recent:\n%s\n\n', excelFile);

    % Delete older matches
    listing(idx) = [];
    fprintf('Deleting old VDP input files:\n');
    for k = 1:numel(listing)
        fpath = fullfile(listing(k).folder, listing(k).name);
        try
            delete(fpath);
            fprintf('  Deleted: %s\n', fpath);
        catch ME
            warning('VDPInputs:DeleteFailed', ...
                'Could not delete file: %s\nReason: %s', ...
                fpath, ME.message);
        end
    end
    fprintf('\n');
end


function excelFile = prompt_for_excel_file()
    [f, p] = uigetfile({'*.xlsx', 'Excel Files (*.xlsx)'}, ...
        'Select VDP input Excel file');

    if isequal(f, 0)
        error('VDPInputs:NoExcelSelected', ...
            'No Excel file was selected.');
    end

    excelFile = fullfile(p, f);

end