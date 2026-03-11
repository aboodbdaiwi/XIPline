function fmt = detectSpreadsheetFormat(filename)

    headers = string(readcell(filename, "Range", "1:1"));
    headers = strip(headers);
    headers = headers(~ismissing(headers) & headers ~= "");

    % Optional normalization
    headers_norm = upper(headers);

    new_markers = ["SELECTED", "SUBJECTID_PO", "RECON", "RUNMODE", "POLARIZER"];
    old_markers = ["SESSION", "PROTON_SERIES_NUMER", "FREQUENCYOFFSET_HZ_", "RUN"];

    nNew = sum(ismember(new_markers, headers_norm));
    nOld = sum(ismember(old_markers, headers_norm));

    if nNew >= 3 && nOld == 0
        fmt = "new";
    elseif nOld >= 2 && nNew == 0
        fmt = "old";
    elseif ismember("SELECTED", headers_norm) && ismember("RUNMODE", headers_norm)
        fmt = "new";
    elseif ismember("SESSION", headers_norm) && ismember("PROTON_SERIES_NUMER", headers_norm)
        fmt = "old";
    else
        fmt = "unknown";
        error("Could not confidently detect spreadsheet schema.");
    end
end