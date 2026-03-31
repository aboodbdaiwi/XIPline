function T = LoadTable(xlsxFile)
    import VentilationFunctions.*

    fmt = detectSpreadsheetFormat(xlsxFile);

    switch fmt
        case "old"
            T = parseOldFormat(xlsxFile);
        case "new"
            T = parseNewFormat(xlsxFile);
        otherwise
            error("Unknown spreadsheet format: %s", fmt);
    end

   
end