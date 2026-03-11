function Tstd = parseOldFormat(xlsxFile)
%PARSEOLDVDPFORMAT Read old spreadsheet and map to canonical schema.
import VentilationFunctions.*
    T = readtable(xlsxFile, "TextType", "string");

    % Rename old/raw names to canonical names
    T = renameIfPresent(T, "Subject",               "SubjectID");
    T = renameIfPresent(T, "ScanDate",              "ScanDate");
    T = renameIfPresent(T, "Session",               "Session");
    T = renameIfPresent(T, "XENON_SERIES_NUMBER",      "XeSeriesNumber");
    T = renameIfPresent(T, "PROTON_SERIES_NUMER",   "ProtonSeriesNumber");
    T = renameIfPresent(T, "FrequencyOffset_Hz_",   "FrequencyOffsetHz");
    T = renameIfPresent(T, "Run",                   "RunMode");
    T = renameIfPresent(T, "SEX",                   "Sex");

    % Add columns that only exist in new format but that downstream code
    % may expect to always exist.
    T = addMissingVar(T, "Selected",          false(height(T),1));
    T = addMissingVar(T, "SubjectID_PO",      strings(height(T),1));
    T = addMissingVar(T, "Recon",             strings(height(T),1));
    T = addMissingVar(T, "Polarizer",         strings(height(T),1));
    T = addMissingVar(T, "Notes",             strings(height(T),1));

    % Type cleanup
    T = coerceExcelTypes(T);

    % Restrict/reorder to canonical schema
    Tstd = enforceCanonicalColumns(T);
end