function Tstd = parseNewFormat(xlsxFile)
%PARSENEWVDPFORMAT Read new spreadsheet and map to canonical schema.
import VDPInputs.*
    T = readtable(xlsxFile, "TextType", "char");

    % Rename new/raw names to canonical names
    T = renameIfPresent(T, 'SubjectID_PO',         'SubjectID_PO');
    T = renameIfPresent(T, 'Selected',             'Selected');
    T = renameIfPresent(T, 'Recon',                'Recon');
    T = renameIfPresent(T, 'RunMode',              'RunMode');
    T = renameIfPresent(T, 'Polarizer',            'Polarizer');
    T = renameIfPresent(T, 'Notes',                'OtherInternalProcessingNotes');
    T = renameIfPresent(T, 'Comments',             'Notes');
    T = renameIfPresent(T, 'PROTON_SERIES_NUMBER', 'ProtonSeriesNumber');
    T = renameIfPresent(T, 'FrequencyOffset_Hz',   'FrequencyOffsetHz');
    T = renameIfPresent(T, 'SEX',                  'Sex');
    T = renameIfPresent(T, 'XENON_SERIES_NUMBER',  'XeSeriesNumber');

    % If old-only fields do not exist, create them anyway if your
    % downstream code wants a stable schema.
    T = addMissingVar(T, "Session", strings(height(T),1));

    % Optional: if SubjectID is your preferred canonical name, derive it
    if ~ismember("SubjectID", T.Properties.VariableNames) && ...
       ismember("SubjectID_PO", T.Properties.VariableNames)
        T.SubjectID = T.SubjectID_PO;
    end

    T = coerceExcelTypes(T);
    T = normalizeStudyColumn(T);
    T = normalizeSexColumn(T);
    T = normalizeRunColumn(T);
    T = normalizeVentPath(T);
    Tstd = enforceCanonicalColumns(T);
end