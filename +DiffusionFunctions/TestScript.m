excelFile = "C:\Users\MCM5BK\OneDrive - cchmc\Documents\03_Data Analysis\02_Data Logs\WorkflowOverhaul2026\main_rawDiff.xlsx";

C = DiffusionFunctions.standardize_CCHMC_DiffInputTable(excelFile);

raw = readcell(excelFile);
data = raw(2:end,:);   % match adapter output rows

Old = struct();

% Original rawDiff column mapping from main script, shifted to data rows
Old.StudyCol     = data(:,2);
Old.SubjectCol   = data(:,3);
Old.ScanDateCol  = data(:,4);
Old.SubNumCol    = data(:,5);
Old.DiffFileCol  = data(:,6);
Old.ACQ_TypeCol  = data(:,7);
Old.ScannerSW    = data(:,9);
Old.ScannerCol   = repmat({''}, size(data,1), 1);
Old.SexCol       = data(:,10);
Old.AgeCol       = data(:,11);
Old.DiseaseCol   = data(:,13);
Old.NoteCol      = data(:,15);
Old.ImageQCol    = data(:,16);
Old.MaskCol      = data(:,18);
Old.RuneCol      = data(:,20);

% rawDiff branch did not originally define SerNum
Old.SerNum       = repmat({NaN}, size(data,1), 1);

fields = fieldnames(C);

result = table('Size', [numel(fields), 4], ...
    'VariableTypes', {'string','double','double','string'}, ...
    'VariableNames', {'Field','NumMismatches','FirstMismatchRow','Example'});

for k = 1:numel(fields)
    f = fields{k};

    a = normalizeForCompare(C.(f));
    b = normalizeForCompare(Old.(f));

    if numel(a) ~= numel(b)
        error("Size mismatch in %s: C has %d rows, Old has %d rows", ...
            f, numel(a), numel(b));
    end

    mismatch = a ~= b;

    result.Field(k) = string(f);
    result.NumMismatches(k) = sum(mismatch);

    if any(mismatch)
        firstIdx = find(mismatch, 1);
        result.FirstMismatchRow(k) = firstIdx;
        result.Example(k) = sprintf("C='%s' | Old='%s'", a(firstIdx), b(firstIdx));
    else
        result.FirstMismatchRow(k) = NaN;
        result.Example(k) = "";
    end
end

disp(result)

function s = normalizeForCompare(x)

    if iscell(x)
        s = strings(size(x));
        for ii = 1:numel(x)
            v = x{ii};

            if ismissing(v)
                s(ii) = "";
            elseif isempty(v)
                s(ii) = "";
            elseif isnumeric(v)
                if isscalar(v) && isnan(v)
                    s(ii) = "";
                else
                    s(ii) = string(v);
                end
            elseif isdatetime(v)
                s(ii) = string(v, "yyyyMMdd");
            else
                s(ii) = strtrim(string(v));
            end
        end
    else
        if isdatetime(x)
            s = string(x, "yyyyMMdd");
        else
            s = strtrim(string(x));
        end
    end

    s(ismissing(s)) = "";
    s = s(:);
end