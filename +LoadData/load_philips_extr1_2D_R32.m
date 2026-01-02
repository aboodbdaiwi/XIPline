function [ImgK,NoiK,kx_oversample_factor] = load_philips_extr1_2D_R32(filename,ch_range)
%% ========================================================================
% Robust Philips .list/.data reader for Cartesian k-space (STD vectors)
% Fixes: header-order dependence + bad numeric parsing (str2double) +
%        fragile fixed-column substring parsing.
%
% INPUTS:
%   filename  : full path WITHOUT extension (e.g., '...\sub-0003_ses-20141231_acq-1')
%   ch_range  : [ch_min ch_max] in 1-based coil indices you want to keep
%               Example: ch_range = [1 32];
%
% OUTPUTS:
%   ImgK      : [ky, kx, sl, ch, extr1] complex single
%   NoiK      : (not populated in this script; returned as [])
%   kx_oversample_factor : parsed from header if present (otherwise [])
% ========================================================================


FileNameData = [filename, '.data'];
FileNameList = [filename, '.list'];

assert(isfile(FileNameList), "Missing .list file: %s", FileNameList);
assert(isfile(FileNameData), "Missing .data file: %s", FileNameData);

% -------------------- Settings -------------------------------------------
DataType           = 'float32';
DataBytesPerType   = 4;        % float32
isLittleEndian     = true;     % typical Philips exports
machinefmt         = ternary(isLittleEndian, 'ieee-le', 'ieee-be');

% -------------------- Read full .list into memory ------------------------
lines = readlines(FileNameList);
lines = string(lines);

% -------------------- Parse header robustly (regex) ----------------------
hdr = struct();

hdr.sl_size     = parseScalarAfterColon(lines, "number_of_locations");
hdr.extr1_size  = parseScalarAfterColon(lines, "number_of_extra_attribute_1_values");

hdr.kx_range    = parseRangeAfterColon(lines, "kx_range");   % [start end]
hdr.ky_range    = parseRangeAfterColon(lines, "ky_range");   % [start end]

% Optional (not required for indexing here)
hdr.kx_oversample_factor = parseScalarAfterColon(lines, "kx_oversample_factor", true);

% Validate required header fields
req = ["sl_size","extr1_size","kx_range","ky_range"];
for k = 1:numel(req)
    if ~isfield(hdr, req(k)) || isempty(hdr.(req(k)))
        error("Header parse failed: missing '%s' in .list header.", req(k));
    end
end

kx_size = hdr.kx_range(2) - hdr.kx_range(1) + 1;
ky_size = hdr.ky_range(2) - hdr.ky_range(1) + 1;

% Channel selection bookkeeping
ch_size = ch_range(2) - ch_range(1) + 1;

% Preallocate output
ImgK = complex(zeros(ky_size, kx_size, hdr.sl_size, ch_size, hdr.extr1_size, 'single'));

% -------------------- Open data file -------------------------------------
fidData = fopen(FileNameData, 'r', machinefmt);
if fidData < 0, error("Could not open .data: %s", FileNameData); end

% -------------------- Find STD vector index section ----------------------
% STD lines begin with "  STD"
isSTD = startsWith(lines, "  STD");

if ~any(isSTD)
    fclose(fidData);
    error("No STD lines found in .list. Nothing to read.");
end

% -------------------- Read STD vectors -----------------------------------
nSTD = 0;
nSkipOutOfRange = 0;
nSkipBadRead    = 0;

for i = find(isSTD).'
    line = lines(i);

    % Parse numbers after the "STD" token
    % Format (per .list header):
    % typ mix dyn card echo loca chan extr1 extr2 ky kz n.a aver sign rf grad enc rtop rr size offset
    nums = sscanf(extractAfter(line, "STD"), '%f').';
    if numel(nums) < 20
        nSkipBadRead = nSkipBadRead + 1;
        continue;
    end

    loca  = nums(5);
    chan  = 0; %nums(2); force to be 1 channel only 
    extr1 = nums(7);
    ky    = nums(9);

    sizeBytes = nums(end-1);
    offset    = nums(end);

    % Convert to 1-based indices
    slIndex    = loca + 1;
    chIndexAbs = chan + 1; % absolute coil index (1-based)
    extr1Index = extr1 + 1;
    kyIndex    = ky - hdr.ky_range(1) + 1;

    % Channel range filter
    if chIndexAbs < ch_range(1) || chIndexAbs > ch_range(2)
        nSkipOutOfRange = nSkipOutOfRange + 1;
        continue;
    end
    chIndex = chIndexAbs - ch_range(1) + 1;

    % Bounds check
    if kyIndex < 1 || kyIndex > ky_size || ...
       slIndex < 1 || slIndex > hdr.sl_size || ...
       extr1Index < 1 || extr1Index > hdr.extr1_size
        nSkipOutOfRange = nSkipOutOfRange + 1;
        continue;
    end

    % Seek & read
    fseek(fidData, offset, 'bof');
    nFloat = sizeBytes / DataBytesPerType;

    temp = fread(fidData, nFloat, DataType);

    % Expect complex interleaved: [Re Im Re Im ...]
    if numel(temp) < 2*kx_size
        nSkipBadRead = nSkipBadRead + 1;
        continue;
    end

    cvec = complex(single(temp(1:2:2*kx_size)), single(temp(2:2:2*kx_size)));
    ImgK(kyIndex, :, slIndex, chIndex, extr1Index) = cvec(:).';  % row into kx

    nSTD = nSTD + 1;
end

fclose(fidData);

% -------------------- Report ---------------------------------------------
fprintf("\n==== Philips .list/.data read summary ====\n");
fprintf("kx_range = [%d %d] => kx_size = %d\n", hdr.kx_range(1), hdr.kx_range(2), kx_size);
fprintf("ky_range = [%d %d] => ky_size = %d\n", hdr.ky_range(1), hdr.ky_range(2), ky_size);
fprintf("sl_size  = %d\n", hdr.sl_size);
fprintf("extr1_size = %d\n", hdr.extr1_size);
fprintf("ch_range = [%d %d] => ch_size = %d\n", ch_range(1), ch_range(2), ch_size);
fprintf("STD lines processed/filled = %d\n", nSTD);
fprintf("Skipped (out of range)     = %d\n", nSkipOutOfRange);
fprintf("Skipped (bad read/parse)   = %d\n\n", nSkipBadRead);

% Outputs requested by function signature
NoiK = [];  % not populated by this script
kx_oversample_factor = hdr.kx_oversample_factor;

% ======================= Helper functions ===============================
function val = parseScalarAfterColon(lines, key, optional)
    if nargin < 3, optional = false; end
    pat = "^\.\s*0\s+0\s+0\s+" + key + "\s*:\s*([-\d\.eE\+]+)\s*$";
    m = regexp(lines, pat, 'tokens', 'once');
    idx = find(~cellfun(@isempty, m), 1, 'first');
    if isempty(idx)
        if optional, val = []; return; end
        val = [];
        return;
    end
    val = str2double(m{idx}{1});
end

function rng = parseRangeAfterColon(lines, key)
    pat = "^\.\s*0\s+0\s+0\s+" + key + "\s*:\s*([-\d]+)\s+([-\d]+)\s*$";
    m = regexp(lines, pat, 'tokens', 'once');
    idx = find(~cellfun(@isempty, m), 1, 'first');
    if isempty(idx)
        rng = [];
        return;
    end
    rng = [str2double(m{idx}{1}), str2double(m{idx}{2})];
end

function out = ternary(cond, a, b)
    if cond, out = a; else, out = b; end
end
end
