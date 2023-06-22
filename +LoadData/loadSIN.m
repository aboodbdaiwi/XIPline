%% LOADSIN    Reads a Philips .SIN file
%
% INFO = LOADSIN(FILENAME)
%
%   FILENAME is a string containing a file prefix or name of the .SIN
%   text file, e.g. RAW_001 or RAW_001.SIN
%
%   INFO is a structure containing details from the .SIN file
%
% Example:
%  info = loadsin('raw_777.SIN');
%
%  See also: LOADLABRAW, LOADPDFXML
%
%  Dependencies: none
%

%% Revision History
% * 2012.03.07    initial version - ryanrobison
% * 2012.08.03    added .(label) to prefixcnt and valcnt - welcheb

function info = loadSIN(filename)

% Allow user to select a file if input FILENAME is not provided or is empty
if nargin < 1 || isempty(filename)
  [fn, pn] = uigetfile({'*.SIN'},'Select a SIN file');
  if fn ~= 0
    filename = sprintf('%s%s',pn,fn);
  else
    disp('LOADSIN cancelled');
    return;
  end
end

% Parse the filename.
% It may be the SIN or just the filename prefix
% Instead of REGEXP, use REGEXPI which igores case
toks = regexpi(filename,'^(.*?)(\.sin)?$','tokens');
prefix = toks{1}{1};
sinname = sprintf('%s.SIN',prefix);
info.filename = sinname;

fid = fopen(sinname,'r');
if fid~=-1,
    sintext = fread(fid,inf,'uint8=>char')';
    fclose(fid);
else
    error( sprintf('cannot open %s for reading', listname) );
end

line_expr = '(\d+)\s(\d+)\s(\d+):';
line_begin = regexp(sintext,line_expr);
cell_expr = '[^:\s]+';
for ind = 1:length(line_begin)
    start_idx = line_begin(ind);
    if ind == length(line_begin)
        end_idx = length(sintext);
    else
        end_idx = line_begin(ind+1)-1;
    end
    line{ind} = sintext(start_idx:end_idx);
    entries{ind} = regexp(line{ind},cell_expr,'match');
    num_entries = length(entries{ind});
    label = char(entries{ind}(4));
    if ~isfield(info,(label))
        info.(label) = [];
    end

    if ~isfield(info.(label),'prefix')
        info.(label).prefix = [];
        prefcnt.(label) = 1;
    end
    info.(label).prefix{prefcnt.(label)} = entries{ind}(1:3);
    prefcnt.(label) = prefcnt.(label) + 1;
    
    if ~isfield(info.(label),'vals')
        info.(label).vals = [];
        valcnt.(label) = 1;
    end
    if(ind > 1)
        info.(label).vals(1:num_entries-4,valcnt.(label)) = str2double(entries{ind}(5:num_entries));
        valcnt.(label) = valcnt.(label) + 1;
    else
        info.(label).vals = entries{ind}(5);
    end
end