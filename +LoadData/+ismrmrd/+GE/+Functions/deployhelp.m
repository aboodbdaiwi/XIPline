function deployhelp(varargin)
%DEPLOYHELP Display help files in matlab compiled code (interactive mode)

pth = [getenv('homedir') '/matlab/'];
fname = strtrim(varargin{1});     % remove whitespace
mfname = [pth fname];
if isempty(regexpi(mfname,'\.m')), mfname = [mfname '.m']; end

if ~exist(mfname,'file'), 
    warning('matlab file not found: help %s',mfname); 
    return;
end
fid = fopen(mfname);

while 1,
    tline = fgetl(fid);
    if ~isempty(tline),
        if (tline(1) == '%'), 
            fprintf('%s\n',tline(2:end));
            break; 
        end
    end
end

while 1,
    tline = fgetl(fid);
    if isempty(tline), fclose(fid); break; end
    if (tline(1) ~= '%'), 
        fclose(fid);
        break;
    end
    fprintf('%s\n',tline(2:end));
end


% %HELP Interface to type help text from /export/home/sdc/fidall/helpfiles
% type(['/export/home/sdc/fidall/helpfiles/' varargin{1} '.help']);
