function [d,h,fnames] = read_all_p(series2inc,what,verb)
%READ_ALL_P  Read multiple p-files and ScanArchives
%[d,h,fnames] = read_all_p(series2inc,what)
% series2inc   Series numbers to include
%       what   Load: 0=pfiles+ScanArchives,1=pfiles,2=ScanArchives      (0)
%       verb   Verbose mode                                          (true)
%          d   Data (in cell array)
%          h   Header structures (in cell array)
%
%  11/2023 Rolf Schulte
if ((nargin<1) && (nargout<1)), help(mfilename); return; end

if ~exist('series2inc','var'), series2inc = []; end
if ~exist('what','var'), what = []; end
if isempty(what),        what = 0; end
if ~exist('verb','var'), verb = []; end
if isempty(verb),        verb = true; end
if isunix
    dsep = '/';
else
    dsep = '\';
end
use_ox = false;


%% read header
switch what
    case 0
        flist = [dir('**/*.h5') ; dir('**/*.7')];
        suffix = '*.h5 & *.7';
    case 1
        flist = dir('**/*.7');
        suffix = '*.7';
    case 2
        flist = dir('**/*.h5');
        suffix = '*.h5';
    otherwise, error('what(=%d) unknown',what);
end
if verb
    fprintf('---\n%d files with suffix %s found\n',length(flist),suffix);
end
iinc = [];
for ll=1:length(flist)
    if isempty(regexpi(flist(ll).name,'noisestatistics','once'))
        iinc = [iinc,ll];
    end
end
if verb && (length(flist)>length(iinc))
    fprintf('---\nRemoving noisestatistics: %d files remaining\n',length(flist)); 
end
flist = flist(iinc);
iinc = [];
for ll=1:length(flist)
    fname = [flist(ll).folder dsep flist(ll).name];
    if exist(fname,'file')
        if verb, fprintf('---\nReading header %s\n',fname); end
        if ~isempty(regexpi(fname,'\.7$','once'))
            h{ll} = read_MR_headers(fname);
        else
            if ~isempty(regexpi(fname,'\.h5$','once'))
                if use_ox
                    h{ll} = read_archive_header(fname);
                else
                    h{ll} = read_arc_header(fname);
                end
            else
                warning('fname(=%s) not ending with .7 or .h5',fname);
            end
        end
        if isempty(h{ll})
            warning('h{%d} empty: could not read %s',ll,fname);
        else
            iinc = [iinc,ll];
        end
    else
        warning('strange input'); disp(fname);
    end
end
if ~exist('h','var'), error('no files found'); end

%% remove empty headers
if length(iinc)<length(flist)
    if verb
        fprintf('---\nRemoving empty headers: %d files remaining\n',length(iinc));
    end
    flist = flist(iinc);
    for ll=1:length(iinc)
        htmp{ll} = h{iinc(ll)};
    end
    h = htmp;
    clear htmp
end


%% sort according to time when data was acquired
if verb, fprintf('---\nSorting files according to acquired time\n'); end
for ll=1:length(h), tt(ll) = h{ll}.image.im_datetime; end
[~,ii] = sort(tt);
flist = flist(ii);
for ll=1:length(h), tmp{ll} = h{ii(ll)}; end
h = tmp;


%% check for data to include
fnames{1} = '';
inc = 0;
do_inc = true;
for ll=1:length(h)
    if ~isempty(series2inc)
        series  = h{ll}.exam.ex_numseries;
        fprintf('---\nl=%d series=%d\n',ll,series);
        if any(series==series2inc)
            do_inc = true;
        else
            do_inc = false;
        end
    end
    if do_inc
        inc = inc+1;
        fnames{inc} = [flist(ll).folder dsep flist(ll).name];
    end
end
if inc==0, warning('no series found'); end
if verb, fprintf('---\n%d files included\n',length(fnames)); end


%% load data
if nargout>1
    if verb, fprintf('---\nReading %d files\n',length(fnames)); end
    clear h
    for ll=1:length(fnames)
        [d{ll},h{ll}] = read_p(fnames{ll});
    end
end

end      % get_pfname.m
