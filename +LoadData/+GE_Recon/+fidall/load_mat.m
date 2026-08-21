function mm = load_mat(fdir)
%LOAD_MAT  Search for mat file in directory and load it
%  mm = load_mat()
% 11/2022 Rolf Schulte
% if nargin<1, help(mfilename); return; end

if isunix, dirsep = '/';
else,      dirsep = '\'; end
if ~exist('fdir','var'), fdir = []; end
if isempty(fdir)
    xx = dir('*.mat');
else
    if strcmp(fdir(end),'/') || strcmp(fdir(end),'\')
        fdir = fdir(1:end-1); 
    end
    xx = dir([fdir dirsep '*.mat']);
end

switch length(xx)
    case 0, error('no files found');
    case 1
    otherwise
        warning('multiple mat files found; loading first one');
end
fname = [xx(1).folder dirsep xx(1).name];
fprintf('Loading file\n''%s''\n',fname);


%% checks
if ~exist(fname,'file')
    error('file(=''%s'') not found',fname);
end
if nargout<1, warning('no output argument'); end


%% load mat file
mm = load(fname);


end      % load_mat.m
