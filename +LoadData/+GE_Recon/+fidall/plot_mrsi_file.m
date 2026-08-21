function plot_mrsi_file(fname,do_uiwait)
%PLOT_MRSI_FILE  Load MRSI data from file and start plot_mrsi
% plot_mrsi_file(fname,do_uiwait)
%     fname  MRSI mat Filename (spec,par)
% do_uiwait  Wait for figure to close before continuing    (true)
%
% 7/2020 Rolf Schulte
if nargin<1, help(mfilename); return; end

%% input variables
if ~ischar(fname),            error('~ischar(fname)'); end
if ~exist(fname,'file'),      error('fname(=''%s'') not found',fname); end
if ~exist('do_uiwait','var'), do_uiwait = []; end
if isempty(do_uiwait),        do_uiwait = true; end


%% loading MRSI mat file
xx = load(fname);
if ~isfield(xx,'spec')
    error('spec not found'); 
else
    spec = xx.spec;
end
if ~isfield(xx,'par')
    warning('par structure not found');
    par = struct;
else
    par = xx.par;
end


%% execute plotting
h = plot_mrsi(spec,par);

%% wait for figure to close before continuing
if do_uiwait, uiwait(h); end


end
