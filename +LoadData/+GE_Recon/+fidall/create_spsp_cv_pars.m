function create_spsp_cv_pars(freq,flip,fname)
%CREATE_SPSP_CV_PARS  Generate + save cv_parsXX.cvs file
% create_spsp_cv_pars(freq,flip,fname)
%   freq  Frequencies                                 [Hz]      ([])
%   flip  Flip angles                                 [deg]     ([])
%  fname  File name   (if [] -> stdout)                         ([])
%
% 11/2020 Rolf Schulte
if nargin<1, help(mfilename); return; end

%% input parameters
if ~exist('freq','var'),  freq = []; end
if ~exist('flip','var'),  flip = []; end
if ~exist('fname','var'), fname = []; end


%% open file or revert to stdout
if ~isempty(fname)
    fid = fopen(fname,'w');
else
    fid = 1;
end


%% print spsp parameters
if ~isempty(freq)
    fprintf(fid,'spsp_nfreqs %d\n',length(freq));
    for l=1:length(freq)
        fprintf(fid,'spsp_freq%d %g\n',l,freq(l));
    end
end
if ~isempty(flip)
    fprintf(fid,'spsp_nflips %d\n',length(flip));
    for l=1:length(flip)
        fprintf(fid,'spsp_flip%d %g\n',l,flip(l));
    end
end


%% close file
if ~isempty(fname), fclose(fid); end


end      % create_spsp_cv_pars.m
