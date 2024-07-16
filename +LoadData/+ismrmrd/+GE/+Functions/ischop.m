function chp = ischop(h,verb)
%ISCHOP  Returns true if data is chopped and false otherwise
%  chp = ischop(h,verb)
%    h   header structure
% verb   verbose                                           (default=false)
%
% 1/2021 Rolf Schulte
if nargin<1, help(mfilename); return; end

if ~exist('verb','var'), verb = []; end
if isempty(verb),        verb = false; end


%% check for fields
if isfield(h,'rdb_hdr')
    rdb_hdr = h.rdb_hdr; 
else
    if isfield(h,'RawHeader')
        rdb_hdr = h.RawHeader;
    else
        warning('h.rdb_hdr or h.RawHeader not found');
        rdb_hdr.data_collect_type = 2; % default = chopped
    end
end
if isfield(rdb_hdr,'data_collect_type')
    data_collect_type = rdb_hdr.data_collect_type;
end

%% actual determination of chopping
if sub_isodd(data_collect_type)        % chopping
    if verb, fprintf('Data not chopped\n'); end
    chp = false;
else
    if verb, fprintf('Data chopped\n'); end
    chp = true;
end

end      % ischop.m


%% sub-functions
function bool = sub_isodd(number)
if ~isfloat(number), number = double(number); end
if number-floor(number) ~= 0
  warning('Input number not integer')
end

bool = logical(number/2-floor(number/2) ~= 0);
end
