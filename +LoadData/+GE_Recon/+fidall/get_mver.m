function mver = get_mver
%GET_MVER Extract Matlab Version as double
% mver = get_mver
%  5/2023 Rolf Schulte

mver = version;
ii = regexp(mver,'\.');
if length(ii)<2, warning('length(ii)(=%g)<2',length(ii)); end
mver = mver(1:(ii(2)-1));
if isempty(mver),  warning('isempty(mver)'); end
mver = str2double(mver);
if isempty(mver),  warning('isempty(mver)'); end

end      % get_mver.m
