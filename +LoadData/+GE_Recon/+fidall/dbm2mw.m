function mW = dbm2mw(dBm)
%DBM2WATT Convert dBm values to mW
% mW = dbm2mw(dBm)
% Formula: mW = 10^(dBm/10);
% 12/202 Rolf Schulte
if (nargin<1), help(mfilename); return; end
mW = 10.^(dBm/10);

end      % dbm2mw.m
