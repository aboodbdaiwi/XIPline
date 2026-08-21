function [f,f_str,f_ppm] = freqs2h(b0)
%FREQS2H  Frequencies of [6-6'-2H2]glucose and metabolites
%         HDO, [6,6'-2H2]Glc, [4-2H]Glx, [3-2H]Lac
%[f,f_str,f_ppm] = freqs2h(b0)
%     b0  Main magnetic field strength          [T]  (3)
%      f  Frequencies                           [Hz]
%  f_str  Metabolite names
%  f_ppm  Frequencies                           [ppm]
%
%  4/2024 Rolf Schulte
if ~exist('b0','var'), b0 = []; end
if isempty(b0), b0 = 3; end

f_ppm = [4.7 3.7 2.1 1.3];
f_str = {'HDO','Glc','Glx','Lac'};
f = gyrogamma(2)*b0/2/pi*1d-6*f_ppm;

end      % freqs2h.m
