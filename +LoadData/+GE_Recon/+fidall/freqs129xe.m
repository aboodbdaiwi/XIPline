function [f,f_str,f_ppm] = freqs129xe(b0)
%FREQS129XE  Frequencies of dissolved-phase Xenon (sensitive to oxygenation)
%[f,f_str,f_ppm] = freqs13c(b0)
%     b0  Main magnetic field strength          [T]  (3)
%      f  Frequencies                           [Hz]
%  f_str  Metabolite names
%  f_ppm  Frequencies                           [ppm]
%
%  2/2023 Rolf Schulte
if ~exist('b0','var'), b0 = []; end
if isempty(b0), b0 = 3; end

f = [35331165 35338869 35338162];
f_ppm = f/abs(gyrogamma(129))/3*2*pi*1d6;
f_ppm = f_ppm-f_ppm(1);
f = (f-f(1))/3*b0;
f_str = {'lung','blood','tissue'};

end      % freqs129xe.m
