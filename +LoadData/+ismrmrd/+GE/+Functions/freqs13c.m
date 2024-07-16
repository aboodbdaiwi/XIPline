function [f,f_str] = freqs13c(ind)
%FREQS13C  Frequencies of 13C1: Lac, Pyr-Hyd, Ala, Pyr, Bi-Carbonate
% Experimentally determined (values from pnas_merrit deviate too much from
% actual values)
% [f,f_str] = freqs13c(ind)
%       f  Freqs [Hz]
%     ind  Index
%6/2010 Rolf Schulte
f = [0 -125 -215 -392 -716];
f_str = {'Lac','Pyr-Hyd','Ala','Pyr','Bi-Carb'};
if exist('ind','var'), 
    f = f(ind); 
    f_str = f_str(ind);
end
