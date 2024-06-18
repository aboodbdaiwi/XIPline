function phi = ernst_angle(TR,T1)
% ERNST_ANGLE Calculate SNR optimal excitation angle 
%  formula:  phi = arccos(exp(-TR/T1));
%
%  phi = ernst_angle(TR,T1)
%  phi  Optimal flip angle of excitation pulse   [deg]
%   TR  Repetition time                          [s]
%   T1  T1 relaxation time                       [s]
%
%  11/2018 Rolf Schulte
if (nargin<1), help(mfilename); return; end

phi = acos(exp(-TR/T1))*180/pi;
