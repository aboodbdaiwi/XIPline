function phi = ssfp_angle(T1,T2)
% SSFP_ANGLE Calculate SNR optimal excitation angle for balanced SSFP 
%  formula:  phi = arccos((T1-T2)/(T1+T2));
% Assumption: TR<<T2
%
%  phi = ssfp_angle(T1,T2)
%  phi  Optimal flip angle of excitation pulse   [deg]
%   T1  T1 relaxation time                       [s]
%   T2  T2 relaxation time                       [s]
%
%  7/2020 Rolf Schulte
if (nargin<1), help(mfilename); return; end

T1 = T1(:);
T2 = T2(:).';
n1 = length(T1);
n2 = length(T2);

if (n1>1) || (n2>1)
    T1 = repmat(T1,[1 n2]);
    T2 = repmat(T2,[n1 1]);
end

phi = acos((T1-T2)./(T1+T2))*180/pi;

end      % ssfp_angle.m
