function filter_vec = fermi_filter_1D(window_size,Wf)
%   Inputs:
%      
%   Outputs:
%                   
%   Package: https://github.com/aboodbdaiwi/HP129Xe_Analysis_App
%
%   Author: Abdullah S. Bdaiwi
%   Work email: abdullah.bdaiwi@cchmc.org
%   Personal email: abdaiwi89@gmail.com
%   Website: https://www.cincinnatichildrens.org/research/divisions/c/cpir
%
%   Please add updates at the end. Ex: 3/10/24 - ASB: update .... 

%GE default for Wf is 10
% Centering corrected to window_size/2 on 10/16/17 by ZIC.

filter_vec = ones(window_size,1);
rf = window_size/2; %GE Defauly
v = 1:rf;
Hv = 1./(1+exp((v-rf)/Wf)); % make general for all window sizes 10/16/17, ZIC

filter_vec(1:window_size/2)=flip(Hv,2);
filter_vec(window_size/2+1:end)=  Hv;

end


