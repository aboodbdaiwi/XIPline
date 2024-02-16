function [Lm,Hk]=LmSEMXe2(DDC,alpha,D,ldd,Bi,Ci)
%   Inputs:
%      
%   Outputs:
%                   
%   Package: https://github.com/aboodbdaiwi/HP129Xe_Analysis_App
%
%   Please add updates at the end. Ex: 3/10/24 - ASB: update .... 

% Calculates Lm and diffusion lenght distribution Hk for a pixel with
% parameters DDC and alpha.
% Use after running LmSEMXeStart outside loops for speed

if alpha>0.9
    alpha=0.9;
    %disp('alpha is larger than 0.9, estimating Lm for alpha=0.9 instead')
elseif alpha<0.5
    alpha=0.5;
    %disp('alpha is less than 0.1, estimating Lm for alpha=0.1 instead')
end


t0=1./DDC;
delta=alpha.*(alpha-0.5)./(1-alpha);
fk=1+Ci(round((alpha-0.1).*100)+1).*(D.*t0).^delta;
Hk=(t0.*Bi(round((alpha-0.1).*100+1))./(D*t0).^((1-alpha/2)/(1-alpha))).*exp(-((1-alpha).*alpha.^(alpha./(1-alpha)))./(D.*t0).^(alpha./(1-alpha))).*fk;
Hk=Hk/sum(Hk);
Lm=sum(Hk.*ldd);


