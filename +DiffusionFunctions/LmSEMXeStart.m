function [D,ldd,Bi,Ci]=LmSEMXeStart(Do,diftime)
%   Inputs:
%      
%   Outputs:
%                   
%   Package: https://github.com/aboodbdaiwi/HP129Xe_Analysis_App
%
%   Please add updates at the end. Ex: 3/10/24 - ASB: update .... 

% calculates coefficients, diffusion times and diffusion lenghts before
% looping to speed up fit of image data (avoid unnecessary recalculations)

if nargin==1,
    diftime=3.5e-3;
end
if nargin<1,
    Do=0.14;
    diftime=3.5e-3;
end

D=0.001:0.001:Do;
ldd=1e4*sqrt(2*D*diftime);

betavect=[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9];
betai=0.1:0.01:0.9;
B=[0.145 0.197 0.243 0.285 0.382 0.306 0.360 0.435 0.7];
C=[0.89 0.50 0.35 0.25 0 0.13 0.22 0.4015 0.33];
Bi=interp1(betavect,B,betai,'cubic');
Ci=interp1(betavect,C,betai,'cubic');



