function ga = golden_angle(what,do_deg)
%GOLDEN_ANGLE  Calculate golden angle for pseudo random rotations
%    ga = golden_angle(what,do_deg)
%  what   0=Wikipedia,1=normal golden angle, >1->tiny golden angle
%do_deg   Convert output angle ga from [rad] to [deg]
%
% Literature: https://en.wikipedia.org/wiki/Golden_angle
%             itomi34_1262_2015_wundrak.pdf
% 12/2020 Rolf Schulte
if (nargin<1), help(mfilename); return; end


%% input parameters
if ~exist('do_deg','var'), do_deg = []; end
if isempty(do_deg), do_deg = false; end
if abs(what-floor(what))
    warning('what(=%g) not a natural number; rounding to %g\n',what,round(what));
    what=round(what);
end


%% calculate different golden angles
if what<0, error('what(=%d)<0',what); end
if what==0
    % https://en.wikipedia.org/wiki/Golden_angle
    ga = pi*(3-sqrt(5));
else
    % itomi34_1262_2015_wundrak.pdf
    % ga = pi/((sqrt(5)+1)/2);  
    % 1=normal golden angle
    % what>1: tiny golden angle
    ga = pi./((sqrt(5)+2*what-1)/2);
end


%% convert from [rad] to [deg]
if do_deg, ga = ga*180/pi; end


end      % golden_angle.m
