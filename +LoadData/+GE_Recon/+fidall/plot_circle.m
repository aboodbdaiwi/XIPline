function plot_circle(x,y,r,clr)
%PLOT_CIRCLE  Overlay circle on imagesc
% plot_circle(x,y,r,clr)
%    x  x coordinates
%    x  x coordinates
%    r  radius
%  clr  colour
%
% 10/2017 Rolf Schulte (adapted from mathworks homepage)
if (nargin<1), help(mfilename); return; end;

hold on
th = 0:pi/50:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
plot(xunit,yunit,clr);
hold off

end   % main function plot_circle.m
