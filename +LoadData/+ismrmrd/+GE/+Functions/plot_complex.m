function plot_complex(w,z)
% PLOT_COMPLEX  Plots a complex function in Cartesian coordinates.
%    PLOT_COMPLEX(w,z);
%
%    abs = red
%    Re  = blue
%    Im  = green dashed
%
%    See also PLOT, PLOT_POLAR.

if nargin<1, help(mfilename); return; end

if nargin<2,
  z = w; 
  w = (1:numel(z));
end
if isreal(z), warning('Real input'); end

% clf reset;
plot(w,abs(z),'r-', w,imag(z),'g--', w,real(z),'b');
axis tight;

title('abs = r   re = b   im = g--');

