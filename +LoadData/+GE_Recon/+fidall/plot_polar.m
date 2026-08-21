function plot_polar(varargin)
% PLOT_PLOAR  Plots a complex function in polar coordinates.
%    PLOT_POLAR(w,z);
%
%    See also PLOT, PLOT_COMPLEX.
if nargin<1, help(mfilename); return; end
% error(nargchk(1,2,nargin));
switch nargin,
    case 1, 
        z{1} = varargin{1};
        w{1} = [1:length(z{1})];
        c{1} = '';
        n = 1;
    case 2,
        w{1} = varargin{1};
        z{1} = varargin{2};
        c{1} = '';
        n = 1;
    otherwise,
        if rem(nargin,3)~=0, error('nargin not multiple of 3'); end
        n = nargin/3;
        for l=1:n,
            w{l} = varargin{(l-1)*3+1};
            z{l} = varargin{(l-1)*3+2};
            c{l} = varargin{(l-1)*3+3};
        end
end
if isreal(z{1}), warning('Real input'); end

% clf reset;
subplot(2,1,1);
plot(w{1},abs(z{1}),c{1});
for l=2:n,
    hold on; plot(w{l},abs(z{l}),c{l}); hold off;
end
ylabel('abs');
axis tight;
subplot(2,1,2);
plot(w{1},unwrap(angle(z{1})),c{1});
for l=2:n,
    hold on; plot(w{l},unwrap(angle(z{l})),c{l}); hold off;
end
ylabel('phase');
axis tight;
