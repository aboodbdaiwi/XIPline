function fil = lb_fun(n,bw,lbs,t0,plt)
%LB_FUN  Line broadening function (Lorentzian + Gaussian)
%   fil = lb_fun(n,bw,lbs,t0)
%     n   length
%    bw   bandwidth/spectral width        [Hz]
%   lbs   filter widths                   [Hz]
%         [Lorentzian Gaussian]
%    t0   shift filter by this time        [s]
%         'h' = filter reference in middle
%         default=0
%   plt   plotting (default=false)
% Literature: (1) NMR Data Processing; Hoch,Stern; 1996
% 6/2007 Rolf Schulte
% 2/2015 RFS: fixed bug in Gaussian linewidth: 2->4 in exp denominator
if nargin<3, help(mfilename); return; end

%% parse input
if ~exist('t0','var'), t0 = []; end
if isempty(t0),        t0 = 0; end
if ~exist('plt','var'), plt = []; end
if isempty(plt),        plt = false; end

t = (0:n-1)/bw;                                  % time
if ischar(t0), t0 = n/2/bw; end                  % filter reference on half
if (t0>t(end)), warning('t0 (=%g) > t(end) (%g)',t0,t(end)); end

fil = ones(1,n);
%% generate filter functions
for l=1:length(lbs)
    if lbs(l)~=0
        switch l
            case 1      % exponential filter; Ref 1 p 32
                fil = fil.*exp(-pi*lbs(l)*abs(t-t0));
            case 2      % Gaussian filter; Ref 1 p 33
                fil = fil.*exp(-(pi*lbs(l)*(t-t0)).^2/(4*log(2)));
            otherwise, error('unknown filter function: lbs too long');
        end
    end
end

%% plot filter
if plt, plot(t,fil); end

end      % main function lb_fun.m
