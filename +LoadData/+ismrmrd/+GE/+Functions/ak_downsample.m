function d = ak_downsample(d,dec,bw,plt)
%AK_DOWNSAMPLE Downsampling and filtering of arbitrary k-space data 
%   Cartesian iFFT, filtering, discarding & FFT
%
% d = ak_downsample(d,dec,bw)
%   d  data (filter and decimate along 2nd dim)
% dec  downsample (take every dec'th point) dec=2,3,4... 
%      filter with bw<=1/dec first
%  bw  Filter sampling BW (fft -> bw -> ifft): 0<bw<1
%      (default = 0.95/dec)
%
% 10/2008 Rolf Schulte
% See also RAW2CG.
if (nargin<1), help(mfilename); return; end;
if ~exist('plt','var'), plt = true; end
si = size(d);

if dec<1, error('dec<1'); end
if ~exist('bw','var'), bw = []; end
if isempty(bw), bw = 0.95/dec; end
if ((bw>1)||(bw<0)), error('bw>1 or bw<0'); end

if length(si)<2, error('length(si)<2'); end
ns = si(2);

% lowpass filter data; reduce (noise) effective sampling BW
spec = fftshift(ifft(d,[],2),2);
xa = linspace(-1,1,ns);
if plt, figure; subplot(3,1,1); plot_complex(xa,spec(1,:,1)); end

xx = floor(bw*ns/2+0.5);
ind = floor(ns/2)+(-xx:xx-1)+1;
lbs = zeros(1,ns); lbs(ind)=1;
lbs = real(ifft2(ifftshift(gausswin(ns,ns/70).'.*fftshift(fft2(lbs)))));

si_rm = si; si_rm(2) = 1;
spec = spec.*repmat(lbs,si_rm);

lowbo = ceil(ns/2-ns/(2*dec)+1);
upbo  = ceil(ns/2+ns/(2*dec));
d = fft(ifftshift(spec(:,lowbo:upbo,:),2),[],2);
if length(size(d))~=length(si),
    fprintf('length(size(d)): orig=%g\tnew=%g\n',length(si),length(size(d)));
    si_rs = si; si_rs(2) = upbo-lowbo+1;
    d = reshape(d,si_rs);
    fprintf('\treshaping to %gD\n',length(size(d)));
end

if plt,
    subplot(3,1,2);
    plot(xa,lbs,'b',[1 1]/dec,[-0.1 1.1],'r',-[1 1]/dec,[-0.1 1.1],'r'); 
    axis([-1 1 -0.1 1.1])
    subplot(3,1,3); plot_complex(xa,spec(1,:,1));
end
