function [snr,signal,noise,f_max,ind_noise] = mrs_snr(fid,spec,bw,signal_bound,ind_noise,verb)
%MRS_SNR  Calculate SNR of MRS
%         noise=std of complex spectrum outside spectral ROI
%         signal = max of 3 methods: 
%                  1=max(abs(spec))
%                  2=zero-crossing of imag(spec) near max
%                  3=slow FT with high discretisation around max
%
%[snr,signal,noise,f_max,ind_noise] = mrs_snr(fid,spec,bw,signal_bound,ind_noise,verb)
%         fid  time-domain data (2D; sampling along 2)
%        spec  freq-domain data (optional; default=[])
%              ifempty: spec=ifftshift(fft(conj(fid),[],2),2);
%          bw  sampling bandwidth                     [Hz]
%signal_bound  take max below (default=[]=all)        [Hz]
%   ind_noise  index list for noise calculation; default
%              =(((hz<-1050)&(hz>hz(1)+20))|((hz>600)&((hz<hz(end)-20))));
%        verb  verbose mode (default=false)
%         snr  signal-to-noise ratio
%      signal  maximum signal
%       noise  noise level
%       f_max  frequency at maximum signal            [Hz]
%
% 3/2015 Rolf Schulte
if nargin<1, help(mfilename); return; end

%% default parameters
if ~exist('bw','var'), error('input BW'); end
if ~exist('signal_bound','var'), signal_bound = []; end
if isempty(signal_bound), signal_bound = bw; end
if ~exist('ind_noise','var'), ind_noise = []; end
if ~exist('verb','var'),     verb = []; end
if isempty(verb),            verb = false; end


%% input checks
if length(size(fid))~=2, warning('fid not 2D'); end
[nf,ns] = size(fid);
if nf~=1, warning('size(fid,1) (=%g) ~= 1',nf); end
if isempty(spec),spec = ifftshift(fft(conj(fid),[],2),2); end
if length(size(spec))~=2, warning('spec not 2D'); end
[nf2,zf] = size(spec);
if (nf~=nf2), warning('nf(=%g) ~= nf2(=%g)',nf,nf2); end
if length(bw)~=1, warning('length(bw)(=%g) ~= 1',length(bw)); end


%% axes + pars
hz  = -(-zf/2:zf/2-1)*bw/zf;            % frequency axis
time = (0:ns-1)/bw;                    % time axis


%% calculate noise
if isempty(ind_noise),
    ind_noise = (((hz<-1050)&(hz>min(hz)+20))|((hz>600)&((hz<max(hz)-20))));
end
noise = std((spec(:,ind_noise)),[],2);
%ind_noiseTD = (ns-400:ns-10);
%noise = std((ddd(1,ind_noiseTD)))*sqrt(zf)


%% calculate signal
f_max  = zeros(nf,3); signal  = zeros(nf,3);

% method 1: maximum signal 
[signal(:,1),sigind1] = max(abs(spec).*repmat((hz<signal_bound),[nf 1]),[],2);
f_max(:,1) = hz(sigind1).';

% method 2: signal at frequency of zero crossing of imaginary part
for l1=1:nf
    f_max(l1,2) = interp1(imag(spec(l1,((1:-1:-1)+sigind1(l1)))),...
        hz((-1:1)+sigind1(l1)),0,'linear');
    if ~isnan(f_max(l1,2)),
        signal(l1,2) = abs(fid(l1,:)*(exp(1i*2*pi*time.'*f_max(l1,2))));
    end

    % method 3: signal around maximum
    ff = f_max(l1,1) + linspace(-bw/zf,bw/zf,64);
    ss = abs(fid(l1,:)*(exp(1i*2*pi*time.'*ff)));     % slow FT

    [signal(l1,3),tmp] = max(ss);
    f_max(l1,3) = ff(tmp);
    if verb,
        fprintf('l1=%g: signal = (%g, %g, %g)\n',l1,signal(l1,:));
        fprintf('        f_max = (%g, %g, %g)\n',f_max(l1,:));
    end
end

[signal,tmp] = max(signal,[],2);
f_max = f_max(:,tmp);
if verb,
    for l1=1:nf,
        fprintf('signal = %g; f_max = %g\n',signal(l1),f_max(l1));
    end
end
snr = signal./noise;
