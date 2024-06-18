function [kPX,rXeff,rat] = memo(b,tm,ind_met,ind_ts,ts0,plt_tdfd)
%MEMO  Metabolic two-side exchange rate modelling in frequency domain
%
% kPX = memo(b,tm,ind_met,imd_ts,ts0,plt_tdfd)
%        b  Metabolic images  (mtx,mtx,#exc,#metabolites,#slice)
%       tm  Metabolite repetition time
%  ind_met  Index to metabolite [P X]
%   ind_ts  Index to time steps
%      ts0  Set first time step to zero
% plt_tdfd  Plot time steps and frequency response of signal at max(X)
%
%      kPX  Metabolic exchange rate images   (mtx mtx #slices)
%    rXeff  Relaxation rate images           (mtx mtx #slices)
%      rat  Ratio X/Pyr (summed => @ f=0)    (mtx mtx #slices)
% 4/2013 Rolf Schulte
if (nargin<1), help(mfilename); return; end;
if ~exist('ind_ts','var'), ind_ts = []; end
if ~exist('ts0','var'), ts0 = false; end
if ~exist('plt_tdfd','var'), plt_tdfd = false; end
if ~isreal(b),
    fprintf('Attention: b complex -> taking abs\n');
    b = abs(b);
end
if isodd(length(ind_ts)), error('length(ind_ts) must be even'); end
tmp = abs(b(:,:,end,:,:)); 
noise = std(tmp(:));                   % noise determination
lambda = noise.^2*3d3;                 % regularisation for kPX inversion
lambda_rat = noise*5d2;
if ~isempty(ind_ts), b = b(:,:,ind_ts,:,:); end; % truncating time steps
if ts0, b(:,:,1,:,:) = 0; end                    % set first time step to 0
[n1,n2,nexc,nmet,nslices] = size(b);
kPX = zeros(n1,n2,nslices);            % turnover maps
rXeff = zeros(n1,n2,nslices);          % relaxation maps
rat = zeros(n1,n2,nslices);            % ratio maps
denom = zeros(n1,n2,nslices);          % denominator

nf = nexc;         % frequency discretisation (zero filling unecessary)

bFD = fftshift(fft(b,nf,3),3);         % FFT along time steps
ff = (-nf/2:nf/2-1)/tm;                % frequency grid
if (ff(nf/2+1)~=0), warning('memo:ff0','(ff(nf/2+1)~=0)'); end
lf0 = nf/2+1; lf1 = lf0+1;

fprintf('noise=%g lambda=%g lambda_rat=%g\n',noise,lambda,lambda_rat);

for lls=1:nslices,                     % loop through slices
    p = bFD(:,:,:,ind_met(1),lls);     % Pyr in frequency domain
    x = bFD(:,:,:,ind_met(2),lls);     % X in frequency domain
    %     p = p.*exp(-1i*repmat(angle(p(:,:,lf0)),[1 1 nf]));
    %     x = x.*exp(-1i*repmat(angle(x(:,:,lf0)),[1 1 nf]));
    % abs(b) seems to yield better results than phasing metabolic images
    
    if max(x(:))>max(p(:)), warning('memo:max','max(x(:))>max(p(:))'); end
    den = p(:,:,lf1).*x(:,:,lf0) - x(:,:,lf1).*p(:,:,lf0);
    denom(:,:,lls) = den;
    den = den + lambda;
    k = (1i*2*pi*ff(lf1)/nf) .* x(:,:,lf1).*x(:,:,lf0)./den;
    % calculate actual turnover
    %   factor (1i*2*pi*ff(lf1)/nf) -> derivative of numerical fft
    kPX(:,:,lls) = k;
    fprintf('kPX real=%.3g imag=%.3g abs=%.3g\n',...
        sqrt(sum(real(k(:)).^2)),sqrt(sum(imag(k(:)).^2)),...
        sqrt(sum(abs(k(:)).^2)));
    rXeff(:,:,lls) = (1i*2*pi*ff(lf1)/nf) .* x(:,:,lf1).*p(:,:,lf0)./den;

    tmp = [p(1,:,lf0) p(end,:,lf0) p(:,1,lf0).' p(:,end,lf0).' ...
        p(2,:,lf0) p(end-1,:,lf0) p(:,2,lf0).' p(:,end-1,lf0).'];
    rat(:,:,lls) = x(:,:,lf0)./(p(:,:,lf0)+lambda_rat);
    
%     figure(lls); clf;
%     subplot(3,2,1); imagesc(abs(p(:,:,lf0)),[0 0.9*max(abs(p(:)))]);
%     title('p(0)'); colorbar; axis image off;
%     subplot(3,2,2); imagesc(abs(x(:,:,lf0)),[0 0.9*max(abs(x(:)))]);
%     title('x(0)'); colorbar; axis image off;    
%     subplot(3,2,3); imagesc(abs(den),[0 0.9*max(abs(den(:)))]);
%     title('(p(1)\cdotx(0)-x(1)\cdotp(0))+\lambda'); colorbar; axis image off
%     subplot(3,2,4); imagesc(abs(kPX(:,:,lls)),[0 0.1]);
%     title('kPX'); colorbar; axis image off
%     subplot(3,2,5); imagesc(abs(rXeff(:,:,lls)),[0 0.1]);
%     title('|rXeff|'); colorbar; axis image off
%     subplot(3,2,6); imagesc(rat,[0 0.5]);
%     % [0 10*mean(abs(rat(:)))]);
%     title('x(0)/p(0)'); colorbar; axis image off
end

% plot frequency response of signal at max(X)
if plt_tdfd,
    [tmp,l1] = max(b(:,:,4,1,1),[],1);
    [tmp,l2] = max(tmp,[],2);
    l1 = l1(l2);
    fprintf('l1=%g l2=%g\n',l1,l2);

    figure(100); clf;
    subplot(4,1,1); plot(squeeze(b(l1,l2,:,:,lls)))
    subplot(4,1,2); plot(ff,squeeze(abs(bFD(l1,l2,:,:,lls))))
    subplot(4,1,3); plot(ff,squeeze(real(bFD(l1,l2,:,:,lls))))
    subplot(4,1,4); plot(ff,squeeze(imag(bFD(l1,l2,:,:,lls))))
    
    xx = reshape(x(l1,l2,:),[1 nf]);
    pp = reshape(p(l1,l2,:),[1 nf]);
    kPX = 1i*2*pi*ff .* xx.*x(l1,l2,(nf/2+1))./...
        (pp.*x(l1,l2,(nf/2+1)) - xx.*p(l1,l2,(nf/2+1)));
    rXeff = 1i*2*pi*ff .* xx.*p(l1,l2,(nf/2+1))./...
        (pp.*x(l1,l2,(nf/2+1)) - xx.*p(l1,l2,(nf/2+1)));
    figure(101); clf; plot_polar(ff,rXeff,'r',ff,kPX,'b');
    subplot(2,1,1); axis([-0.1 0.1 0 0.05]);
end

% plotting
tmp = zeros(n1,n2,nslices,4);
tmp(:,:,:,1) = squeeze(bFD(:,:,lf0,ind_met(1),:));
tmp(:,:,:,2) = squeeze(bFD(:,:,lf0,ind_met(2),:));
tmp(:,:,:,3) = rat;
tmp(:,:,:,4) = abs(denom);
figure(1); clf; imagesc_row(tmp,[],'row');
xlabel('slices');
ylabel('denom    ratio      x(0)      p(0)');

tmp = zeros(n1,n2,nslices,4);
tmp(:,:,:,1) = real(kPX);   tmp(:,:,:,2) = imag(kPX);
tmp(:,:,:,3) = real(rXeff); tmp(:,:,:,4) = imag(rXeff);
figure(2); clf; imagesc_row(tmp,[-0.1 0.1]); colorbar;
xlabel('slices');
ylabel('imag(rXeff)    real(rXeff)    imag(kPX)    real(kPX)');

