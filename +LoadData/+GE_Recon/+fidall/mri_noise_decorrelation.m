function b = mri_noise_decorrelation(b,nper,dim_coil,thresh,verb)
%MRI_NOISE_DECORRELATION  Noise pre-whitening for multi-channel Rx data
%       b = mri_noise_decorrelation(b,nper,dim_coil,verb)
%       b   Reconstructed multi-coil data
%    nper   Percent of data classified as noise                        (15)
%dim_coil   Coil dimension                                              (5)
%  thresh   Threshold to exclude noisy coils (diag(psi)>thresh; 0=off)  (0)
%    verb   Verbose (>0) + plotting (>1 &>2)                            (1)
%
%  1/2024 Rolf Schulte
if nargin<1, help(mfilename); return; end


%% input + default parameters
if ~exist('nper','var'),     nper = []; end
if isempty(nper),            nper = 15; end
if nper<0,   warning('nper(=%g)<0%%',nper); end
if nper>100, warning('nper(=%g)>100%%',nper); end
if ~exist('dim_coil','var'), dim_coil = []; end
if isempty(dim_coil),        dim_coil = 5; end
if ~exist('thresh','var'),   thresh = []; end
if isempty(thresh),          thresh = 0; end
if ~exist('verb','var'),     verb = []; end
if isempty(verb),            verb = 1; end


%% checks + sizes
nn = size(b);
ndim = length(nn);
if ndim>dim_coil
    warning('length(nn)(=%g)>dim_coil(=%g)',ndim,dim_coil);
end
nc = size(b,dim_coil);


%% noise decorrelation/pre-whitening
if ((nper>0) && (nc>1))
    % determine noise covariance psi
    if verb>0, fprintf('Calculating noise covariance matrix psi\n'); end
    [psi,ind_noise] = sub_noise_covariance(b,nper,dim_coil,verb,90);
    if any(isnan(psi(:)))
        warning('isnan(psi) -> skipping noise decorrelation');
        return
    end
    
    % pre-whitening
    % literature: ISMRM2014, sunrise, Michael S Hansen; Nuts & Bolts...
    % http://hansenms.github.io/sunrise/sunrise2014/
    if verb>0, fprintf('Noise pre-whitening/decorrelation\n'); end
    try
        L = chol(psi,'lower');
    catch
        err = lasterror;
        warning('Cholesky factorization failed; skipping noise decorrelation');
        fprintf('%s\n',err.message);
        return
    end
    ii = (1:length(nn));
    ii = ii(ii~=dim_coil);
    b = L\reshape(permute(b,[dim_coil ii]),[nc prod(nn(ii))]);
    % b = bsxfun(@rdivide,b,diag(psi));
    if ((thresh>0) && any(diag(psi)>thresh))
        iex = diag(psi)>thresh;
        warning('Excluding %d/%d coil elements with noise covariance > %g',...
            sum(iex),nc,thresh);
        fprintf('Excluded coils = '); disp(find(iex).');
        b = b(~iex,:);
        nc = nc - sum(iex);
    end
    % b = pinv(L)*reshape(permute(b,[dim_coil ii]),[nc prod(nn(ii))]);
    b = permute(reshape(b,[nc,nn(ii)]),[2:ndim,1]);
    
    
    % determine noise covariance after decorrelation
    if verb>1
        fprintf('Calculating noise covariance matrix after decorrelation\n');
        % sub_noise_covariance(b,nper,dim_coil,verb,91);
        sub_noise_covariance(b,ind_noise,dim_coil,verb,91);
    end
else
    fprintf('Attention: nper(=%g)<=0 or nc(=%g)==1\n',nper,nc);
    fprintf(' -> skipping noise decorrelation\n');
end

end      % mri_noise_decorrelation.m


%% sub-functions
function [psi,ind_noise] = sub_noise_covariance(br,noise_threshold,dim_coil,verb,figid)

%% reshape br to [nc,nrest]
nn = size(br);
ii = (1:length(nn));
ii = ii(ii~=dim_coil);

br = permute(br,[dim_coil ii]);
nc = size(br,1);
br = reshape(br,[nc prod(nn(ii))]);


%% calculate standard-deviation and threshold data to extract noise
if length(noise_threshold)==1
    brms = sqrt(mean(br.*conj(br),1));
    istart = sum(abs(brms)<1d-10)+1;               % exclude zeroes
    iend = ceil(noise_threshold/100*size(br,2)) + istart;
    if iend<1, error('iend(=%g)<1',iend); end
    if iend>size(br,2), error('iend(=%g)>size(br,2)(=%g)',iend,size(br,2)); end

    [~,ind_noise] = sort(brms);
    ind_noise = ind_noise(1,istart:iend);
else
    ind_noise = noise_threshold;
end

noise = br(:,ind_noise);
if verb>0
    fprintf('#noise samples: %.3g%%\n',size(noise,2)/size(br,2)*100);
end


%% calculate noise covariance and scale to mean of diagonal element=1
psi = noise*noise';
psi = psi/mean(diag(abs(psi)));   % normalise to average


%% plotting
if ((verb>1) && (figid>0))
    figure(figid); clf
    imagesc(abs(psi),[0 2]); 
    colorbar;
    title('Noise covariance');
    drawnow
end

end      % sub_coise_covariance
