function [dd,smap] = csi_coil_combine(d,ind_ts,ind_svds,ind_phasing,nper,...
    thresh,nspat)
%CSI_COIL_COMBINE  SVDS Coil combination for CSI/MRS
%   Determines phases and weights of individual coil elements from 
%   FIDs or spectra, phases and weights channels, and then averages them
%
%  [dd,smap] = csi_coil_combine(d,ind_ts,ind_svds,ind_phasing,nper,thresh,nspat)
%          d   Spatially reconstructed data [ns,nx,ny,nz,nt,nc]
%              ns=#spectral pts; nx+ny+nz=#spatial pts; nt=#time steps; nc=#coils
%     ind_ts   Index list for averaging along time steps (dim=5)
%   ind_svds   Index list for FID/spectral points to consider in svds
%              default=128 largest signal
%ind_phasing   Index list for FID/spectral points to consider for phasing
%              if logical: true = 10 largest signal; false = skip phasing
%              default = true
%      nper    Percent of data classified as noise                     (15)
%              >0 -> perform noise decorrelation
%    thresh    Exclude noise coils                                      (2)
%     nspat    Extend by +-nspat pixels                                 (0)
%
%        dd    coil-combined FID/echo       [ns,nx,ny,nz,nt]
%      smap    sensitivity map              [nx,ny,nz,nc]
%
% 11/2024 Rolf Schulte
if nargin<1, help(mfilename); return; end

pts4svds = 128;              % #largest signal used for svds 
% (truncate data for numerical speed)
pts4phasing = 7;             % #largest signal used for phasing


%% default input
nn = size(d);                % ns,nx,ny,nz,nt,nc
if ~exist('ind_ts','var'),   ind_ts = []; end
if isempty(ind_ts),          ind_ts = (1:size(d,5)); end
if ~exist('ind_svds','var'), ind_svds = []; end
if ~exist('ind_phasing','var'), ind_phasing = []; end
if isempty(ind_phasing),     ind_phasing = true; end
do_phasing = false;
if isnumeric(ind_phasing), do_phasing = true; end
if islogical(ind_phasing), if ind_phasing, do_phasing = true; end; end
if ~exist('nper','var'),     nper = []; end
if isempty(nper),            nper = 15; end
if ~exist('thresh','var'),   thresh = []; end
if isempty(thresh),          thresh = 2; end
if ~exist('nspat','var'),    nspat = []; end
if isempty(nspat),           nspat = 0; end


%% check for multi-coil data (coil dim=6)
if length(nn)~=6
    warning('length(size(d))(=%g) ~= 6; skipping coil combination',length(nn));
    dd = d;
    smap = [];
    return
end


%% Noise decorrelation
if nper>0
    fprintf('Noise decorrelation\n');
    d = mri_noise_decorrelation(d,nper,6,thresh,0);
end
nc = size(d,6);


%% create index lists
if isempty(ind_svds) || do_phasing
    % create index list of largest signals for svds and phasing
    % can be used for FID (starting at max) or spectrum (max in middle)
    [~,ind] = sort(sqrt(sum(sum(sum(sum(sum(d.*conj(d),2),3),4),5),6)));
    if isempty(ind_svds)
        if nn(1)>pts4svds
            ind_svds    = sort(ind((nn(1)-pts4svds+1):nn(1)));
        else
            ind_svds = (1:nn(1));
        end
    end
    if do_phasing
        if nn(1)<pts4phasing, error('nn(1)(=%g) < %g',nn(1),pts4phasing); end
        ind_phasing = sort(ind((nn(1)-pts4phasing+1):nn(1)));
    end
end


if any(ind_svds<1),     warning('ind_svds<1');  ind_svds = ind_svds(ind_svds>0); end
if any(ind_svds>nn(1)), warning('ind_svds>ns'); ind_svds = ind_svds(ind_svds<=nn(1)); end

if nn(5)>1
    fprintf('Attention: nn(5)(=%g)>1; averaging along dim 5\n',nn(5));
    dorig = d;
    d = mean(d(:,:,:,:,ind_ts,:),5);
end


%% determine sensitivity map (combination weights)
smap = zeros([1 nn(2:4) 1 nc]);
for l2=1:nn(2)
    for l3=1:nn(3)
        for l4=1:nn(4)
            if nspat>0
                ii2 = l2+(-nspat:nspat); 
                ii2 = ii2((ii2>=1)&(ii2<=nn(2)));
                ii3 = l3+(-nspat:nspat);
                ii3 = ii3((ii3>=1)&(ii3<=nn(3)));
                ii4 = l2+(-nspat:nspat);
                ii4 = ii4((ii4>=1)&(ii4<=nn(4)));
                A = d(ind_svds,ii2,ii3,ii4,1,:);
                nA = size(A);
                A = reshape(A,prod(nA(1:5)),nA(6));
            else
                A = squeeze(d(ind_svds,l2,l3,l4,1,:));
            end
            [~,~,V] = svd(A,0);
            smap(1,l2,l3,l4,1,:) = conj(V(:,1));
        end
    end
end


%% normalisation
smap = bsxfun(@rdivide,smap,sum(abs(smap),6));


%% actual phasing and weighting plus coil combination
if nn(5)>1
    % dd = sum(dorig.*repmat(conj(smap),[nn(1) 1 1 1 nn(5) 1]),6);
    dd = sum(bsxfun(@times,dorig,conj(smap)),6);
else
    % dd = sum(d.*repmat(conj(smap),[nn(1) 1 1 1 1 1]),6);
    dd = sum(bsxfun(@times,d,conj(smap)),6);
end



%% phasing of spectra: set phase of first point of FID to zero
if do_phasing
    pha = mean(dd(ind_phasing,:,:,:,:),1);
    pha = pha./abs(pha);
    % dd = dd.*repmat(conj(pha),[nn(1) 1 1 1 1]);
    dd = bsxfun(@times,dd,conj(pha));
end


end      % csi_coil_combine.m
