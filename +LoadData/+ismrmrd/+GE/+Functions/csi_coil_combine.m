function [dd,smap] = csi_coil_combine(d,ind_ts,ind_svds,ind_phasing)
%CSI_COIL_COMBINE  SVDS Coil combination for CSI/MRS
%   Determines phases and weights of individual coil elements from 
%   FIDs or spectra, phases and weights channels, and then averages them
%
% [dd,smap] = csi_coil_combine(d,method,ind_fid,ind_ts)
%          d   Spatially reconstructed data [ns,nx,ny,nz,nt,nc]
%              ns=#spectral pts; nx+ny+nz=#spatial pts; nt=#time steps; nc=#coils
%     ind_ts   Index list for averaging along time steps (dim=5)
%   ind_svds   Index list for FID/spectral points to consider in svds
%              default=256 largest signal
%ind_phasing   Index list for FID/spectral points to consider for phasing
%              if logical: true = 10 largest signal; false = skip phasing
%              default = true
%
%        dd  coil-combined FID/echo       [ns,nx,ny,nz,nt]
%      smap  sensitivity map              [nx,ny,nz,nc]
%
% 10/2019 Rolf Schulte
if nargin<1, help(mfilename); return; end

pts4svds = 256;              % #largest signal used for svds 
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
if ~strcmpi(class(d),'double')
    warning('class(d)(=%s) not double; converting to double for sdvs',...
        class(d));
    d = double(d);
end

%% check for multi-coil data (coil dim=6)
if length(nn)~=6
    warning('length(size(d))(=%g) ~= 6; skipping coil combination',length(nn));
    dd = d;
    smap = [];
    return
end

%% create index lists
if isempty(ind_svds) || do_phasing
    % create index list of largest signals for svds and phasing
    % can be used for FID (starting at max) or spectrum (max in middle)
    [~,ind] = sort(sqrt(sum(sum(sum(sum(sum(d.*conj(d),2),3),4),5),6)));
    if nn(1)>pts4svds
        ind_svds    = sort(ind((nn(1)-pts4svds+1):nn(1)));
    else
        ind_svds = (1:nn(1));
    end
    if nn(1)<pts4phasing, error('nn(1)(=%g) < %g',nn(1),pts4phasing); end
    ind_phasing = sort(ind((nn(1)-pts4phasing+1):nn(1)));
end


if any(ind_svds<1),     warning('ind_svds<1');  ind_svds = ind_svds(ind_svds>0); end
if any(ind_svds>nn(1)), warning('ind_svds>ns'); ind_svds = ind_svds(ind_svds<=nn(1)); end

if nn(5)>1
    fprintf('Attention: nn(5)(=%g)>1; averaging along dim 5\n',nn(5));
    dorig = d;
    d = mean(d(:,:,:,:,ind_ts,:),5);
end


%% determine sensitivity map (combination weights)
smap = zeros([1 nn(2:4) 1 nn(6)]);
for l2=1:nn(2)
    for l3=1:nn(3)
        for l4=1:nn(4)
            A = squeeze(d(ind_svds,l2,l3,l4,1,:));
            if ~isempty(regexpi(version,'R2015b'))
                [U,S,V,flag] = svds(A,2);
                V = V(:,1);   % fix bug in svds function in Matlab R2015b
            else
                [U,S,V,flag] = svds(A,1);
            end
            smap(1,l2,l3,l4,1,:) = conj(V);
            if flag~=0
                fprintf('l2=%g; l3=%g; l4=%g; flag=%g\n',l2,l3,l4,flag);
            end
        end
    end
end
% normalisation
smap = smap./repmat(sum(abs(smap),6),[1 1 1 1 1 nn(6)]);

% tmp=zeros(8,8,8,2);tmp(:,:,:,1)=squeeze(abs(smap));tmp(:,:,:,2)=squeeze(angle(smap));
% imagesc_row(tmp,[],'row')

%% actual phasing and weighting plus coil combination
if nn(5)>1
    dd = sum(dorig.*repmat(conj(smap),[nn(1) 1 1 1 nn(5) 1]),6);
else
    dd = sum(d.*repmat(conj(smap),[nn(1) 1 1 1 1 1]),6);
end

%% phasing of spectra: set phase of first point of FID to zero
if do_phasing
    %pha = angle(mean(dd(ind_phasing,:,:,:,:),1))+pi/4;
    % pha = mean(angle(dd(ind_phasing,:,:,:,:)),1);
    %dd = dd.*repmat(exp(-1i*pha),[nn(1) 1 1 1 1]);
    
    pha = mean(dd(ind_phasing,:,:,:,:),1);
    pha = pha./abs(pha);
    dd = dd.*repmat(conj(pha),[nn(1) 1 1 1 1]);
    
end

end      % main function csi_coil_combine.m
