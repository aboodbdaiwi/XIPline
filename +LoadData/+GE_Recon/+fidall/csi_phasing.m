function d = csi_phasing(d,ind_phasing,pts4phasing)
%CSI_PHASING  Phase correct CSI voxel
%          dd = csi_phasing(d,ind_phasing,pts4phasing)
%           d   Spatially reconstructed, coil-combined data [ns,nx,ny,nz,nt]
%               ns=#spectral pts; nx+ny+nz=#spatial pts; nt=#time steps
% ind_phasing   Index list for FID/spectral points                     ([])
% pts4phasing   #largest signal used for phasing                        (7)
%
% 10/2024 Rolf Schulte
if nargin<1, help(mfilename); return; end



%% default input
if ~exist('ind_phasing','var'), ind_phasing = []; end
if ~exist('pts4phasing','var'), pts4phasing = []; end
if isempty(pts4phasing),        pts4phasing = 7; end


%% check for multi-coil data (coil dim=6)
nn = size(d);                % ns,nx,ny,nz,nt
if length(nn)>5
    warning('length(size(d))(=%g)>5; skipping phasing',length(nn));
    return;
end
nn = [nn , ones(1,5-length(nn))];
if nn(5)>1
    fprintf('nn(5)(=%g)>1; averaging along dim 5\n',nn(5));
    dd = mean(d,5);
else
    dd = d;
end


%% create index lists
% create index list of largest signals
% can be used for FID (starting at max) or spectrum (max in middle)
if isempty(ind_phasing)
    [~,ind] = sort(sqrt(sum(sum(sum(dd.*conj(dd),2),3),4)));
    if pts4phasing>nn(1)
        error('pts4phasing(=%g) > nn(1)(=%g)',pts4phasing,nn(1));
    end
    if pts4phasing<1, error('pts4phasing(=%g) < 1',pts4phasing); end
    ind_phasing = sort(ind((nn(1)-pts4phasing+1):nn(1)));
end


%% phasing of spectra/fid
pha = mean(dd(ind_phasing,:,:,:,:),1);
pha = pha./abs(pha);
% d = d.*repmat(conj(pha),[nn(1) 1 1 1 1]);
d = bsxfun(@times,d,conj(pha));


end      % csi_phasing.m
