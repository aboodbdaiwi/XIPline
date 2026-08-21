function bb = mrsi2pseuimage(spec,zf,npts,t_comb,how2ind)
%MRSI2PSEUIMAGE  Calculate pseudo-image from MRSI data
%     bb = mrsi2pseuimage(spec,zf,npts,t_comb,how2ind)
%   spec   MRSI [#s,#x,#y,#z,#t,#c]
%     zf   Zero-filling  [#x,#y,#z]              ([128,128,size(spec,4)])
%   npts   #largest spectral points for rms-averaging                   (4)
% t_comb   Combine time (0=off,1=mean,2=rms)                            (1)
%how2ind   How to do indexing for rms spectral selection                (1)
%          0=skip
%          1=npts-largest individual voxel
%          2=npts-largest global voxel
%          3=central spectral point
%
%    bb   Pseudo images [#x,#y,#z,#t]
%
%  6/2024 Rolf Schulte
if nargin<1, help(mfilename); return; end


%% misc input
if ~exist('zf','var'),       zf = []; end
if isempty(zf),              zf = [128,128]; end
if length(zf)<2,             zf = [zf,zf]; end
if length(zf)<3, zf = [zf(:).',size(spec,4)]; end
if length(zf)<4, zf = [zf(:).',size(spec,5)]; end
for ll=1:4
    if zf(ll)<size(spec,ll+1), zf(ll) = size(spec,ll+1); end
end
% if zf(1)<size(spec,2), zf(1) = size(spec,2); end
% if zf(2)<size(spec,3), zf(2) = size(spec,3); end
if ~exist('npts','var'),     npts = []; end
if isempty(npts),            npts = 4; end
if ~exist('t_comb','var'),   t_comb = []; end
if isempty(t_comb),          t_comb = 1; end
if ~exist('how2ind','var'),  how2ind = []; end
if isempty(how2ind),         how2ind = 1; end


%% combine time
if size(spec,5)>1
    switch t_comb
        case 0
        case 1, spec = mean(spec,5); zf(4) = 1;
        case 2, spec = sqrt(mean(spec.*conj(spec),5)); zf(4) = 1;
        otherwise, warning('t_comb (=%g) undefined',t_comb);
    end
end


%% calculate RMS signal intensities
% spec = ifft(fftshift(spec,1),[],1);
bb = abs(spec.*conj(spec));                 % square
if size(bb,6)>1, bb = mean(bb,6); end       % coil combination 
switch how2ind
    case 0
        fprintf('Skipping pseudo image generation\n');
        bb = [];
        return
    case 1    % rms of largest individual voxel
        bb = sort(bb,'descend');
        bb = sqrt(mean(bb(1:npts,:,:,:,:),1));
        fprintf('RMS of %d-largest individual voxel\n',npts);
    case 2    % rms of same spectral index points
        bs = mean(mean(mean(mean(bb,2),3),4),5);
        [~,ii] = sort(bs,'descend');
        bb = sqrt(mean(bb(ii(1:npts),:,:,:,:),1));
        fprintf('RMS of %d-largest global spectral points: ind=%s\n',...
            npts,num2str(ii(1:npts)));
    case 3    % rms of central spectral point
        ii = ceil(size(bb,1)/2+0.5);
        bb = sqrt(mean(bb(ii,:,:,:,:),1));
        fprintf('RMS of central spectral point: ind=%d\n',ii);
    otherwise
        error('how2ind(=%g)~=1,2 or 3',how2ind);
end
bb = permute(bb,[2 3 4 5 1]);               % remove spectral dimension


%% spatial interpolation
if any(zf>0)
    bb = truma(bb,true,zf);
    bb(bb<0) = 0;
end

end      % mrsi2rmsimage.m
