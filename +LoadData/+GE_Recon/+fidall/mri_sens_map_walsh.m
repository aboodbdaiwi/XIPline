function smap = mri_sens_map_walsh(b,filtsize,verb)
%MRI_COIL_COMBINE_WALSH  MRI coil combination 
%[b,smap] = mri_sens_map_walsh(b,filtsize,verb)
%       b   Reconstructed 5D data (nx,ny,nz,nt,nc)
%filtsize   Filter size for spatial noise mask                          (4)
%    verb   Verbose                                                  (true)
%       b   Coil-combined images
%    smap   Sensitivity map
%
% Literature: mrm43_682_2000_walsh.pdf
% 12/2022 Rolf Schulte
if nargin<1, help(mfilename); return; end


%% input parameters
if ~exist('verb','var'), verb = []; end
if isempty(verb),        verb = true; end
if ~exist('filtsize','var'), filtsize = []; end
if isempty(filtsize),    filtsize = 4; end
if verb>0, timerVal = tic; end         % record execution time


%% image dimensions
if length(size(b))~=5, warning('data(dim=%g) ~= 5D',length(size(b))); end
[nx,ny,nz,nt,nc] = size(b);
nn = nx*ny*nz;
if verb, fprintf('Calculating sensitivity maps via Walsh method: #coils=%g\n',nc); end
if nc==1, warning('#coils==1'); end
% if nt~=1, warning('#timesteps(=%d)~=1',nt); end




%% correlation matrix
if verb, fprintf('Calculating correlation matrix\n'); end
Rs = complex(zeros(nc,nc,nx,ny,nz,nt,class(b)));
onz = ones(filtsize,filtsize,filtsize,class(b));
for lc1=1:nc
    for lc2=1:nc
        if verb>1, fprintf('.'); end
        % Rs(lc1,lc2,:,:,:) = permute(convn(b(:,:,:,1,lc1).*conj(b(:,:,:,1,lc2)),...
        %     onz,'same'),[4 5 1 2 3]);
        for lt=1:nt
            Rs(lc1,lc2,:,:,:,lt) = permute(convn(b(:,:,:,lt,lc1).*conj(b(:,:,:,lt,lc2)),...
                onz,'same'),[4 5 1 2 3]);
        end
    end
end
if verb>1, fprintf('\n'); end
if nt~=1
    if verb>0
        fprintf('nt(=%d) > 1: averaging correlation matrix\n',nt);
    end
    Rs = mean(Rs,6);
end


%% compute and apply filter at each voxel
if verb, fprintf('Calculating coil sensitivity maps\n'); end
Rs   = Rs(:,:,:);
smap = complex(zeros(nc,nn,class(b)));
for ll = 1:nn
    [U,~] = svd(Rs(:,:,ll),0);
    smap(:,ll) = U(:,1);
end
clear Rs
if any(any(isnan(smap)))
    warning('isnan(smap)): zeroing');
    smap(isnan(smap)) = 0;
end


%% normalise: sum-over-coils=1
smap = bsxfun(@rdivide,smap,sqrt(sum(smap.*conj(smap),1)));


%% reshaping to image size
smap = reshape(smap',[nx,ny,nz,1,nc]);


%% print info
if verb>0
    fprintf('mri_sens_map_walsh: runtime = %.2f [s]\n',toc(timerVal));
end

end      % mri_sens_map_walsh.m
