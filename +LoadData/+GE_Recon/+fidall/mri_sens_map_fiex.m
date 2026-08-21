function smap = mri_sens_map_fiex(bc,bref,zf,dim,do_fiex,mtx,mask,verb)
%MRI_SENS_MAP  Extract 2D or 3D receive coil sensitivity maps 
%     smap = mri_sens_map_fiex(bc,bref,zf,dim,do_fiex,mtx,mask,verb)
%       bc   Surface coil image        [mtx(1),mtx(2),mtx(3),nt,nc]
%     bref   Body coil image (opt; otherwise root-mean-square of bc)
%       zf   Desired final smap matrix size                            ([])
%      dim   Data dimension: 2 or 3
%  do_fiex   Piece-wise polynomial smoothing and extrapolation       (true)
%            3D data can be computationally expensive
%      mtx   Resample to matrix size for calculations                  (32)
%     mask   Calculate over masked area                             (false)
%     verb   Verbose (>0) + plotting (>1 &>2)                           (1)
%
%     smap   sensitivity maps          [zf(1),zf(2),zf(3),1,nc]
%
%  3/2023 Rolf Schulte
if (nargin<1), help(mfilename); return; end


%% input parameters
if ~exist('bref','var'),       bref = []; end
if ~exist('zf','var'),         zf = []; end
if ~exist('dim','var'),        dim = []; end
if ~exist('do_fiex','var'),    do_fiex = []; end
if isempty(do_fiex),           do_fiex = true; end
if ~exist('mtx','var'),        mtx = []; end
if ~exist('mask','var'),       mask = []; end
if ~exist('verb','var'),       verb = []; end
if isempty(verb),              verb = true; end


%% checks
if verb>0
    fprintf('Calculating Rx coil sensitivities: polynomial smoothing\n');
end
nn = size(bc);
if  length(nn)~=5, error('length(size(bc))(=%g)~=5', length(nn)); end
nt = nn(4);                            % #time steps
nc = nn(5);                            % #coils
if nc==1, error('nc(=%g)==1',nc); end
if isempty(dim)
    if nn(3)>5
        dim = 3;
    else
        dim = 2;
    end
    warning('dim empty assuming %D data',dim);
end
if ~any(dim==[2 3]), error('dim(=%g)~=2 or 3',dim); end
if isempty(mtx),  mtx = 32; end
if length(mtx)<2, mtx = [mtx mtx(1)]; end
if length(mtx)<3
    if dim==2
        mtx = [mtx nn(3)];
    else
        mtx = [mtx mtx(1)];
    end
end
for ll=1:3
    if mtx(ll)>nn(ll)
        warning('mtx(%d)(=%d)>nn(%d)(=%d): capping',ll,mtx(ll),ll,nn(ll));
        mtx(ll) = nn(ll);
    end
end
if ~isempty(zf)
    if length(zf)==dim
        warning('length(zf)(=%g)==dim(=%g)',length(zf),dim);
    end
else
    zf = nn(1:dim);
end
if dim==2, zf(3) = nn(3); end
if isreal(bc), warning('bc is real'); end
if nt>1
    warning('nt(=%d)>1: summing data along time');
    bc = mean(bc,4);
end
if ~isempty(bref)
    if size(bref,4)>1, bref = mean(bref,4); end
    if any(size(bref)~=nn(1:3))
        warning('size bc and bref mismatch');
        disp(size(bc)); disp(size(bref));
    end
end


%% calculate mask
do_mask = false;
if (size(mask,1)==nn(1)) &&  (size(mask,2)==nn(2)) && (size(mask,3)==nn(3))
    do_mask = true;
end
if length(mask)==1
    if mask(1,1,1)
        do_mask = true;
        if dim==3  % 3D
            mask = calc_mask(sqrt(mean(mean(bc.*conj(bc),4),5)),1,[],[],[],verb>0);
        else       % 2D
            mask = false(nn(1:3));
            for l3=1:mtx(3)
                mask(:,:,l3) = calc_mask(sqrt(mean(mean(bc(:,:,l3,:,:).* ...
                    conj(bc(:,:,l3,:,:)),4),5)),1,[],[],[],verb>0);
            end
        end
    end
end


%% down-sampling
if ~isempty(mtx)
    if verb>0
        fprintf('Resampling from (%d,%d,%d,1,%d) to (%d,%d,%d,1,%d)\n',...
            nn(1:3),nc,mtx,nc); 
    end
    bx = complex(zeros(mtx(1),mtx(2),mtx(3),1,nc,class(bc)));
    for lc=1:nc
        bx(:,:,:,1,lc) = truma(bc(:,:,:,1,lc),true,mtx);
    end
    bc = bx;
    clear bx;
    if ~isempty(bref)
        bref = truma(bref,true,mtx);
    end
    if do_mask
        mask_lowres = truma(single(mask),true,mtx)>0.2;
    end
end


%% calculate sensitivity maps
if verb>0, fprintf('Calculating sensitivity maps\n'); end
if isempty(bref)
    bref = sqrt(sum(bc.*conj(bc),5));
else
    bref = abs(bref);
end
smap = bsxfun(@rdivide,conj(bc),bref);
if do_mask
    smap = bsxfun(@times,smap,mask_lowres);
end


%% checks
if any(isnan(smap(:)))
    warning('isnan(smap) -> 0');
    smap(isnan(smap)) = 0;
end
if verb>1, sub_plotting(smap,dim,111); end


%% polynomial smoothing
if do_fiex>0
    if verb>0, fprintf('%dD Polynomial smoothing (fiex3d)\n',dim); end
    if do_mask && (do_fiex>1)
        ex_rounds = 1;
    else
        ex_rounds = 0;
    end
    if dim==3
        for lc=1:nc
            smap(:,:,:,1,lc) = fiex3d(smap(:,:,:,1,lc),1,4,'','','','',...
                ex_rounds,(verb>1));
        end
    else
        for l3=1:mtx(3)
            for lc=1:nc
                smap(:,:,l3,1,lc) = fiex3d(smap(:,:,l3,1,lc),1,4,'','','','',...
                    ex_rounds,(verb>1));
            end
        end
    end
    if verb>1, sub_plotting(smap,dim,113); end
end
smap(isnan(smap)) = 0;


%% up-sampling
if ~isempty(zf)
    sx = complex(zeros(zf(1),zf(2),zf(3),1,nc,class(bc)));
    for lc=1:nc
        sx(:,:,:,1,lc) = truma(smap(:,:,:,1,lc),true,zf);
    end
    smap = sx;
    clear sx;
    if do_mask
        if zf~=nn(1:3)
           mask = truma(single(mask),true,zf)>0.3;
        end
        smap = bsxfun(@times,smap,mask);
    end
end


%% normalise: sum-over-coils=1
smap = bsxfun(@rdivide,smap,sqrt(mean(smap.*conj(smap),5)));
smap(isnan(smap)) = 0;
if verb>1, sub_plotting(smap,dim,115); end


%% final checks
if any(isnan(smap(:)))
    warning('isnan(smap) -> 0');
    smap(isnan(smap)) = 0;
end


end      % mri_sens_map_fiex.m


%% sub-functions
% plotting
function sub_plotting(smap,dim,figid)
figure(figid); clf; set(gcf,'name','abs(smap)');
if dim>2
    imagesc_ind3d(squeeze(abs(smap)));
else
    imagesc_row(squeeze(abs(smap)),'','','',true);
end
% title('sensitivity map (magnitude)')

figure(figid+1); clf; set(gcf,'name','angle(smap)');
if dim>2
    imagesc_ind3d(squeeze(angle(smap)));
else
    imagesc_row(squeeze(angle(smap)),'','','',true);
end
% title('sensitivity map (phase)')
drawnow
end      % sub_plotting
