function smap = mri_sens_map_espirit(dd,zf,dim,ksize,ncal,do_vcc,rm_ph0,verb)
%MRI_SENS_MAP_ESPIRIT  Extract 2D or 3D receive coil sensitivity maps 
%    smap = mri_sens_map_espirit(bc,zf,dim,ksize,ncal,do_vcc,rm_ph0,verb)
%      bc   Surface coil image    [mtx(1),mtx(2),mtx(3),nt,nc]
%      zf   Desired final smap matrix size                            (mtx)
%     dim   Data dimension: 2 or 3
%   ksize   k-space kernel size                                         (6)
%    ncal   Size of calibration area                                   ([])
%  do_vcc   Virtual conjugate coils: remove object phase             (true)
%  rm_ph0   Remove average coil phase to avoid phase wraps           (true)
%    verb   Verbose (>0) + plotting (>1 &>2)                            (1)
%    smap   sensitivity maps      [zf(1),zf(2),zf(3),1,nc]
%
% Literature: mrm71_990_2014_uecker.pdf, mrm77_1201_2017_uecker.pdf
% 12/2022 Rolf Schulte
if (nargin<1), help(mfilename); return; end

nsvd = 1;          % #sens maps


%% input parameters
if ~exist('zf','var'),       zf = []; end
if ~exist('dim','var'),      dim = []; end
if ~exist('ksize','var'),    ksize = []; end
if isempty(ksize),           ksize = 6; end
if ~exist('ncal','var'),     ncal = []; end
if length(ncal)==1,          ncal = repmat(ncal,[1 3]); end
if ~exist('do_vcc','var'),   do_vcc = []; end
if isempty(do_vcc),          do_vcc = true; end
if ~exist('rm_ph0','var'),   rm_ph0 = []; end
if isempty(rm_ph0),          rm_ph0 = true; end
if ~exist('verb','var'),     verb = []; end
if isempty(verb),            verb = true; end


%% checks
if verb>0
    fprintf('Calculating Rx coil sensitivities: eSPIRIT\n');
end
nn = size(dd);                         % image size
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
if ~isempty(zf)
    if length(zf)==dim
        warning('length(zf)(=%g)==dim(=%g)',length(zf),dim);
    end
else
    zf = nn(1:dim);
end
if dim==2
    if isempty(ncal), ncal = [32 32 nn(3)]; end
    ncal(3) = nn(3);
    zf(3) = nn(3);
else
    if isempty(ncal), ncal = [16 16 16]; end
end
if isreal(dd), warning('bc is real'); end
if nt>1
    fprintf('nt(=%d) > 1\n',nt);
    if dim==3
        warning('nt>1 && 3D; averaging over time');
        dd = mean(dd,4);
    end
end
for l=1:3
    if ncal(l)>nn(l)
        warning('ncal(%d)(=%d)>nn(%d)(=%d): capping',...
        l,ncal(l),l,nn(l));
        ncal(l) = nn(l);
    end
end


%% remove average phase from individual coils to avoid phase wraps
if rm_ph0
    ph0 = zeros(1,1,1,nt,nc);
    bbabs = sqrt(mean(mean(dd.*conj(dd),5),4));
    mask = bbabs>0.7*mean(bbabs(:));
    for lc=1:nc
        for lt=1:nt
            ph = angle(dd(:,:,:,lt,lc));
            ph = ph(mask);
            ph0(1,1,1,lt,lc) = (mean(ph(ph>0)) - mean(ph(ph<0)))/2;
        end
        ph0 = mean(ph0,4);
        if verb>1
            fprintf('mean(phase) coil%d = %g [rad]\n',lc,ph0(1,1,1,1,lc));
        end
    end
    clear ph bbabs mask
    dd = bsxfun(@times,dd,exp(-1i*ph0));
end


%% virtual conjugate coil: half the phase
if do_vcc
    phafu = dd./abs(dd);
    dd = dd.*conj(sqrt(phafu));
end


%% fft image to k-space
for ld=1:dim
    dd = ifftshift(fft(fftshift(dd,ld),[],ld),ld);
end
dd = truma(dd,false,[ncal nt nc],[],true);   % truncate to central k-space


%% virtual conjugate coil: actual doubling+conjugation of coils
if do_vcc
    if dim==2
        dd2 = conj(dd([1,end:-1:2],[1,end:-1:2],:,:,:));
    else
        dd2 = conj(dd([1,end:-1:2],[1,end:-1:2],[1,end:-1:2],:,:));
    end
end


n1 = (ncal(1)-ksize);
n2 = (ncal(2)-ksize);
if dim==2, n3 = nn(3);
else,      n3 = (ncal(3)-ksize); end
mtx = ncal;

if do_vcc, nc2 = 2*nc;
else,      nc2 = nc;
end
smap = complex(zeros(mtx(1),mtx(2),mtx(3),1,nc2,nsvd,class(dd)));


%% actual kernel creation, fft to spatial domain, svd fit
if dim==2     % 2D
    for l3=1:n3
        gq = complex(zeros(n1*n2*nt,nc2,mtx(1)*mtx(2),class(dd)));

        for l2=1:n2
            for l1=1:n1
                ll = l1 + (l2-1)*n1;
                % fprintf('l1=%d l2=%d ll=%d\n',l1,l2,ll);
                ii1 = l1-1+(1:ksize);
                ii2 = l2-1+(1:ksize);
                gkern = reshape(dd(ii1,ii2,l3,:,:),[ksize ksize nt nc]);
                if do_vcc
                    gkvcc = reshape(dd2(ii1,ii2,l3,:,:),[ksize ksize nt nc]);
                    gkern = cat(4,gkern,gkvcc);
                end
                gkern = truma(gkern,false,[mtx(1) mtx(2) nt nc2]);
                gki = fftshift(fftshift(ifft(ifft(ifftshift(ifftshift(gkern,...
                    2),1),[],2),[],1),2),1);
                ii = ll+(0:(nt-1))*n1*n2;
                gq(ii,:,:) = reshape(permute(gki,[3 4 1 2]),[nt nc2 mtx(1)*mtx(2)]);

            end
        end
        for l2=1:mtx(2)
            for l1=1:mtx(1)
                ll = l1 + (l2-1)*mtx(1);
                [~,~,V] = svd(gq(:,:,ll),0);
                smap(l1,l2,l3,1,:,:) = reshape(V(:,1:nsvd),[1 1 1 1 nc2 nsvd]);
            end
        end
    end
else          % 3D
    gq = complex(zeros(n1*n2*n3,nc2,mtx(1)*mtx(2)*mtx(3),class(dd)));
    for l3=1:n3
        for l2=1:n2
            for l1=1:n1
                ll = l1 + (l2-1)*n1 + (l3-1)*n1*n2;
                % fprintf('l1=%d l2=%d l3=%d ll=%d\n',l1,l2,l3,ll);
                ii1 = l1-1+(1:ksize);
                ii2 = l2-1+(1:ksize);
                ii3 = l3-1+(1:ksize);
                gkern = reshape(dd(ii1,ii2,ii3,1,:),[ksize ksize ksize nc]);
                if do_vcc
                    gkvcc = reshape(dd2(ii1,ii2,ii3,1,:),[ksize ksize ksize nc]);
                    gkern = cat(4,gkern,gkvcc);
                end
                gkern = truma(gkern,false,[mtx(1) mtx(2) mtx(3) nc2]);
                gki = fftshift(fftshift(fftshift(ifft(ifft(ifft(...
                    ifftshift(ifftshift(ifftshift(gkern,...
                    3),2),1),[],3),[],2),[],1),3),2),1);
                gq(ll,:,:) = reshape(permute(gki,[4 1 2 3]),[1 nc2 mtx(1)*mtx(2)*mtx(3)]);
            end
        end
    end
    for l3=1:mtx(3)
        for l2=1:mtx(2)
            for l1=1:mtx(1)
                ll = l1 + (l2-1)*mtx(1) + (l3-1)*mtx(1)*mtx(2);
                [~,~,V] = svd(gq(:,:,ll),0);
                smap(l1,l2,l3,1,:,:) = reshape(V(:,1:nsvd),[1 1 1 1 nc2 nsvd]);
                % [~,~,V] = svds(gq(:,:,ll),nsvd);
                % smap(l1,l2,:,:) = reshape(V,[1 1 nc nsvd]);
            end
        end
    end
end      % if dim==2


%% remove phase of coil #1 (=0 or pi) to remove phase ambiguity
smap = bsxfun(@times,smap,conj(smap(:,:,:,:,1,:)./abs(smap(:,:,:,:,1,:))));


%% adjust phase to yield real-only image
if do_vcc
    % phase difference physical-virtual coils: 
    % same as calculating and subtracting phi in Eq. 5 in mrm77_1201_2017_uecker.pdf
    smap = smap(:,:,:,:,1:nc,:).*conj(smap(:,:,:,:,nc+(1:nc),:))./...
        abs(smap(:,:,:,:,1:nc,:));
end


%% checks
sum_isnan = sum(isnan(smap(:)));
if sum_isnan>0
    warning('isnan(smap)(=%d) -> 0',sum_isnan);
    smap(isnan(smap)) = 0;
end
if verb>1, sub_plotting(smap,dim,111); end


%% up-sampling
if ~isempty(zf)
    smap = truma(smap,true,[zf 1 nc nsvd]);
end

%% add initial average coil phase again
if rm_ph0
    smap = bsxfun(@times,smap,exp(-1i*(ph0-ph0(1,1,1,1,1))));
    % smap = bsxfun(@times,smap,exp(-1i*(ph0)));
end


%% normalise: sum-over-coils=1
smap = bsxfun(@rdivide,smap,sqrt(mean(smap.*conj(smap),5)));
smap(isnan(smap)) = 0;


%% final checks
if any(isnan(smap(:)))
    warning('isnan(smap) -> 0');
    smap(isnan(smap)) = 0;
end
if verb>1, sub_plotting(smap,dim,113); end


end      % mri_sens_map_espirit.m


%% sub-functions
% plotting
function sub_plotting(smap,dim,figid)
rshp = false;
if size(smap,6)>1, smap = smap(:,:,:,:,:); end
if (dim==2) && length(size(squeeze(smap)))<4, rshp = true; end
figure(figid); clf
if dim>2
    imagesc_ind3d(squeeze(abs(smap)));
else
    imagesc_row(squeeze(abs(smap)),'','',rshp,true,[],300);
end
set(figid,'name','abs(smap)')

figure(figid+1); clf
if dim>2
    imagesc_ind3d(squeeze(angle(smap)));
else
    imagesc_row(squeeze(angle(smap)),'','',rshp,true,[],300);
end
set(figid+1,'name','angle(smap)')
drawnow

end      % sub_plotting


