function [bb,smap] = mri_coil_combine(b,noise_thresh,time_combine,...
    do_fiex,mtx,verb)
%MRI_COIL_COMBINE  SVDS Coil combination for MRI
%   Determines sensitivity maps of individual coil elements, and optimally
%   combines individual Rx channels
%
%[bb,smap] = mri_coil_combine(b,noise_thresh,time_combine,do_fiex,mtx,verb)
%           b   Spatially reconstructed data 
%               2D=(nx,ny,nt,nm,ns,nc)   recon_spiral
%               3D=(nx,ny,nz,nt,nc)      recon_grid3d
%               spatial=nx,ny,nz/ns; #time steps=nt, #metabolite=nm
%noise_thresh   Percent of data classified as noise             (20)
%               if 0 -> skip
%time_combine   Average time steps:                             ('all')
%               'all'=mean over all time steps
%               'max'=weight and sum according to spatial maximum
%               [time steps;intensities]
%     do_fiex   Polynomial smoothing of sensitivity maps        (true)
%         mtx   Resample to matrix size for calculations        ([])
%        verb   Verbose (>0) + plotting (>1 &>2)                (1)
%
%          bb   Coil-combined images (same size as b w/ nc=1)
%         smap  Sensitivity map      (same size as b w/ nt=nm=nc=1)
%
% 10/2019 Rolf Schulte
if nargin<1, help(mfilename); return; end

%% input + default parameters
if ~exist('noise_thresh','var'); noise_thresh = []; end
if isempty(noise_thresh),        noise_thresh = 20; end
if noise_thresh<0,   warning('noise_thresh(=%g)<0%%',noise_thresh); end
if noise_thresh>100, warning('noise_thresh(=%g)>100%%',noise_thresh); end
if ~exist('time_combine','var'); time_combine = ''; end
if isempty(time_combine),        time_combine = 'all'; end
if ~exist('do_fiex','var');      do_fiex = []; end
if isempty(do_fiex),             do_fiex = true; end
if ~exist('mtx','var');          mtx = []; end
if ~exist('verb','var');         verb = []; end
if isempty(verb),                verb = 1; end
bb = b;       % copy original data


%% misc pars + reshaping
switch length(size(b))
    case 5     % 3D data reconstructed with recon_grid3d
        dim = 3;
        [nx,ny,nz,nt,nc] = size(b);
        nm = 1;
        dim_coil = 5;
    case 6     % 2D data reconstructed with recon_spiral
        dim = 2;
        [nx,ny,nt,nm,nz,nc] = size(b);
        dim_coil = 6;
    otherwise, error('length(size(b))(=%d)~=5 or 6',length(size(b)));
end


%% time averaging
if time_combine(1)>0
    if verb>0, fprintf('Averaging along time dimension\n'); end
    if ischar(time_combine)
        switch time_combine
            case 'all'
                time_combine = [(1:nt);ones(1,nt)/nt];
            case 'max'
                tmp = sqrt(mean(b.*conj(b),dim_coil));
                if dim==2
                    tmp = max(max(max(max(tmp,[],1),[],2),[],4),[],5);
                else
                    tmp = max(max(max(tmp,[],1),[],2),[],3);
                end
                time_combine = [(1:nt);squeeze(tmp).'];
            otherwise
                error('time_combine(=%s): unkown option',time_combine);
        end
    end
    if size(time_combine,1)==1
        tmp = size(time_combine,2);
        time_combine(2,:) = ones(1,tmp)/tmp;
    end
    if verb>0
        tmp = num2str(time_combine,'%.3g  ');
        fprintf('time_combine=\n\t%s\n\t%s\n',tmp(1,:),tmp(2,:)); 
    end
    if dim==2
        ww = permute(time_combine(2,:),[1 3 2]); 
        b = sum(bsxfun(@times,b(:,:,time_combine(1,:),:,:,:),ww),3);
    else
        ww = permute(time_combine(2,:),[1 3 4 2]); 
        b = sum(bsxfun(@times,b(:,:,:,time_combine(1,:),:),ww),4);        
    end
    nt = 1;
    
    if verb>2
        figure
        bbabs=sqrt(mean(b.*conj(b),dim_coil));
        imagesc_row(squeeze(bbabs),'','ind','',true);
    end
end


%% reshaping data
if dim==2
    if verb>0, fprintf('2D input data: reshaping\n'); end
    br = complex(zeros(nx,ny,nz,nt*nm,nc,class(b)));
    for lz=1:nz
        for lt=1:nt
            for lm=1:nm
                for lc=1:nc
                    ii = lt+(lm-1)*nt;
                    br(:,:,lz,ii,lc) = b(:,:,lt,lm,lz,lc);
                end
            end
        end
    end
    b = br;
    clear br;
else
    if verb>0, fprintf('3D input data\n'); end
end


%% noise decorrelation/pre-whitening
if noise_thresh>0
    % determine noise covariance psi
    if verb>0, fprintf('Calculating noise covariance matrix psi\n'); end
    psi = sub_noise_covariance(bb,noise_thresh,dim_coil,90*(verb>2));
    
    
    % pre-whitening
    % literature: ISMRM2014, sunrise, Michael S Hansen; Nuts & Bolts...
    % http://hansenms.github.io/sunrise/sunrise2014/
    if verb>0, fprintf('Noise pre-whitening/decorrelation\n'); end
    L = chol(psi,'lower');
    L_inv = inv(L);
    nn = size(b);
    b = L_inv*reshape(permute(b,[5 1 2 3 4]),[nc prod(nn(1:4))]);
    b = permute(reshape(b,[nc,nn(1:4)]),[2 3 4 5 1]);
    
    if dim==2
        nn = size(bb);
        bb = L_inv*reshape(permute(bb,[6 1 2 3 4 5]),[nc prod(nn(1:5))]);
        bb = permute(reshape(bb,[nc,nn(1:5)]),[2 3 4 5 6 1]);
    else
        nn = size(bb);
        bb = L_inv*reshape(permute(bb,[5 1 2 3 4]),[nc prod(nn(1:4))]);
        bb = permute(reshape(bb,[nc,nn(1:4)]),[2 3 4 5 1]);
    end
    
    % determine noise covariance after decorrelation
    if verb>0
        fprintf('Calculating noise covariance matrix after decorrelation\n'); 
    end
    sub_noise_covariance(bb,noise_thresh,dim_coil,91*(verb>2));
end


%% down-sampling
if ~isempty(mtx)
    if verb>0, fprintf('Resampling from (%d,%d,%d) to ',nx,ny,nz); end
    nx1 = nx; ny1 = ny; nz1 = nz;
    mtx = mtx(:).';
    nx = mtx(1,1);
    if size(mtx,2)>1, ny = mtx(1,2); end
    if size(mtx,2)>2, nz = mtx(1,3); end
    if verb>0, fprintf('(%d,%d,%d)\n',nx,ny,nz); end
    bx = complex(zeros(nx,ny,nz,nt*nm,nc));
    for lt=1:(nt*nm)
        for lc=1:nc
            bx(:,:,:,lt,lc) = truma(b(:,:,:,lt,lc),true,[nx ny nz]);
        end
    end
    b = bx;
    clear bx;
end


%% determine sensitivity map (combination weights)
if verb>0, fprintf('Calculating sensitivity maps via svds\n'); end
smap = complex(zeros(nx,ny,nz,1,nc));       % init matrix
ss = zeros(nx,ny,nz);
for lx=1:nx
    for ly=1:ny
        for lz=1:nz
            A = double(squeeze(b(lx,ly,lz,:,:)));
            [~,S,V,flag] = svds(A,1);
            smap(lx,ly,lz,1,:) = permute(V,[2 3 4 5 1]);
            ss(lx,ly,lz) = S;
            
            if flag~=0
                fprintf('lx=%d; ly=%d; lz=%d; flag=%d\n',lx,ly,lz,flag);
            end
        end
    end
end


%% remove average phase
phafu = sum(smap,5); phafu = phafu./abs(phafu);
smap = smap.*conj(phafu);


%% masking data
if false
    mask_thresh = 1.3
    mask = ss>mask_thresh;
    smap = bsxfun(@times,smap,mask);
    ex_rounds = 2;
else
    ex_rounds = 0;
end
if verb>1
    figure(100); clf
    imagesc_row(squeeze(real(smap)),'','','',true);
    drawnow
end


%% polynomial smoothing
if do_fiex
    for lc=1:nc
        smap(:,:,:,1,lc) = fiex3d(smap(:,:,:,1,lc),1,4,'','','','',ex_rounds);
    end
    if verb>1
        figure(101); clf
        imagesc_row(squeeze(real(smap)),'','','',true);
        drawnow
    end
end


%% normalise: sum-over-coils=1
% smap = smap./repmat(sum(abs(smap),5),[1 1 1 1 nc]); % normalise
smap = bsxfun(@rdivide,smap,sqrt(mean(smap.*conj(smap),5)));
if verb>1
    figure(102); clf
    imagesc_row(squeeze(real(smap)),'','','',true);
    drawnow
end


%% up-sampling
if ~isempty(mtx)
    nx = nx1; ny = ny1; nz = nz1;
    sx = complex(zeros(nx,ny,nz,1,nc));
    for lc=1:nc
        sx(:,:,:,1,lc) = truma(smap(:,:,:,1,lc),true,[nx ny nz]);
    end
    smap = sx;
    clear sx;
end


%% actual phasing and weighting plus coil combination
if verb>0, fprintf('Actual phased and weighted coil combination\n'); end
if dim==2
    smap = reshape(smap,[nx,ny,1,1,nz,nc]);
end
bb = sum(bsxfun(@times,bb,smap),dim_coil);


end      % main function mri_coil_combine.m


%% sub-functions
function psi = sub_noise_covariance(br,noise_threshold,dim_coil,figid)

%% reshape br to [nc,nrest]
nn = size(br);
ii = (1:length(nn));
ii = ii(ii~=dim_coil);

br = permute(br,[dim_coil ii]);
nc = size(br,1);
br = reshape(br,[nc prod(nn(ii))]);


%% calculate standard-deviation and threshold data to extract noise
% stdb = std(br(:));
% ii = sqrt(mean(br.*conj(br),1))<noise_threshold*stdb;
% fprintf('#noise samples: %.3g%%\n',sum(ii)/size(ii,2)*100);

iend = ceil(noise_threshold/100*size(br,2));
if iend<1, error('iend(=%g)<1',iend); end
if iend>size(br,2), error('iend(=%g)>size(br,2)(=%g)',iend,size(br,2)); end

[~,ii] = sort(sqrt(mean(br.*conj(br),1)));
ii = ii(1,1:iend);


noise = br(:,ii);
fprintf('#noise samples: %.3g%%\n',size(noise,2)/size(br,2)*100);

%% calculate noise covariance and scale to mean of diagonal element=1
psi = abs(noise*noise');
psi = psi/mean(diag(psi));


%% plotting
if figid>0
    figure(figid); clf
    imagesc(psi,[0 2]); 
    colorbar;
    title('Noise covariance');
    drawnow
end

end      % sub-function sub_coise_covariance
