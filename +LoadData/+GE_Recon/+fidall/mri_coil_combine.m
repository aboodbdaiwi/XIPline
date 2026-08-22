function [bb,smap] = mri_coil_combine(b,nper,tcomb,method,mtx,mask,verb,dim)
%MRI_COIL_COMBINE  SVDS Coil combination for MRI
%   Determines sensitivity maps of individual coil elements, and optimally
%   combines individual Rx channels
%
%[bb,smap] = mri_coil_combine(b,nper,tcomb,method,mtx,mask,verb,dim)
%            b   Spatially reconstructed data 
%                2D=(nx,ny,nt,nm,ns,nc)   recon_spiral
%                2D/3D=(nx,ny,nz,nt,nc)   recon_grid3d/recon_qti
%            spatial=nx,ny,nz/ns; #time steps=nt, #metabolite=nm
%     nper   Noise decorrelation - percent of data class. as noise     (10)
%            if 0 -> skip
%t_combine   Average time steps:                                    ('all')
%            'all'=mean over all time steps
%            'max'=weight and sum according to spatial maximum
%            'first'=take only first time step
%            [time steps;intensities]
%   method   Rx coil sensitivity map methods                            (2)
%            0-2: Lowres + Polynomial smoothing of sensitivity maps:
%                 0=fiex off; 1=inter; 2=inter+extrapolate
%            3=eSPIRIT
%            4=Walsh
%            5=RMS
%            or smap as direct input
%      mtx   Resample to matrix size for calculations                  ([])
%     mask   Calculate over masked area                                ([])
%     verb   Verbose (>0) + plotting (>1 &>2)                           (1)
%      dim   Dimensions: 2 or 3; default from size(b)
%
%       bb   Coil-combined images (same size as b w/ nc=1)
%     smap   Sensitivity map      (same size as b w/ nt=nm=nc=1)
%
%  9/2024 Rolf Schulte
if nargin<1, help(mfilename); return; end


%% input + default parameters
if ~exist('nper','var'),    nper = []; end
if isempty(nper),           nper = 10; end
if nper<0,                  warning('nper(=%g)<0%%',nper);   end
if nper>100,                warning('nper(=%g)>100%%',nper); end
if ~exist('tcomb','var'),   tcomb = ''; end
if isempty(tcomb),          tcomb = 'all'; end
if ~exist('method','var'),  method = []; end
if isempty(method),         method = 2; end
if ~exist('mtx','var'),     mtx = []; end
if ~exist('mask','var'),    mask = []; end
if ~exist('verb','var'),    verb = []; end
if isempty(verb),           verb = 1; end
if ~exist('dim','var'),     dim = []; end
bb = b;       % copy original data
if verb>0, timerVal = tic; end         % record runtime


%% misc pars + reshaping
if ~isempty(mtx)
    if (mtx(1)>1) && (mtx(1)<8)
        warning('mtx(1)(=%d)<8',mtx(1));
    end
    if (mtx(1)<2), mtx = []; end
end
switch length(size(b))
    case 5     % 3D data recon_grid3d; or 2D/3D recon_qti
        [nx,ny,nz,nt,nc] = size(b);
        nm = 1;
        dim_coil = 5;
        if isempty(dim)
            if nz<10
                dim = 2;
            else
                dim = 3;
            end
        end
    case 6     % 2D data reconstructed with recon_spiral
        [nx,ny,nt,nm,nz,nc] = size(b);
        dim_coil = 6;
        if isempty(dim), dim = 2; end
    otherwise, error('length(size(b))(=%d)~=5 or 6',length(size(b)));
end


%% check is smap provided as input
if ((size(method,1)>1) && (size(method,2)>1))
    if verb>0
        fprintf('size(method)(=%s)>1 1; assuming smap input\n',...
            num2str(size(method)));
    end
    smap = method;
    method = -99;
    tcomb(1) = 0;
    switch dim_coil
        case 5
            dimint = [1 2 3];
            dimchk = 5;
            dimone = 4;
        case 6
            dimint = [1 2];
            dimchk = [5 6];
            dimone = [3 4];
    end
    for ll=dimone
        if size(smap,ll)~=1
            error('size(smap,%d)(=%d)~=1',ll,size(smap,ll));
        end
    end
    for ll=dimchk
        if size(smap,ll)~=size(b,ll)
            error('size(smap,%d)(=%d)~=size(b,%d)(=%d)',ll,size(smap,ll),...
                ll,size(b,ll));
        end
    end
    nb = size(b);
    nintp = size(smap);
    if any(nb(dimint)~=nintp(dimint))
        nintp(dimint) = nb(dimint);
        if verb>0, fprintf('interpolating smap to %s\n',num2str(nintp)); end
        smap = truma(smap,true,nintp);
    end
end


%% RMS coil combination
if method==5
    if verb>0, fprintf('RMS coil combination\n'); end
    bb = sqrt(mean(bb.*conj(bb),dim_coil));
    smap = [];
    return
end


%% time averaging
if tcomb(1)>0
    if verb>0, fprintf('Averaging along time dimension\n'); end
    if ischar(tcomb)
        switch tcomb(1:3)
            case 'all'
                tcomb = [(1:nt);ones(1,nt)/nt];
            case 'max'
                tmp = sqrt(mean(b.*conj(b),dim_coil));
                if dim_coil==6
                    tmp = max(max(max(max(tmp,[],1),[],2),[],4),[],5);
                else
                    tmp = max(max(max(tmp,[],1),[],2),[],3);
                end
                tcomb = [(1:nt);squeeze(tmp).'];
            case 'fir'
                tcomb = [1;1];
            otherwise
                error('tcomb(=%s): unkown option',tcomb);
        end
    end
    if size(tcomb,1)==1
        tmp = size(tcomb,2);
        tcomb(2,:) = ones(1,tmp)/tmp;
    end
    if verb>0
        tmp = num2str(tcomb,'%.3g  ');
        fprintf('tcomb = [%s ; %s] [steps;intensities]\n',tmp(1,:),tmp(2,:)); 
    end
    if dim_coil==6
        ww = permute(tcomb(2,:),[1 3 2]); 
        b = sum(bsxfun(@times,b(:,:,tcomb(1,:),:,:,:),ww),3);
    else
        ww = permute(tcomb(2,:),[1 3 4 2]); 
        b = sum(bsxfun(@times,b(:,:,:,tcomb(1,:),:),ww),4);        
    end
    nt = 1;
    
    if verb>2
        figure(98); clf; set(gcf,'name','bb - rms coil combined');
        bbabs=sqrt(mean(b.*conj(b),dim_coil));
        if dim>2
            imagesc_ind3d(squeeze(bbabs));
        else
            imagesc_row(squeeze(bbabs),'','ind','',true);
        end
    end
end


%% reshaping data
if dim_coil==6
    if verb>0, fprintf('2D input data: reshaping\n'); end
    br = complex(zeros(nx,ny,nz,nt*nm,nc,class(bb)));
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
    if verb>0, fprintf('%dD input data; not reshaping\n',dim); end
end


%% noise decorrelation
if nper>0
    if verb, fprintf('Noise decorrelation\n'); end
    b = mri_noise_decorrelation(b,nper,5,0,verb);
end


%% calculate sensitivity maps
switch method
    case -99       % smap provided as input
    case {0,1,2}   % low-res plus polynomial smoothing
        smap = LoadData.GE_Recon.fidall.mri_sens_map_fiex(b,[],[],dim,method,mtx,mask,verb);
    case 3         % eSPIRIT
        smap = LoadData.GE_Recon.fidall.mri_sens_map_espirit(b,[],dim,[],mtx,true,true,verb);
    case 4         % Walsh
        smap =LoadData.GE_Recon.fidall. mri_sens_map_walsh(b);
    otherwise
        error('method(=%g) unknown',method);
end


%% Applying mask
if method>2
    if isempty(mask), mask = false; end
    if (size(mask,1)==size(smap,1)) && ...
            (size(mask,2)==size(smap,2)) && ...
            (size(mask,3)==size(smap,3))
        smap = bsxfun(@times,smap,mask);
    end
end
if dim_coil==6
    smap = reshape(smap,[nx,ny,1,1,nz,nc]);
end


%% actual phasing and weighting plus coil combination
if verb>0, fprintf('Actual phased and weighted coil combination\n'); end
bb = sum(bsxfun(@times,bb,smap),dim_coil);


%% plotting
if verb>2
    figure(99); clf; set(gcf,'name','bb - smap coil combined');
    if dim>2
        imagesc_ind3d(squeeze(abs(bb)));
    else
        imagesc_row(squeeze(abs(bb)),'','inm','',true);
    end
end


%% total execution time
if verb>0
    fprintf('mri_coil_combine: runtime = %.2f [s]\n',toc(timerVal));
end

end      % mri_coil_combine.m
