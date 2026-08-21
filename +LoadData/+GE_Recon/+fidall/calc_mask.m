function mask = calc_mask(bb,method,par,seed,denoise,verb,area)
%CALC_MASK  Calculate mask based on dynamic thresholding
%   mask = calc_mask(bb,method,par,seed,denoise,verb,area)
%     bb   Images
% method   1=thresholding; 2=regiongrowing; 3=activecontour             (1)
%    par   Thresholding parameters
%          method==1: scale=mean of par(1) to 1 largest data
%                     threshold = par(2)*scale                  ([0.9 0.2])
%          method==2: rgdiff=par(1)*(mean+median)/2                   (0.3)
%          method==3: [#iter,SmoothFactor,ContractionBias] ([300 0.3 -0.1])
%   seed   regiongrowing: Seed point or #seed points                    (2)
%          activecontour: #lines along x-y-z,low-high      ([0 0;0 0;0 10])
%denoise   Shrink + dilate to remove small patches                   (true)
%   verb   Verbose: 0=off, 1=verbose, 2=verbose+plotting                (1)
%   area   Range of mask area [%%]: (min max)                          ([])
%
%  9/2024 Rolf Schulte
if nargin<1, help(mfilename); return; end


%% default parameters
if ~exist('method','var'), method = []; end
if isempty(method),        method = 1; end
if ~exist('par','var'),    par = []; end
if ~exist('seed','var'),   seed = []; end
if ~exist('denoise','var'), denoise = []; end
if isempty(denoise),       denoise = true; end
if ~exist('verb','var'),   verb = []; end
if isempty(verb),          verb = 1; end
if length(verb)<2,         verb = [verb 1]; end
if ~exist('area','var'),   area = []; end


%% misc parameters + checks
nn = size(bb);
dim = length(nn);
if dim>3, warning('dim(=%d)>3',dim); end
if dim<2, warning('dim(=%d)<2',dim); end
if ~isreal(bb)
    warning('bb complex; taking magnitude');
    bb = abs(bb);
end
    

%% aply threshold
if ~isempty(par)
    if par(1)<1d-10, method = 0; end
end
switch method
    case 0
        fprintf('Setting mask=true(nn)\n');
        mask = true(nn);
        return;
    case 1
        fprintf('Calculating mask via threshold\n');
        if isempty(par), par = [0.9 0.2]; end
        % calculate threshold
        bs = sort(bb(:));
        ii = round(prod(nn)*par(1)):prod(nn);
        scale = mean(bs(ii));
        if isnan(scale)
            warning('isnan(scale): setting to max');
            scale = max(bb(:));
        end
        thresh = par(2)*scale;
        % actual thresholding
        % mask = bb>thresh;
        mask = imfill(bb>thresh,'holes');
    case 2
        fprintf('Calculating mask via regiongrowing\n');
        if isempty(seed), seed = 2; end
        if isscalar(seed)
            nseed = seed;
            seed = ones(1,dim);
            if nseed>1, seed = [seed;nn]; end
            if nseed>2, warning('nseed(=%d)>2; not implemented',nseed); end
        end
        if isempty(par), par = 0.3; end
        rgdiff = par(1)*(mean(bb(:))+median(bb(:)))/2;
        mask = ~regiongrowing(bb,rgdiff,seed(1,1:dim));
        for l1=2:size(seed,1)
            mask = mask & ~regiongrowing(bb,rgdiff,seed(l1,1:dim));
        end
    case 3
        fprintf('Calculating mask via activecontour\n');
        if isempty(par), par = 300; end
        if length(par)<2, par(2) = 0.3; end
        if length(par)<3, par(3) = -0.1; end
        mask_init = false(nn);
        if isempty(seed)
            switch dim
                case 2, seed = [0,0;0,10]; 
                case 3, seed = [0,0;0,0;0,10];
            end
        end
        for ll=1:dim, ii{ll} = [1:seed(ll,1), ((-seed(ll,2)+1):0)+nn(ll)]; end
        mask_init(ii{1},:,:) = true;
        mask_init(:,ii{2},:) = true;
        if dim>2, mask_init(:,:,ii{3}) = true; end
        mask = ~activecontour(bb,mask_init,par(1),'Chan-Vese',...
            'SmoothFactor',par(2),'ContractionBias',par(3));
    otherwise
        warning('method(=%g) unkown',method);
        mask = true(nn);
        return;
end


%% shrink + dilate to remove small patches
if denoise
    mask = imdilate(~mask,strel_sphere(2,dim));
    mask = ~logical(mask);
    mask = imdilate(mask,strel_sphere(2,dim));
end


%% print info
per_mask = sum(mask(:))/prod(nn(1:dim));
if verb(1)>0
    fprintf('mask area = %g [%%]\n',per_mask*100);
end


%% recalculate, if mask area not within given range
if ~isempty(area)
    if length(area)<2, area = [area 100]; end
    if (area(2)<area(1)) || (area(1)<0) || (area(2)>100)
        error('area(=%g %g)',area(1),area(2));
    end
    lp = 3-method;
    if length(par)<3, par(3) = 0; end
    if (per_mask*100>area(2)) && (par(3)~=2)
        if verb
            fprintf('Increasing par(%d) from %g to %g\n',lp,par(lp),par(lp)*1.1);
        end
        par(lp) = par(lp)*1.1;
        par(3) = 1;
        mask = calc_mask(bb,method,par,seed,denoise,verb,area);
    end
    if (per_mask*100<area(1)) && (par(3)~=1)
        if verb
            fprintf('Reducing par(%d) from %g to %g\n',lp,par(lp),par(lp)*0.9);
        end
        par(lp) = par(lp)*0.9;
        par(3) = 2;
        mask = calc_mask(bb,method,par,seed,denoise,verb,area);
    end
end
per_mask = sum(mask(:))/prod(nn(1:dim));
if per_mask>0.9, warning('mask area (=%g) > 90%',per_mask*100); end
if per_mask<0.1, warning('mask area (=%g) < 10%',per_mask*100); end


%% plotting
if verb(1)>1
    figure(20); clf;
    if dim==2
        imagesc(mask);
        axis image
        set(gca,'Position',[0 0 1 1]);
    else
        if verb(2)==1
            pltseed = [];
            if (verb(1)>2) && (method==2)
                pltseed = seed(1,:);
            end
            imagesc_ind3d(mask,pltseed,'','',~isempty(pltseed),true);
        else
            imagesc_row(mask,[],[],true);
        end
    end
    drawnow
end


end      % calc_mask.m
