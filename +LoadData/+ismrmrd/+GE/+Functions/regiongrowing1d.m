function mask = regiongrowing1d(d,dim,thresh,seed,do_normalise)
%  mask = regiongrowing1d(d,dim,thresh,seed,do_normalise)
%                                                               (default)
%           d  Input data (real-valued; arbitrary dimensions)
%         dim  Dimension along which to perform region growing  (1)
%      thresh  Threshold for stopping                           (0.1)
%        seed  Seed point (0=use maximum)                       (0)
%do_normalise  0=off; 1=each row individually; 2=global         (1)
%        mask  Output mask (logical)
%
% 1/2019  Rolf Schulte
if (nargin<1), help(mfilename); return; end


%% default input parameters + checks
if ~exist('dim','var'),    dim = []; end
if isempty(dim),           dim = 1; end
if length(dim)~=1,         error('length(dim)~=1'); end
if ~exist('thresh','var'), thresh = []; end
if isempty(thresh),        thresh = 0.1; end
if length(thresh)~=1,      error('length(thresh)~=1'); end
if ~exist('seed','var'),   seed = []; end
if isempty(seed),          seed = 0; end
if length(seed)~=1,        error('length(seed)~=1'); end
do_max = false;
if seed==0, do_max = true; end
if ~exist('do_normalise','var'),  do_normalise = []; end
if isempty(do_normalise),         do_normalise = 1; end
if length(do_normalise)~=1,       error('length(do_normalise)~=1'); end

if ~isreal(d)
    warning('input data d complex; taking abs(d)');
    d = abs(d);
end
n = size(d);
if dim>length(n), error('dim(=%g)>length(size(d))(=%g)',dim,length(n)); end
if dim<1, error('dim(=%g)<1',dim); end


%% permute data -> always perform search along dim1
if dim>1
    prmt_vec = 1:length(n);
    prmt_vec(1) = dim;
    prmt_vec(dim) = 1;
    d = permute(d,prmt_vec);
end


%% normalisation
switch do_normalise
    case 0    % skip normalisation
    case 1    % normalise each row individually
        dmax = max(d,[],1);
        d = bsxfun(@rdivide,d,dmax);
    case 2    % normalise data globally
        d = d/max(d(:));
end

%%
mask = false(size(d));
[n1,n2] = size(d);
for l2=1:n2
    if do_max, [~,seed] = max(d(:,l2),[],1); end
    for l1=seed:n1
        if d(l1,l2)<thresh, break; end
        mask(l1,l2) = true;
    end
    for l1=(seed-1):-1:1
        if d(l1,l2)<thresh, break; end
        mask(l1,l2) = true;
    end
end


%% permute mask back 
if dim>1
    mask = permute(mask,prmt_vec);
end

%% final size check
for l=1:length(n)
    if (n(l)~=size(mask,l)) && (l~=dim)
        warning('n(%g)(=%g)~=size(mask,%g)(=%g)',n(l),size(mask,l));
    end
end


end      % main function regiongrowing1d.m
