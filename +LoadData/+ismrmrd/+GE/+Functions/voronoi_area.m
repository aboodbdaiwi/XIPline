function dcf_out = voronoi_area(k,verb,out_corr,dcf_thresh)
%VORONOI_AREA  Calculate 2D/3D k-space density via Voronoi triangulation
%       dcf = voronoi_area(k,verb,out_corr,dcf_thresh)
%
%         k   k-space trajectory  ([-0.5 0.5]*mtx)
%             size(k)=(1,#kpts); kx=real(k);ky=imag(k);  OR
%             size(k)=(1,#kpts,dim)
%      verb   Verbose mode (0=off;1=fprintf,2=fprintf+plotting)      (1)
%  out_corr   Correction method for points at edge                   (2)
%             0=set to 0; 1=use radius; 2=use last point
%             3=nearest neigbour; 4=replace outlier with k
%dcf_thresh   Threshold: set dcf>dcf_thresh to dcf_thresh            ([])
%       dcf   area assosiated with each k-space point
%
% 7/2019 Rolf Schulte
% See also VORONOIN, CONVHULL, UNIQUE.
if (nargin<1), help(mfilename); return; end

fufa = 0.6;        % fudge factor to reduce outside areas

if ~exist('verb','var'),       verb = []; end
if isempty(verb),              verb = 1; end
if ~exist('out_corr','var'),   out_corr = []; end
if isempty(out_corr),          out_corr = 2; end
if ~exist('dcf_thresh','var'), dcf_thresh = []; end


%% input checks
if size(k,1)~=1, error('size(k,1)~=1'); end
if isreal(k)
    if length(size(k))==2
        error('isreal(k) && length(size(k))(=%g)==2',length(size(k)));
    end
else
    if length(size(k))>2
        error('~isreal(k) && length(size(k))(=%g)>2',length(size(k)));
    end
    tmp = complex(zeros(1,size(k,2),2,class(k)));
    tmp(1,:,1) = real(k);
    tmp(1,:,2) = imag(k);
    k = tmp;
    clear tmp;
end
k = squeeze(k);         % size(k)=(#kpts,dim)
nk = size(k,1);


%% preprocess data
[ku,k_ind1,k_ind2] = unique(k,'rows','stable');

nki1 = size(k_ind1,1);                 % #k-points used for voronoin
nupts = nk-nki1;                       % #non-unique, excluded points
if ((nupts>0) && (verb>0))
    fprintf('%g non-unique points\n',nupts);
end

%% make an index list of duplicate points
lnu = 0;
ind = {};
ii = find(accumarray(k_ind2,1)>1);     % preselect repeated data
for lki1=ii.'
    tmp = k_ind2==lki1;
    if sum(tmp,1)>1
        lnu = lnu+1;
        ind{lnu} = find(tmp);
    end
end



%% calculate voronoi vertices
[V,C] = voronoin(ku);
dcf  = zeros(1,size(ku,1));
maxk = max(abs(ku),[],1);
lout = 0;


%% loop through all triangles and determine area
for i1 = 1:length(C)
    x = V(C{i1},:);
    
    if size(x,1)>0
        
        if out_corr==4
            % replace edge points (inf) with k-space
            ii = any(isinf(x),2) | any(abs(x)>maxk,2);
            if any(ii)
                x(ii,:) = repmat(ku(i1,:),[sum(ii),1]);
            end
        end
        
        % check for points outside the sampled area
        % if  any(isinf(x(:))) || any(any(abs(x)>maxk,2))
        if  any(isinf(x(:))) || any(any(abs(x)>repmat(maxk,[size(x,1) 1]),2))
            % outside: radius -> approximate area
            lout = lout+1;
            switch out_corr
                case 1       % approximate via radius
                    radii = sort(sqrt(sum(bsxfun(@minus,x,ku(i1,:)).^2,2)));
                    dcf(1,i1) = radii(1,1)*radii(2,1)*fufa;
                    if dcf(1,i1)>1
                        warning('dcf(1,i1)>1');
                    end
                case 2       % use last one
                    if i1>1
                        dcf(1,i1) = dcf(1,i1-1);
                    else
                        dcf(1,i1) = NaN;
                    end
                case {3,4}   % use nearest neighbour
                    dcf(1,i1) = NaN;
                otherwise % do nothing, i.e. dcf=0
            end
        else
            % inside: correct area
            % dcf(1,i1) = polyarea(x(:,1),x(:,2));    % 2D only
            [~,dcf(1,i1)] = convhull(x);
        end
    else
        warning('size(x,1)(=%g)==0',size(x,1));
        dcf(1,i1) = NaN;      % replace later with nearest neighbour
    end
end


%% checks
if (lout>0) && (verb>0), fprintf('points outside = %g\n',lout); end
if any(isinf(dcf(:)))
    dcf(isinf(dcf)) = 0;
    warning('dcf contains inf -> setting to zero');
end
if any(isnan(dcf(:)))        % replace NaNs with nearest neighbour
    ii = find(isnan(dcf));
    kurem = ku(~isnan(dcf),:);
    dcfrem = dcf(1,~isnan(dcf));
    
    for ln=1:length(ii)
        [~,kri] = min(sum(bsxfun(@minus,kurem,ku(ii(ln),:)).^2,2));
        dcf(1,ii(ln))=dcfrem(1,kri);
    end
    

end
if any(isnan(dcf))
    warning('dcf contains NaN -> setting to zero');
    dcf(isnan(dcf)) = 0;
end


%% include non-unique points in final dcf
dcf_out = zeros(1,nk);
dcf_out(1,k_ind1) = dcf;
for li=1:length(ind)
    ii = ind{li};                      % index list to non-unique points
    ni = length(ii);                   % #non-unique points
    val = dcf_out(1,ii(1))/ni;         % scale 1st non-unique point
    for lii=1:ni, dcf_out(1,ii(lii)) = val; end
end


%% check for large dcf values
if ~isempty(dcf_thresh)
    ii = dcf_out>dcf_thresh;
    if verb>0
        fprintf('Capping %g dcf values at dcf_thresh=%g\n',sum(ii),dcf_thresh);
    end
    dcf_out(1,ii) = dcf_thresh;
end
if max(dcf_out)>2
    warning('max(dcf_out)(=%g)>2',max(dcf_out));
end
mean_dcf_unsamp = mean(dcf_out(dcf_out>1));
if mean_dcf_unsamp>1
    warning('mean(dcf(dcf>1))=%g',mean_dcf_unsamp);
end


%% plotting
if verb>1, clf,plot(dcf_out); end


end      % main function voronoi_area.m
