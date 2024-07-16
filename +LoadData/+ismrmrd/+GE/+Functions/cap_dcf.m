function dcf = cap_dcf(dcf,k,mtx,rad_cap)
%CAP_DCF  Cap values of dcf at/to certain radius
%    dcf = cap_dcf(dcf,k,rad)
%    dcf   Density compensation function
%      k   k-Space normalised to +-0.5 (1,#k-pts,dim)
%    mtx   Matrix size
%rad_cap   Radius at which to cap dcf
%
%  7/2019 Rolf Schulte
if (nargin<1), help(mfilename); return; end


%% input parameters
if ~exist('rad_cap','var'), rad_cap = []; end
if isempty(rad_cap),        rad_cap = 0.5-1.1/mtx; end

dpre = 2/mtx;                          % distance for pre-selection
infac = 1-2/mtx;                       % factor for inner radius
dim = size(k,3);                       % k-space dimensions
nk = size(k,2);                        % #k-pts


%% checks
if ~isreal(k),     error('~isreal(k)'); end
if ~isreal(dcf),   error('~isreal(dcf)'); end
if size(k,1)~=1,   error('size(k,1)(=%g)~=1',size(k,1)); end
if size(dcf,1)~=1, error('size(dcf,1)(=%g)~=1',size(dcf,1)); end
if size(dcf,2)~=nk
    error('size(dcf,2)(=%g)~=size(k,2)(=%g)',size(dcf,2),nk); 
end
if ((dim<2)||(dim>3)), error('size(k,3)(=%g)~=2 or 3',dim); end
if max(abs(k(:)))>0.5
    warning('max(abs(k(:)))(=%g)>0.5',max(abs(k(:)))); 
end
if rad_cap>0.5, warning('rad_cap(=%g)>0.5',rad_cap); end
if rad_cap<0,   warning('rad_cap(=%g)<0',rad_cap); end
if length(rad_cap)~=1
    error('length(rad_cap)(=%g)~=1',length(rad_cap)); 
end


%% cap data: search for nearest points inside rad_cap; set dcf to it
rad = sqrt(sum(k.^2,3));               % radius
icap = rad>rad_cap;                    % index for capping
irem = ~icap & (rad>infac*rad_cap);    % index to remaining k+dcf
dcfrem = dcf(1,irem);                  % reduce data size
k = single(k); dpre = single(dpre);    % slight speed improvement
krem = k(1,irem,:);                    % reduce data size

% % numerically too slow for large data
% for lc=find(icap)
%     % [~,ii] = min(sqrt(sum(bsxfun(@minus,krem,k(1,lc,:)).^2,3)));
%     [~,ii] = min(sum(bsxfun(@minus,krem,k(1,lc,:)).^2,3));
%     
%     dcf(1,lc) = dcfrem(1,ii);
% end


%% loop through outside data points
for lc=find(icap)
    % pre-select data surrounding point
    switch dim
        case 2
            ipre = find((abs(krem(1,:,1)-k(1,lc,1))<dpre) & ...
                (abs(krem(1,:,2)-k(1,lc,2))<dpre));
        case 3 
            ipre = find((abs(krem(1,:,1)-k(1,lc,1))<dpre) & ...
                (abs(krem(1,:,2)-k(1,lc,2))<dpre) & ...
                (abs(krem(1,:,3)-k(1,lc,3))<dpre));
    end
    
    % search for closest data point
    % [~,ii] = min(sqrt(sum(bsxfun(@minus,krem,k(1,lc,:)).^2,3)));
    [~,ii] = min(sum(bsxfun(@minus,krem(1,ipre,:),k(1,lc,:)).^2,3));
    if isempty(ipre)
        warning('isempty(ipre); lc=%g',lc); 
    end
    if isempty(ii),   warning('isempty(ii); lc=%g',lc);   end
    
    % set to value of closest data point
    dcf(1,lc) = dcfrem(1,ipre(ii));
end


end      % main function cap_dcf.m
