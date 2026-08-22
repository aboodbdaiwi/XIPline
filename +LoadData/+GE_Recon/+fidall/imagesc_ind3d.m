function bplt = imagesc_ind3d(b,ind,cm,zm,plt_cross,tght,scale,rot,cb,...
    do_repelem)
%IMAGESC_IND3D Plot slices through 3D images in 3 orientations
%  bplt = imagesc_ind3d(b,ind,cm,zm,plt_cross,tght,scale,rot,cb,do_repelem)
%
%        b   4D image matrix to display (x,y,z,#rows1)
%            #rows1 = 1 or nrows
%      ind   Index to slices+crosses    (#rows2,[ix iy iz])         (mtx/2)
%            #rows2 = 1 or nrows
%       cm   Color axis scaling                                 ([min max])
%       zm   Zoom: pixels around centre to plot                        ([])
%plt_cross   Plot cross of index                                    (false)
%     tght   Minimise white space around figure                      (true)
%    scale   individual scaling imagesc_row      ([]->off)             ([])
%      rot   Rotate slices via rot90 [x,y,z]                           ([])
%       cb   Display colorbar                                       (false)
%do_repelem  Fill full image squares (if mtx not isotropic)          (true)
%
%     bplt  Cross-sections to plot      (max(x,y,z),max(x,y,z),3)
%
% See also IMAGESC_ROW.
%  5/2024 Rolf Schulte
if (nargin<1), help(mfilename); return; end


%% default parameters + checks
if ~exist('ind','var'),  ind = []; end
if ~exist('cm','var'),   cm = []; end
if ~exist('zm','var'),   zm = []; end
if ~exist('plt_cross','var'), plt_cross = []; end
if isempty(plt_cross),   plt_cross = false; end
if ~exist('tght','var'), tght = []; end
if isempty(tght),        tght = true; end
if ~exist('scale','var'),scale = []; end
if ~exist('rot','var'),  rot = []; end
if ~exist('cb','var'),   cb = []; end
if ~exist('do_repelem','var'), do_repelem = []; end
if isempty(do_repelem),  do_repelem = true; end

if ~isreal(b(:)),     fprintf('b not real: taking abs\n'); b = abs(b); end
if any(isnan(b(:))),  warning('b isnan'); end
if any(isinf(b(:))),  warning('b isinf'); end
if length(size(b))>4, warning('b>4D'); end
if length(size(b))<3, warning('b<3D'); end


%% misc parameters
[m1,m2,m3,m4] = size(b);          % m1,m2,m3=image matrix; m4=nrows or 1
if isempty(ind), ind = floor([m1 m2 m3]/2+0.5); end   % default=centre
mm = max([m1,m2,m3]);             % display (zerofill to) max of all 3 dims
if size(ind,2)~=3, warning('size(ind,2)(=%g)~=3',size(ind,2)); end
nind1 = size(ind,1);              % 1 or nrows


%% calculate nrows (max([m4 nind1])) and repmat b & ind
if m4~=nind1
    if ((m4==1)&&(nind1>1))
        nrows = nind1;
        b = repmat(b,[1 1 1 nrows]);
    end
    if ((m4>1)&&(nind1==1))
        nrows = m4;
        ind = repmat(ind,[nrows 1]);
    end
    if ((m4>1)&&(nind1>1))
        warning('size(b,4)(%g)~=size(ind,1)(=%g)',m4,nind1); 
    end
else
    nrows = m4;
end
if ~exist('nrows','var'), error('nrows not existing'); end


%% construct matrix for imagesc
ind_round = floor(ind+0.5);
if isempty(zm)
    if do_repelem
        rr = floor(max([m1 m2 m3])./[m1 m2 m3]);
    else
        rr = [1 1 1];
    end
    i1 = 1:(rr(1)*m1); 
    i2 = 1:(rr(2)*m2); 
    i3 = 1:(rr(3)*m3);
    bplt = zeros(mm,mm,3,nrows);           % zerofill for differing mtx sizes
    for l4=1:nrows
        bplt(i1,i2,1,l4) = repelem(squeeze(b(:,:,ind_round(l4,3),l4)),rr(1),rr(2));
        bplt(i1,i3,2,l4) = repelem(squeeze(b(:,ind_round(l4,2),:,l4)),rr(1),rr(3));
        bplt(i2,i3,3,l4) = repelem(squeeze(b(ind_round(l4,1),:,:,l4)),rr(2),rr(3));
    end
else
    %% zoom in
    if length(zm)~=1, error('length(zm)~=1'); end
    mm = 2*zm+1;
    bplt = zeros(mm,mm,3,nrows);
    for l4=1:nrows
        for l=1:3, ii{l} = (-zm:zm) + ind_round(l4,l); end
        bplt(:,:,1,l4) = b(ii{1},ii{2},ind_round(l4,3),l4);
        bplt(:,:,2,l4) = b(ii{1},ind_round(l4,2),ii{3},l4);
        bplt(:,:,3,l4) = b(ind_round(l4,1),ii{2},ii{3},l4);
    end
    ind = ind-ind_round+zm+1; 
end


%% rot90
if ~isempty(rot)
    if length(rot)==3
        for l4=1:nrows
            for l3=1:3
                bplt(:,:,l3,l4) = rot90(bplt(:,:,l3,l4),rot(l3));
            end
        end
    else
        warning('length(rot)(=%d)~==3; skipping rotation',length(rot));
    end
end


%% actual plotting
LoadData.GE_Recon.fidall.imagesc_row(bplt,cm,scale,'',tght,[],[],cb);


%% plot crosses of index in all slices
if plt_cross
    hold on;
    pln = {[1 2],[1 3],[2 3]};
    for l4=1:nrows
        for l3=1:3
            plt_ind = ind(l4,pln{l3}) + ...
                [(l4-1)*mm (l3-1)*mm];
            plot(plt_ind(2),plt_ind(1),...
                'xy','MarkerSize',6,'LineWidth',1);
        end
    end
end

if nargout==0, clear bplt; end

end      % end main function (imagesc_ind3d.m)
