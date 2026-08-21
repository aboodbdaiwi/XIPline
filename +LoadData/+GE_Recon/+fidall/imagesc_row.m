function imagesc_row(b,cm,scale,rshp,tght,do_clf,possize,cb)
%IMAGESC_ROW  Multiple 2D imagesc in rows: dim1+2=imagesc; dim3=rows
% imagesc_row(b,cm,scale,rshp,tght,do_clf,possize,cb)
%        b   Real-valued images (n1,n2,n3,n4)
%       cm   colormap (min max)                  (min+max of image) -> ([])
%    scale   individual scaling (n3,n4)                (no scaling) -> ([])
%            'row' scale rows
%            'col' scale columns
%            'ind' scale subimages individually to max
%            'inm' scale individually to mean of top 3% values
%            'gl3' scale globally to mean of top 3% values
%     rshp   Reshape data to 4D for best display                    (false)
%     tght   Minimise white space around figure                      (true)
%   do_clf   Do clf                                                  (true)
%  possize   Size of figure [pts]                                     (800)
%       cb   Display colorbar                                       (false)
%
%  2/2025 Rolf Schulte
if (nargin<1), help(mfilename); return; end


off2 = 25;              % offset from lower boundary
bound2 = 70;            % outer boundary
inm_fac = 0.97;         % scale=='inm': mean of 1-inm_fac times inm_fac


%% misc parameters + checks
if ~isreal(b),       warning('b must be real; taking abs'); b = abs(b); end
if any(isnan(b(:))), warning('b isnan; zeroing'); b(isnan(b)) = 0; end
if any(isinf(b(:))), warning('b isinf; zeroing'); b(isinf(b)) = 0; end
if isempty(b),       warning('isempty(b): returning'); return; end

if ~exist('cm','var'),    cm = []; end
if any(isinf(cm(:))),warning('cm isinf; setting to empty'); cm = []; end
if ~exist('scale','var'), scale = []; end
if ~exist('rshp','var'),  rshp = []; end
if isempty(rshp),         rshp = false; end
if ~exist('tght','var'),  tght = []; end
if isempty(tght),         tght = true; end
if ~exist('do_clf','var'),do_clf = []; end
if isempty(do_clf),       do_clf = true; end
if do_clf, clf; end
if ~exist('possize','var'),possize = []; end
if isempty(possize),      possize = 800; end
if ~exist('cb','var'),    cb = []; end
if isempty(cb),           cb = false; end


%% reshape data
if rshp
    if length(size(b))~=3
        fprintf('imagesc_row: length(size(b))=%g\n',length(size(b)));
    end
    [n1,n2,n3] = size(b);
    nn3 = ceil(sqrt(n3));
    nn4 = ceil(n3/nn3);
    b = reshape(b,[n1 n2 n3]);
    b1 = zeros(n1,n2,nn3,nn4);
    for l4=1:nn4
        for l3=1:nn3
            ll = l3 + (l4-1)*nn3;
            if ll<=n3, b1(:,:,l3,l4) = b(:,:,ll); end
        end
    end
    b = b1;
end


%% check length
n = size(b);
if (length(n)~=2 && length(n)~=3 && length(n)~=4)
    warning('length(size(b))=%g; b should be 2D, 3D or 4D',length(n));
    fprintf('trying squeeze\n');
    b = squeeze(b);
    n = size(b);
    if length(n)>4, error('squeeze(b) failed'); end
end
[n1,n2,n3,n4] = size(b);


%% scaling
if max(abs(b(:)))==0
    warning('max(abs(b(:)))==0');
    cm = [0 1];
end
if strcmpi(scale,'row')
    scale = 1./(ones(n3,1)*squeeze(max(max(max(b,[],1),[],2),[],3)).'); 
end
if strcmpi(scale,'col')
    scale = 1./(squeeze(max(max(max(b,[],1),[],2),[],4))*ones(1,n4)); 
end
if strcmpi(scale,'ind')
    scale = reshape(1./(max(max(b,[],1),[],2)),[n3 n4]);
end
if strcmpi(scale,'inm')
    scale = ones(n3,n4);
    for l3=1:n3
        for l4=1:n4
            bs = abs(b(:,:,l3,l4));
            bs = sort(bs(:));
            ii = round(n1*n2*inm_fac):(n1*n2);
            scale(l3,l4) = inm_fac./mean(bs(ii));
        end
    end
    if isempty(cm)
        if any(b(:)<0)
            cm = [-1 1];
        else
            cm = [0 1];
        end
    end
end
if strcmpi(scale,'gl3')
    bs = sort(abs(b(:)),'descend');
    mm = 0.97*mean(bs(1:ceil(0.03*length(bs))));
    scale = ones(n3,n4)/mm;
    if isempty(cm)
        if any(b(:)<0)
            cm = [-1 1];
        else
            cm = [0 1];
        end
    end
end

if ~isempty(scale)
    ns = size(scale);
    if ns(1)~=n3, error('size(scale,1)=%g ~= size(b,3)=%g',ns(1),n3); end
    if ns(2)~=n4, error('size(scale,2)=%g ~= size(b,4)=%g',ns(2),n4); end
end


%% concatenate + scale data
% apply scaling directly to images -> colorbar only valid for global scaling
bb = zeros(n1*n4,n2*n3);
for l3=1:n3
    for l4=1:n4
        ind1 = (1:n1) + (l4-1)*n1;
        ind2 = (1:n2) + (l3-1)*n2;
        if isempty(scale)
            bb(ind1,ind2) = b(:,:,l3,l4);
        else
            bb(ind1,ind2) = scale(l3,l4)*b(:,:,l3,l4);
        end
    end
end


%% default color axis
if isempty(cm)
    % val = sort(abs(bb(:)));
    % max_val = min([1.1*val(floor(0.99*length(val))),max(val)]);
    max_val = max(abs(bb(:)));
    if any(b(:)<0)
        cm = max_val*[-1 1];
    else
        cm = [0 max_val];
    end
end

if any(isnan(cm))
    warning('cm contains NaN; setting to 1');
    cm(isnan(cm)) = 1;
end


%% actual plotting
% concatenated data via imagesc; boundaries added by plot
imagesc(bb);
hold on
for l3=(0:n3)+1/n2/2, plot(l3*n2*[1 1],[0 1]*n1*n4+0.5,'w'); end
for l4=(0:n4)+1/n1/2, plot([0 1]*n2*n3+0.5,l4*n1*[1 1],'w'); end
hold off
caxis(cm);
axis image
set(gca,'XTick',-1,'YTick',-1,'XColor',[1 1 1],'YColor',[1 1 1]);
xlabel('','Color',[0 0 0])
ylabel('','Color',[0 0 0])
box off

if tght
    set(gca,'Position',[0 0 1 1]);          % minimise white space
    posf = get(gcf,'position');             % get figure position
    ss = get(0,'ScreenSize');
    if strcmpi(get(gcf,'Visible'),'on')
        scsi = get(0,'ScreenSize');
        if (possize/n3*n4)<100, possize = round(100*n3/n4); end
        if possize>scsi(3), possize=round(0.9*scsi(3)); end
        if (possize/n3*n4)>(scsi(4)-off2-bound2)
            possize = 0.9*(scsi(4)-off2-bound2)*n3/n4;
        end
    end
    possize2 = possize/n3*n4*n1/n2;
    if (posf(2)+possize2+100)>ss(4), posf(2) = ss(4)-(possize2+100); end
    % set(gcf,'position',[posf(1) off2 possize possize/n3*n4]);
    set(gcf,'position',[posf(1:2) possize possize2]);
    
    
    set(gcf,'PaperUnits','points','PaperPositionMode','manual');
    ppos = get(gcf,'PaperPosition');
    set(gcf,'PaperPosition',[ppos(1) ppos(2) possize possize/n3*n4]);
end
if cb
    colorbar;
    if tght, set(gca,'Position',[0 0 0.92 1]); end
end

end      % imagesc_row.m
