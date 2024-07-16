function [dd] = epsi_reshape(d,npts,ind,epsi_intlv,epsidir,phafu,...
    cor_intlv,epsi_mode,verb)
%EPSI_RESHAPE  Reshaping EPSI data
%     [dd] = epsi_reshape(d,npts,ind,epsi_intlv,epsidir,phafu,...
%                         cor_intlv,epsi_mode,verb)
%         d  Raw data                     [#acqs,#samples,1,#t,#slices,#c]
%      npts  Reshaping matrix                          [#s,#x,#y,#z,#t,#c]
%            for multi-slices #z=1 and size(d,5)>1
%       ind  Index list for selecting sampled data (dim=2)
%epsi_intlv  #EPSI interleaves
%   epsidir  Direction of spatial EPSI encoding   (1=x,2=y,3=z)
%     phafu  Phase function (for epsi_iso==1)                   default=[]
% cor_intlv  Correct phase of interleaves (via svds)            def=false
% epsi_mode  1=flyback;2=bidirectional                          default=1
%      verb  Verbose mode (0=off;1=verbose;2=verb+plotting)     default=0
%
%        dd  Reshaped EPSI data                   [#s,#x,#y,#z,#t,#c]
%            multi-slice in #z
%
% 1/2018 Rolf Schulte
if nargin<1, help(mfilename); return; end
if ~exist('phafu','var'),     phafu = []; end
if ~exist('cor_intlv','var'), cor_intlv = []; end
if isempty(cor_intlv),        cor_intlv = false; end
if ~exist('epsi_mode','var'), epsi_mode = []; end
if isempty(epsi_mode),        epsi_mode = 1; end
if ~exist('verb','var'),      verb = []; end
if isempty(verb),             verb = 0; end


%% misc input + pars
[nacq,ns,n3,n4,n5,nc] = size(d);
if n3~=1, warning('size(d,3)(=%g) ~= 1',n3); end
if (n5>1)
    if (size(d,5)~=npts(4)), 
        warning('size(d,5)(=%g) ~= npts(4)(=%g)',n5,npts(4)); 
        fprintf('   setting npts(4)=%g\n',n5);
    end
    nslices = n5;
    nz = 1;
else
    nz = npts(4);
    nslices = 1;
end
 
if ~any(epsidir==(1:3)), error('epsidir(=%g) not 1,2 or 3',epsidir); end
npts = npts(:).';
npts = [npts ones(1,5-length(npts))];
if length(npts)<6, error('length(npts)<6'); end
nk = sum(ind);
% if nk~=npts(1)*npts(epsidir+1)*npts(5)/epsi_intlv, 
%     warning('sum(ind)(=%g) ~= npts(1)(=%g)*npts(%g)(=%g)*npts(5)(=%g)/epsi_intlv(=%g)=%g',...
%         nk,npts(1),epsidir+1,npts(epsidir+1),npts(5),epsi_intlv,...
%         npts(1)*npts(epsidir+1)*npts(5)/epsi_intlv); 
%     npts(1) = npts(1)+1;
% end
trunc = false;
if nk>npts(1)*npts(epsidir+1)/epsi_intlv,
    if nk<(npts(1)+1)*npts(epsidir+1)/epsi_intlv,
        warning('sum(ind)(=%g) ~= npts(1)(=%g)*npts(%g)(=%g)/epsi_intlv(=%g)=%g',...
            nk,npts(1),epsidir+1,npts(epsidir+1),epsi_intlv,...
            npts(1)*npts(epsidir+1)/epsi_intlv); 
    else
        trunc = true;
        npts1_orig = npts(1);
        tmp = nk*epsi_intlv/npts(epsidir+1);
        if abs(tmp-floor(tmp))>1d-10, 
            warning('rounding ns (%g->%g)',tmp,floor(tmp+0.5));
            tmp = floor(tmp+0.5);
        end
        if verb>0, fprintf('npts(1)=%g -> %g\n',npts(1),tmp); end
        npts(1) = tmp;
    end
end


%% reshape data w/o epsi
nn = npts;
nn(1) = nk*epsi_intlv; nn(epsidir+1) = 1;
ned = npts(epsidir+1);
d1 = zeros(nn);
if verb>0, 
    fprintf('nn =   '); disp(nn);
    fprintf('npts = '); disp(npts)
end
ii = false(epsi_intlv,nn(1));
for li=1:epsi_intlv,
    for lspec = 1:(npts(1)/epsi_intlv),
        ii(li,((li-1)*ned+(1:ned)+epsi_intlv*(lspec-1)*ned)) = true;
    end
end
for l2=1:nn(2),
    for l3=1:nn(3),
        for l4=1:nz, % nn(4),
            for li=1:epsi_intlv,
                ll1 = l2+(l3-1)*nn(2)+(l4-1)*nn(3)*nn(2);
                ll1 = epsi_intlv*(ll1-1)+li;
                if verb>0, fprintf('l2=%g l3=%g l4=%g li=%g ll1=%g\n',l2,l3,l4,li,ll1); end
                for l5=1:nn(5),
                    for l6=1:nn(6),
                        for lsli=1:nslices,
                            d1(ii(li,:),l2,l3,l4+lsli-1,l5,l6) = d(ll1,ind,1,l5,lsli,l6).';
                            % d1(ii(li,end:-1:1),l2,l3,l4,l5,l6) = d(ll1,ind,1,l5,1,l6).';
                        end
                    end
                end
            end
        end
    end
end


%% compensate phase errors (off-isocentre; see header2epsi.m)
if ~isempty(phafu),
    if verb>0, fprintf('Applying phafu\n'); end
    if size(phafu,1)~=nn(1), 
        warning('size(phafu,1)(=%g) ~= nn(1)(=%g)',size(phafu,1),nn(1));
    end
    if size(phafu,2)~=1,
        warning('size(phafu,1)(=%g) ~= 1',size(phafu,2));
    end
    if length(size(phafu))~=2, warning('length(size(phafu))~=2'); end
    if any(abs(abs(phafu)-1)>1d-7), warning('abs(phafu) ~= 1'); end
    
    d1 = bsxfun(@times,d1,phafu);
end


%% reshape data: epsi part
dd = zeros(npts);
for l_spec=1:npts(1),
    for l_spat=1:npts(epsidir+1),
        % l = (l_spec-1)*npts(epsidir+1)+l_spat;
        l = (l_spec-1)*npts(epsidir+1)+npts(epsidir+1)-l_spat+1;
        % fprintf('l2=%g; l1=%g, l=%g\n',l2,l1,l);
        %         dd(:,l2,l1,:,:,:) = d(:,l,1,:,:,:);
        switch epsidir
            case 1, dd(l_spec,l_spat,:,:,:,:) = d1(l,1,:,:,:,:);
            case 2, dd(l_spec,:,l_spat,:,:,:) = d1(l,:,1,:,:,:);
            case 3, dd(l_spec,:,:,l_spat,:,:) = d1(l,:,:,1,:,:);
        end
    end
end


%% reverse direction for bidirectional EPSI: every other spec pts
if epsi_mode>1
    switch epsidir
        case 1, dd(2:2:end,:,:,:,:,:) = dd(2:2:end,end:-1:1,:,:,:,:);
        case 2, dd(2:2:end,:,:,:,:,:) = dd(2:2:end,:,end:-1:1,:,:,:);
        case 3, dd(2:2:end,:,:,:,:,:) = dd(2:2:end,:,:,end:-1:1,:,:);
    end
end


%% correct phase difference between interleaves
if ((cor_intlv) && (epsi_intlv>1)),
    if nslices>1,
        for l4=1:nslices,
            tmp = abs(dd(:,:,:,l4,:,:));
            [maxval,maxind] = max(tmp(:));
            [i1,i2,i3,i4,i5,i6] = ind2sub(size(tmp),maxind);
        
            mind_svds = min([50 npts(1)]);
            mind_svds = floor(mind_svds/epsi_intlv)*epsi_intlv;
            if i1>mind_svds,
                warning('max(abs(dd),[],1) index=%g > mind_svds=%g',...
                    i1,mind_svds);
            end
            A = zeros(mind_svds/epsi_intlv,epsi_intlv);
            for l=1:epsi_intlv,
                A(:,l) = dd(l:epsi_intlv:mind_svds,i2,i3,l4,i5,i6);
            end
            % [U,S,V,flag] = svds(A,1);
            % phi = angle(V);
            % fix for bug in matlab R2015b svds function
            [U,S,V,flag] = svds(A,2);
            phi = angle(V(:,1));
            phi = exp(1i*(phi-phi(1)));
            phi = repmat(phi,[npts(1)/epsi_intlv 1]);
            dd(:,:,:,l4,:,:) = bsxfun(@times,dd(:,:,:,l4,:,:),phi);
        end
    else
        [maxval,maxind] = max(abs(dd(:)));
        [i1,i2,i3,i4,i5,i6] = ind2sub(size(dd),maxind);
        
        mind_svds = min([50 npts(1)]);
        mind_svds = floor(mind_svds/epsi_intlv)*epsi_intlv;
        if i1>mind_svds,
            warning('max(abs(dd),[],1) index=%g > mind_svds=%g',...
                i1,mind_svds);
        end
        A = zeros(mind_svds/epsi_intlv,epsi_intlv);
        for l=1:epsi_intlv,
            A(:,l) = dd(l:epsi_intlv:mind_svds,i2,i3,i4,i5,i6);
        end
        % [U,S,V,flag] = svds(A,1);
        % phi = angle(V);
        % fix for bug in matlab R2015b svds function
        [U,S,V,flag] = svds(A,2);
        phi = angle(V(:,1));
        phi = exp(1i*(phi-phi(1)));
        phi = repmat(phi,[npts(1)/epsi_intlv 1]);
        dd = bsxfun(@times,dd,phi);
    end
end


%% truncate spectral dimension -> compensate rounding error from epsi_intlv
if trunc
    npts(1) = npts1_orig;
    dd = dd(1:npts1_orig,:,:,:,:,:);
end


%% check size consistency
for l=1:6,
    if size(dd,l)~=npts(l), 
        warning('size mismatch: size(dd,%g)(=%g) ~= npts(%g)(=%g)',...
            l,size(dd,l),l,npts(l)); 
    end
end


%% plotting re-ordered fid for debugging
if verb>1,
    tt=(0:size(d,2)-1);
    ii = size(d,1)/2+(1:epsi_intlv);
    figure(74); clf;
    
    subplot(3,1,1);
    plot(tt,abs(d(ii,:,1,1,1,1)).','',tt(ind),abs(d(ii,ind,1,1,1,1)).','x'); 
    ylabel('abs(d)');
    
    subplot(3,1,2);
    plot(tt,real(d(ii,:,1,1,1,1)).',''); 
    ylabel('real(d)');
    
    subplot(3,1,3);
    t = (0:size(dd,1)-1);
    ix = floor(npts(2)/2)+1;
    iy = floor(npts(3)/2)+1;
    plot(t,abs(dd(:,ix,iy,1,1,1)),'b',...
        t,real(dd(:,ix,iy,1,1,1)),'g',...)
        t,imag(dd(:,ix,iy,1,1,1)),'r');
    ylabel('dd');
    drawnow;
end

