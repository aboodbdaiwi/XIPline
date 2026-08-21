function [a_arr,lb,freq,fid_sup,fid_fit,ind_arr] = hsvd(fid_arr,p,bound,bw,verb,nmax)
%HSVD  Hankel singular value decomposition.
% Decompose a time domain signal into decaying sinusoids (ie Lorentzians).
% The singular matrix yields decay constants and frequencies. A subsequent
% linear least squares fit yields (complex) amplitudes, ie concentrations
% and phases of peaks. HSVD uses the svds function (Lanczos algorithm too
% imprecise for noisy data).
%
%[alpha,lb,freq,fid_sup,fid_fit,ind] = hsvd(fid,p,bound,bw,verb,nmax)
%                                                               (default)
%      fid  time-domain data: spectral dim/hsvd along dim=1
%        p  number of decaying sinusoids                        (25)
%    bound  boundaries for fitting (normalised to bw)           ([-1 1]/2)
%       bw  bandwidth                                           (1)
%     verb  verbose: 0=off; 1=limited; 2=extensive; 3=extensive+plotting)
%                                                               (0)
%     nmax  maximum matrix size for SVD     (default ~ 5*T2*)
%
%    alpha  complex concentrations
%       lb  decays                          (normalised to bw)
%     freq  frequencies                     (normalised to bw)
%  fid_sup  suppressed FID (fid-fid_fit)
%  fid_fit  generated FID of fit (in time domain for truncation)
%      ind  index list: lb>0 and freq inside bound
%
%  Literature:
%    1) H Barkhuijsen, et al. JMR 73, 553-557 (1987)
%    2) "NMR Data Processing", JC Hoch and AS Stern; p89-92.
%
% 11/2020 Rolf Schulte
% See also SVDS.
if nargin<1, help(mfilename); return; end

plw = [1 0.5];   % linewidths for plotting


%% inputs
si = size(fid_arr);
[n1,n2] = size(fid_arr);                         % if >2D -> loop over n2
if ~exist('p','var'),     p = []; end            % #decaying sinusoids
if isempty(p),            p = 25; end
if ~exist('bound','var'), bound = []; end        % boundary for fitting
if ~exist('bw','var'),    bw = [];    end        % bandwidth
if isempty(bw),           bw = 1; end
if isempty(bound),        bound = [-0.5 0.5]*bw; end
if ~exist('verb','var'),  verb = [];    end      % verbose mode
if isempty(verb),         verb = false; end
if ~exist('nmax','var'),  nmax = []; end         % #data points
if verb>2     % plotting
    plt = true; 
else
    plt = false; 
end


%% misc checks
if si(1)==1,         error('size(spec,1)==1'); end
if length(bw)~=1,    error('bw not scalar value'); end
if length(bound)~=2, error('length(bound)~=2'); end
if bound(1)>bound(2)
    warning('bound(1)>bound(2): reordering\n');
    bound = bound(end:-1:1);
end
bound = bound/bw;                             % normalise to +/-0.5


%% determine nmax by looking where FID decayed
fida = sum(abs(reshape(fid_arr,[n1 n2])),2);
if isempty(nmax)
    nnend = 20;
    thresh = mean(fida(n1-(0:nnend-1))) + 5*std(fida(n1-(0:nnend-1)));
    nmax = ceil(1.3*find(fida>thresh,1,'last'));
    max_nmax = false;
    if nmax<128
        warning('nmax(=%g)<128',nmax);
        max_nmax = true;
    end
    if 10*mean(fida(n1-(0:nnend-1))) > mean(fida(1:5))
        warning('10*mean(=%g) > mean(fid(1:5))(=%g)',...
            10*mean(fida(n1-(0:nnend-1))),mean(fida(1:5)));
        max_nmax = true;
    end
    if mean(fida(n1-(0:nnend-1))) > 5*std(fida(n1-(0:nnend-1)))
        warning('mean(fida(n1-(0:nnend-1)))(=%g) > 5*std(fida(n1-(0:nnend-1)))(=%g)',...
            mean(fida(n1-(0:nnend-1))),5*std(fida(n1-(0:nnend-1))));
        max_nmax = true;
    end
    if max_nmax
        fprintf('Warning: setting nmax(=%g) to max(=%g)\n',nmax,n1);
        nmax = n1;
        plot((0:n1-1),fida,'b',[0 n1-1],thresh*[1 1],'r',nmax,0,'mx');
        drawnow; pause(1);
    end
end

tmp = diff(fida(1:20));
if tmp(1)>10*abs(mean(tmp(2:end)))
    warning('1st point of FID significantly smaller; excluding');
    tmp = zeros(size(fid_arr));
    tmp(1:end-1,:) = fid_arr(2:end,:);
    fid_arr = tmp;
end


%% misc parameters/initialisations
n = min([nmax n1]);
m = floor(2*n/3);            % filter order
l = n-m+1;                     % #equations
if l<p
    warning('l<p: p(=%g) too big; n(=%g) too small: decreasing p to %g',p,n,l-1);
    p = l-1;
end
if verb>0, fprintf('n=%g; m=%g; l=%g\n',n,m,l); end
fid_fit = zeros(si);
a_arr   = zeros([p,si(2:end)]);
lb      = zeros([p,si(2:end)]);
freq    = zeros([p,si(2:end)]);
ind_arr = false([p,si(2:end)]);

if plt
    clf;
    x_ax = (-n1/2:(n1/2-1))/n1*bw;
    spec = ifftshift(fft(fid_arr,[],1),1);
    t = (0:n1-1)/bw;
end

n_pts = floor(diff(bound)*n);
if n_pts<5, warning('#pts within bound region=%g\n',n_pts); end


%% loop through FIDs with actual HSVD
for ll=1:n2
    fid = fid_arr(:,ll);
    
    % construct data matrix
    D = repmat((1:m),[l 1])+repmat((1:l).',[1 m])-1;
    D = fid(D);
    
    if verb>1
        if n2>1, fprintf('%g/%g: \t',ll,n2); end
        % fprintf('cond(D) = %g; rank(D) = %g\n',cond(D),rank(D));
    end
    
    % svd for only the p strongest signals
    [U,S] = eig(D*D'); U = U(:,[l:-1:l-p+1]);
    % [U,S] = eigs(D*D',p);                      % 2nd most expensive
    % [U,S,V] = svd(D); U = U(:,1:p);            % 3rd most expensive
    % [U,S,V,flag] = svds(D,p);                  % most expensive
    
    if size(U,2)~=p
        warning('size(U,2)(=%g) ~= p(=%g); filling up U with zeros',size(U,2),p);
        U = [U , zeros(size(U,1),p-size(U,2))];
        size(U)
    end
    ubr = U(end,:);			% bottom row of U
    Z   = (eye(p)+ubr'*ubr/(1-ubr*ubr'))*...
        U(1:end-1,:)'*U(2:end,:);
    % least squares fit (left side is pseudo-inverse of right side)
    b  = eig(Z);			    % eigenvalues yield the decays/frequencies
    
    lb_dt = -log(abs(b))/pi;	  % (normalised) line widths
    f_dt  = angle(b)/(2*pi);	  % (normalised) frequencies
    
    % determine (complex) concentrations
    bm = repmat(b.',[n1 1]).^repmat((0:n1-1).',[1 p]);
    % a = pinv(bm)*fid;           % LSQ for concentrations
    a = bm\fid;                   % LSQ (slightly less expensive)
    
    % least square fit of amplitudes and phases
    % a = exp(log(b)*(0:n1-1)).'\fid;
    
    % create index with sensible values and within bound
    ind = logical((lb_dt>0)&(f_dt>bound(1))&(f_dt<bound(2)));
    if sum(ind)==0, warning('index empty'); end
    if verb>1, fprintf('Excluding: %g/%g\n',sum(~ind),p); end
    
    % generate fitted FID
    if (nargout>3) || plt
        % fid_fit_tmp = (a(ind)*ones(1,n1)).*bm(ind,:);
        if sum(ind)>0
            fid_fit(:,ll) = sum((repmat(a(ind).',[n1 1])).*bm(:,ind),2);
        end
    end
    
    % write output into array
    a_arr(:,ll) = a.';
    lb(:,ll)    = lb_dt.'*bw;
    freq(:,ll)  = f_dt.'*bw;
    ind_arr(:,ll) = ind.';
    
    
    %% plotting
    if plt
        % data processing
        spec_fit = ifftshift(fft(fid_fit(:,ll),[],1),1);
        ff = freq(:,ll); ii = ind_arr(:,ll); ff = ff(ii);
        mas = max(abs(spec(:,ll)));
        maf = max(abs(fid(:)));
        
        % actual plotting
        % plot FID
        subplot(2,1,1);
        plot(t,abs(fid),'b',t,abs(fid_fit(:,ll)),'g',...
            t,abs(fid-fid_fit(:,ll)),'r','LineWidth',plw(1));
        hold on;
        if n<n1
            plot([t(n) t(n)],[-0.1*maf 1.1*maf],'c:','LineWidth',plw(2));
        end
        plot([t(1) t(end)],[0 0],'k');
        hold off;
        axis([0 t(end) -0.01*maf maf]);
        xlabel('time'); ylabel('FID');
        title(sprintf('HSVD: b=original; g=fit; r=residual (%g/%g)',ll,n2));
        
        % plot spectrum
        subplot(2,1,2);
        if ~isempty(ff)
            plot(repmat(ff.',[2 1]),[-2;2]*ones(1,length(ff))*mas,'c:',...
                'LineWidth',plw(2));
        end
        hold on
        plot(x_ax,abs(spec(:,ll)),'b',x_ax,abs(spec_fit),'g',...
            x_ax,abs(spec(:,ll)-spec_fit),'r','LineWidth',plw(1));
        plot([x_ax(1) x_ax(end)],[0 0],'k');
        hold off
        
        xlabel('frequency'); ylabel('spectrum');
        axis([min(x_ax) max(x_ax) -0.01*mas mas]);
        drawnow;
    end
end


%% suppressed FID
if (nargout>3), fid_sup = fid_arr - fid_fit; end


end      % main function hsvd.m
