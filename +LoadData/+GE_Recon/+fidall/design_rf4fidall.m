function [rf,par]=design_rf4fidall(flip,t,bw,ftw,sv,plt,ref,reg,wght,...
    rftype,oversamp,b1max_ref,dt)
%DESIGN_RF4FIDALL  Design selective RF excitation pulse
%  Small-tip angle approximation -> linear least-squares
% [rf,par] = design_rf4fidall(flip,t,bw,ftw,sv,plt,ref,reg,wght,rftype,...
%                             oversamp,b1max_ref,dt)
%
%     flip   Flip angle                                    [deg] (90)
%        t   Duration of RF pulse                          [s]   (5d-3)
%       bw   Bandwidth                                     [Hz]  (2.1d3)
%            slice bw=gyrogamma(nucleus)/2/pi*gmax*thick
%       alternatively: f_pass=passband frequencies (#pb,2) [Hz]
%      ftw   Fractional-transition width                         ([])
%       alternatively: f_stop=stopband frequencies (#sb,2) [Hz]
%       sv   Save RF waveform (Create5xxRfExtFile)               (false)
%            logical or fname
%      plt   Plotting (2=with profiles)                          (2)
%      ref   Time reference [-1..1]                              (0.6)
%            0=midle=regular linear phase pule
%            1=self-refocused
%      reg   Tikhonov regularisation for linear LSQ fit          (1d-2)
%     wght   Weights on stopband to reduce error                 (1)
%   rftype   'spat'=spatially; 'freq'=frequency-selective        ('spat')
% oversamp   Oversampling of spectral domain for fit             (10)
%b1max_ref   B1max reference of fidall                     [G]   (0.1281)
%       dt   Sampling dwell time                           [s]   (4d-6)
%
%       rf   RF waveform
%      par   Parameter structure
%
% See also Create5xxRfExtFile
% 11/2022 Rolf Schulte
if nargin<1, help(mfilename); return; end

if ~exist('flip','var'),     flip = []; end
if isempty(flip),            flip = 90; end
if ~exist('t','var'),        t = []; end
if isempty(t),               t = 5d-3; end
if ~exist('bw','var'),       bw = []; end
if isempty(bw),              bw = 2.1d3; end
if ~exist('ftw','var'),      ftw = []; end
% if isempty(ftw),             ftw = 0.15; end
if ~exist('ref','var'),      ref = []; end
if isempty(ref),             ref = 0.6; end
if ~exist('reg','var'),      reg = []; end
if isempty(reg),             reg = 1d-2; end
if ~exist('wght','var'),     wght = []; end
if isempty(wght),            wght = 1; end
if length(wght)~=1, warning('length(wght)(=%g)~=1',length(wght)); end
if ~exist('sv','var'),       sv = []; end
if isempty(sv),              sv = false; end
if ~exist('plt','var'),      plt = []; end
if isempty(plt),             plt = 2; end
if ~exist('rftype','var'),   rftype = ''; end
if isempty(rftype),          rftype = 'spat'; end
if ~exist('oversamp','var'), oversamp = []; end
if isempty(oversamp),        oversamp = 10; end
if ~exist('b1max_ref','var'),b1max_ref = []; end
if isempty(b1max_ref),       b1max_ref = 0.1281; end
if ~exist('dt','var'),       dt = []; end
if isempty(dt),              dt = 4d-6; end


%% extract pass- and stopband frequencies
if (size(bw,2)~=size(ftw,2)) && ~isempty(ftw)
    error(' size(bw,2)(=%g)~=size(ftw,2)(=%g)', size(bw,2),size(ftw,2));
end
if size(bw,2)==1
    if size(bw,1)~=1,  error('size(bw,1)(=%g)~=1', size(bw,1)); end
    if size(ftw,1)~=1, error('size(ftw,1)(=%g)~=1',size(ftw,1)); end
    f_pass = bw*(1-ftw)*[-1 1]/2;
    f_stop(1,:) = [-1/dt/2 -bw*(1+ftw)/2];
    f_stop(2,:) = [bw*(1+ftw)/2 1/dt/2];
    plt_bw = 2*ceil(2*bw*(1+ftw)/10)*10;
    bwcs = bw;
else
    f_pass = bw;
    f_stop = ftw;
    plt_bw = 2*ceil(2*max(abs([f_stop(:) ; f_pass(:)]))/10)*10;
    bwcs = diff(f_pass(1,:));
end
if size(f_pass,2)~=2, error('size(f_pass,2)(=%g)~=2',size(f_pass,2)); end
if ~isempty(f_stop) && size(f_stop,2)~=2
    error('size(f_stop,2)(=%g)~=2',size(f_stop,2)); 
end
if any(diff(f_pass,2)<0), error('f_pass not increasing'); end
if any(diff(f_stop,2)<0), error('f_stop not increasing'); end
if f_pass(1,1)>0, warning('f_pass(1,1)(=%g)>0',f_pass(1,1)); end
if f_pass(1,2)<0, warning('f_pass(1,2)(=%g)<0',f_pass(1,2)); end


%% misc pars
n = ceil(t/dt/2-1d-10)*2;         % #rf waveform points
if abs(t-dt*n)>1d-10
    fprintf('Attention: increasing t(=%g[us]) to %g[us]\n',t*1d6,dt*n*1d6);
    t = dt*n;
end


%% print pulse info
% flip,t,bw,ftw,sv,plt,ref,reg,wght,rftype,oversamp,b1max_ref,dt
fprintf('flip=%g[deg]; t=%g[ms]; ref=%g; reg=%g; ',flip(1),t*1d3,ref,reg);
fprintf('wght=%g; rftype=''%s''; oversamp=%g\n',wght,rftype,oversamp);


%% discretisation
isodelay = 0.5-ref/2;             % time from iso-ref to end of rf1
kf = ((-n/2:n/2-1)-ref*n/2)*dt;   % spectral k-space (aka time) [s]

nf = ceil(n*oversamp);            % #spectral points
ff = (-nf/2:nf/2-1).'/nf/dt;      % freq grid


%% setup target excitation profile b and mask
b = zeros(nf,1);                  % target profile
mask = false(nf,1);               % points inside stop and passband
W = ones(nf,1);                   % weighting
if length(flip)==1
    target = ones(size(f_pass,1),1);
else
    target = sind(flip(:))/sind(flip(1));
    if any(target>1), warning('target>1'); end
end
for l1=1:size(f_pass,1)
    ii = (ff>=f_pass(l1,1))&(ff<=f_pass(l1,2));
    b(ii,1) = target(l1);
    mask(ii,1) = true;
end
for l1=1:size(f_stop,1)
    ii = (ff>=f_stop(l1,1))&(ff<=f_stop(l1,2));
    mask(ii,1) = true;
    W(ii,1) = wght;
end
b = b(mask,1);
W = W(mask,1);
% fprintf('sum(mask) = %g\n',sum(mask));


%% setup encoding matrix A
A = exp(-2*pi*1i*ff(mask,1)*kf);
fprintf('encoding matrix: size(A)=(%g,%g); cond=',size(A));


%% actual linear LSQ fit
% pinv numerically less efficient but more stable
if abs(wght-1)>1d-10
    % weighted linear LSQ
    AW = bsxfun(@times,A',W.');
    if reg>0
        % Tikhonov regularisation -> reduce peak B1 + improve conditioning
        % rf = pinv(AW*A+reg*eye(size(A,2)))*(AW*b);
        [pinvAW,condA] = sub_pinv_cond(AW*A+reg*eye(size(A,2)));
    else
        % rf = pinv(AW*A)*(AW*b);
        [pinvAW,condA] = sub_pinv_cond(AW*A);
    end
    rf = pinvAW*(AW*b);
else
    % linear LSQ
    if reg>0
        % Tikhonov regularisation
        % rf = pinv(A'*A+reg*eye(size(A,2)))*(A'*b);
        [pinvA,condA] = sub_pinv_cond(A'*A+reg*eye(size(A,2)));
        rf = pinvA*(A'*b);
    else
        % rf = pinv(A)*b;
        [pinvA,condA] = sub_pinv_cond(A);
        rf = pinvA*b;
    end
end
fprintf('%g\n',condA);


%% determine if real-valued rf coefficient
sarrf = sum(abs(real(rf)));
sairf = sum(abs(imag(rf)));
if sarrf>100*sairf, fprintf('Real pulse\n'); rf = real(rf); end


%% print info
b1max = pi*flip(1)*n*max(abs(rf)) * 1d4 / (180 * t *sum(real(rf)) * gyrogamma(1));
fprintf('Isodelay=%g: TE+=%g[ms]\n',isodelay,isodelay*n*dt*1d3);
fprintf('1H B1max=%g[G]=%g[uT]; f1max=%g[Hz]\n',b1max,...
    b1max*100,b1max/1e4*gyrogamma(1)/2/pi);
if b1max>b1max_ref
    fprintf('Warning: b1max(=%g)>b1max_ref(=%g)\n',b1max,b1max_ref);
end


%% create fname
if ischar(sv)
    fname = sv;
    sv = true;
end
if sv
    if ~exist('fname','var')
        % flip,t,bw,ftw,sv,plt,ref,reg,wght,...
        % rftype,oversamp,b1max_ref,dt
        
        fname = sprintf('exrf_%s_flip%g_t%g_ref%g_b1max%.5g_bw%g',...
            rftype,round(flip(1)),t*1d3,ref,b1max*100,round(bwcs));
        if ~isempty(f_stop)
            fname = sprintf('%s_sb%g',fname,round(f_stop(end,1)));
        end
        fname = regexprep(fname,'\.','p');  % replace '.' by 'p'
        fname = regexprep(fname,'-','m');   % replace '-' by 'm'
    end
end


%% plotting
if plt>0
    figure(10); clf;
    tt = (1:n)*dt*1d3;
    if isreal(rf)
        plot(tt,rf,'b');
        ymax = max(rf); ymin = min(rf);
    else
        plot(tt,abs(rf),'b',tt,real(rf),'g:',tt,imag(rf),'r:');
        ymax = max([real(rf) ; imag(rf) ; abs(rf)]);
        ymin = min([real(rf) ; imag(rf)]);
    end
    xlabel('time [ms]'); 
    ylabel('B1')
    grid on;
    axis([0 n*dt*1d3 1.1*ymin 1.1*ymax]);
    title('RF pulse coefficients');
    drawnow;
    if sv, print([fname '_pulse.png'],'-dpng','-r300','-painters'); end
end


%% plotting profile
if plt>1
    intpl = 20;
    ni = n*intpl;
    ax = (-ni/2:ni/2-1)/ni/dt;
    Mxy = fftshift(fft(rf(end:-1:1,1),ni,1));
    phafu = exp(1i*(-ni/2:ni/2-1)*isodelay/pi).';
    Mxy = Mxy.*phafu;             % rewind phase to isodelay reference
    Mxy = Mxy(end:-1:1,1);
    
    figure(12); clf
    plot(ff(mask,1),b,'mx',...
        ax,abs(Mxy),'b:',ax,real(Mxy),'g',ax,imag(Mxy),'r:');
    title('Profile (fft)');
    xlabel('freq [Hz]');
    legend('target','|Mxy|','Mx','My')
    axis([-plt_bw/2 plt_bw/2 -1.1 1.1]);
    grid on;
    if sv, print([fname '_profile.png'],'-dpng','-r300','-painters'); end
else
    Mxy = []; ax = [];
end


%% determining error in pass- and stopbands
if plt==2
    mapb = false(length(ax),1);
    masb = mapb;
    for l1=1:size(f_pass,1)
        ii = (ax>f_pass(l1,1))&(ax<f_pass(l1,2));
        mapb(ii,1) = true;
    end
    for l1=1:size(f_stop,1)
        ii = (ax>f_stop(l1,1))&(ax<f_stop(l1,2));
        masb(ii,1) = true;
    end
    err = [max(abs(Mxy(mapb))-1) ...
        max(abs(Mxy(masb)))];
    fprintf('Maximum error passband=%.4g%%',err(1)*100);
    if any(masb)
        fprintf('; stopband=%.4g%%\n',err(2)*100);
    else
        fprintf('\n');
    end
else
    err = [];
end


%% exporting waveform
if sv
    fprintf('Saving RF to ''%s''\n',fname);
    rf = [rf;0;0];   % waveform must end with 0
    n = n+2;
    
    if ~(strcmp(rftype,'freq') || strcmp(rftype,'spat'))
        error('rftype(=''%s'')~=''freq'' or ''spat''',rftype);
    end
    par.rftype = rftype;
    par.pw = n*dt*1d6;
    par.isodelay = isodelay;
    par.bwcs = bwcs;
    par.maxB1 = b1max;
    par.minwidth = 1;
    par.flip = flip(1);
    if isreal(rf)
        Create5xxRfExtFile([fname '.rho'], rf/max(abs(rf)), par);
    else
        Create5xxRfExtFile([fname '.rho'],   abs(rf)/max(abs(rf)), par);
        Create5xxRfExtFile([fname '.theta'], angle(rf)/pi, par);
    end
    ax = single(ax); Mxy = single(Mxy);
    save([fname '.mat'],'flip','t','bw','ftw','f_pass','f_stop','ref',...
        'reg','wght','rftype','oversamp','b1max_ref','dt',...
        'par','rf','isodelay','b1max','Mxy','ax','err','condA');
else
    par = [];
end


end   % main function design_freqselrf.m


function [X,condA] = sub_pinv_cond(A,tol)
%PINV   Pseudoinverse.
%   X = PINV(A) produces a matrix X of the same dimensions
%   as A' so that A*X*A = A, X*A*X = X and A*X and X*A
%   are Hermitian. The computation is based on SVD(A) and any
%   singular values less than a tolerance are treated as zero.
%
%   PINV(A,TOL) treats all singular values of A that are less than TOL as
%   zero. By default, TOL = max(size(A)) * eps(norm(A)).
%
%   Class support for input A: 
%      float: double, single
%
%   See also RANK.
 
%   Copyright 1984-2015 The MathWorks, Inc. 

[U,S,V] = svd(A,'econ');
s = diag(S);
if nargout>1
    condA = max(s)/min(s);
end
if nargin < 2 
    tol = max(size(A)) * eps(norm(s,inf));
end
r1 = sum(s > tol)+1;
V(:,r1:end) = [];
U(:,r1:end) = [];
s(r1:end) = [];
s = 1./s(:);
X = (V.*s.')*U';


end      % sub-function sub_pinv_norm

