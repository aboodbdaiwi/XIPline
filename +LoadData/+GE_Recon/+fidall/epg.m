function sig = epg(T1,T2,te,tr,fa,ph,spl,acq,nstates,nreps,verb)
%EPG Extended Phase Graph simulations
%     sig = epg(T1,T2,te,tr,fa,ph,spl,acq,nstates,nreps,verb)
%      T1   T1 relaxation                                      [s]
%      T2   T2 relaxation                                      [s]
%      te   Echo time                                          [s]
%      tr   Repetition time                                    [s]
%      fa   Flip angle                                         [rad]
%      ph   RF phase                                           [rad]
%     spl   Spoiler type: 0=bSSFP,1=grd-spoiled SSFP,2=complete
%     acq   Acquire time step
% nstates   #EPG states to simulate and average                        (20)
%   nreps   #repetitions of whole train                                 (2)
%    verb   Verbose mode                                             (true)
%
%     sig   Signal evolution (#exc,#sim)
%
% Literature: jmri41_266_2015_weigel.pdf
%  Initial version from Guido Buonincontri
% 12/2024 Rolf Schulte
if nargin<1, help(mfilename); return; end


%% input parameters
acq_ph = -pi/2;                   % constant acquisition phase offset
T1 = T1(:).'; T2 = T2(:).';
te = te(:); tr = tr(:); fa = fa(:); 
ph = ph(:); spl = spl(:); acq = acq(:);

if isempty(T1), error('isempty(T1)'); end
if isempty(T2), error('isempty(T2)'); end
if isempty(te), error('isempty(te)'); end
if isempty(tr), error('isempty(tr)'); end
if isempty(fa), error('isempty(fa)'); end
if isempty(ph), error('isempty(ph)'); end
if isempty(spl), error('isempty(spl)'); end
if isempty(acq), error('isempty(acq)'); end

if ~exist('nstates','var'),  nstates = []; end
if isempty(nstates),         nstates = 20; end
if ~exist('nreps','var'),    nreps = []; end
if isempty(nreps),           nreps = 2; end
if ~exist('verb','var'),     verb = []; end
if isempty(verb),            verb = true; end


%% expand parameters
nexc = max([size(te,1),size(tr,1),size(fa,1),size(ph,1),size(spl,1),size(acq,1)]);
te = sub_expand(te,nexc);
tr = sub_expand(tr,nexc);
fa = sub_expand(fa,nexc);
ph = sub_expand(ph,nexc);
spl = sub_expand(spl,nexc);
acq = logical(sub_expand(acq,nexc));
nsims = size(T1,2);


%% checks
if max(abs(tr))>0.1
    warning('max(abs(tr))(=%g)>0.1: unit second!?',max(abs(tr))); 
end
if max(abs(fa))>(pi+0.001)
    warning('max(abs(fa))(=%g)>pi: unit rad!?',max(abs(fa)));
end


%% record execution time
if verb>0
    timerVal = tic;                              % record execution time
    fprintf('Starting EPG simulations (nreps=%g)\n',nreps);
    sub_check_input_pars(T1,T2,te,tr,fa,ph,spl,acq);
end


%% initialise values
Fp = complex(zeros(nstates,nsims,'single'));
Fm = complex(zeros(nstates,nsims,'single'));
Zk = zeros(nstates,nsims,'single');
Zk(1,:) = 1;                      % equilibrium magnetisation
sig = complex(zeros(nexc,nsims,'single'));


%% actual EPG simulation
for lrep = 1:nreps
    TE1old = inf; TZ1old = inf;
    for ll = 1:nexc
        % get current timings
        TE1 = te(ll);             % echo time
        TZ1 = tr(ll)-TE1;         % TR-TE -> remaining time

        % calc relaxation matrices
        if (ll==1) || abs(TE1-TE1old)>1d-10
            E1a = exp(-TE1./T1);
            E2a = exp(-TE1./T2);
        end
        if (ll==1) || abs(TZ1-TZ1old)>1d-10
            E1b = exp(-TZ1./T1);
            E2b = exp(-TZ1./T2);
        end
        
        % apply RF pulse
        [Fp,Fm,Zk] = sub_rot(Fp,Fm,Zk,fa(ll),ph(ll));
        
        % relaxation and dephasing until TE
        [Fp,Fm,Zk] = sub_free_precession(Fp,Fm,Zk,E1a,E2a);
        
        % acquire signal
        sig(ll,:) = Fp(1,:);
        
        % relaxation and dephasing until next TR
        [Fp,Fm,Zk] = sub_free_precession(Fp,Fm,Zk,E1b,E2b);
        
        % Apply epg shift due to gradient
        [Fp,Fm] = sub_epg_shift(Fp,Fm,spl(ll));
        
        TE1old = TE1;
        TZ1old = TZ1;
    end
end


%% acquisition phase correction + acq selection
sig = bsxfun(@times,sig,exp(-1i*(ph+acq_ph)));
if all((sum(abs(imag(sig)),1)./sum(abs(real(sig)),1))<1d-10)
    if verb, fprintf('Real data; discarding imaginary part\n'); end
    sig = real(sig);
end
sig = sig(acq,:);                 % discard not acquired data


%% finishing up
if verb>0, fprintf('epg: runtime = %.2f [s]\n',toc(timerVal)); end


end      % epg.m


%% sub-function
% RF pulse rotation
function [fPost1,fPost2,fPost3] = sub_rot(fPre1,fPre2,fPre3,fa,phi)
% Eq.[18]
ca = cos(fa);
sa = sin(fa);

% calculate T matrix
T11 = 0.5*(1+ca);
T12 = 0.5*(1-ca) * exp(1i*2*phi);
T13 = -1i*sa * exp(1i*phi);
T31 = -0.5*1i*sa * exp(-1i*phi);
T33 = ca;

% apply T matrix
fPost1 = T11.*fPre1 + T12.*fPre2 + T13.*fPre3;
fPost2 = conj(T12).*fPre1 + conj(T11).*fPre2 + conj(T13).*fPre3;
fPost3 = T31.*fPre1 + conj(T31).*fPre2 + T33.*fPre3;

end      % sub_rot


%% free precession
function [Fp,Fm,Zk] = sub_free_precession(Fp,Fm,Zk,E1,E2)
% Eq.[24]

% relaxation and dephasing
Zk = bsxfun(@times,E1,Zk);             % T1 relaxation for Z with k > 0
Zk(1,:)= Zk(1,:) + (1-E1);             % T1 recovery for Z with k = 0
Fp = bsxfun(@times,E2,Fp);             % T2 relaxation for F+ and F-
Fm = bsxfun(@times,E2,Fm);             % T2 relaxation for F+ and F-

end      % sub_free_precession


%% spoiling
function [Fp,Fm] = sub_epg_shift(Fp,Fm,type)
% Eq.[22] (type==1)

switch type
    case 0    % bSSFP
        % warning('requires off freq');
    case 1    % gradient-spoiled SSFP
        Fp(2:end,:)   = Fp(1:end-1,:);           % Dephase F+
        Fm(1:end-1,:) = Fm(2:end,:);             % Dephase F-
        Fm(end,:) = 0;
        Fp(1,:) = conj(Fm(1,:));

    case 2    % complete spoiling
        Fp = 0*Fp;
        Fm = 0*Fm;
        
    otherwise
        error('type(=%g) unkown',type);
end

end      % sub_epg_shift


%% if scalar -> expand parameter to list
function par = sub_expand(par,n)
switch size(par,1)
    case 1, par = repmat(par,[n 1]);
    case 2, par = [par(1,1);repmat(par(2,1),[n-1 1])];
    case n
    otherwise, warning('size(par,1)(=%d)~=1,2,n(=%d)',size(par,1),n);
end

end      % sub_expand


%% check input parameters
function sub_check_input_pars(t1,t2,te,tr,fa,ph,spl,acq)
if any(size(t1)~=size(t2)),  warning('size(t1)~=size(t2)'); end
if size(t1,1)~=1, warning('size(t1,1)~=1'); end
if any(size(fa)~=size(ph)),  warning('size(fa)~=size(ph)'); end
if any(size(fa)~=size(spl)), warning('size(fa)~=size(spl)'); end
if any(size(fa)~=size(te)),  warning('size(fa)~=size(te)'); end
if any(size(fa)~=size(tr)),  warning('size(fa)~=size(tr)'); end
if any(size(fa)~=size(acq)), warning('size(fa)~=size(acq)'); end
if any(tr<te), warning('tr<te'); end
if any(spl>2), warning('spl>2'); end
if any(spl<0), warning('spl<0'); end

if max(fa)>2*pi, warning('max(fa)(=%g[rad])>2*pi',max(fa)); end
if min(fa)<0, warning('min(fa)(=%g[rad])<0',min(fa)); end
if max(abs(ph))>2*pi, warning('max(abs(ph))(=%g[rad])>2*pi',max(abs(ph))); end

end      % sub_check_pars
