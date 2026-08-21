function [nsam,dt,cnd] = NSA(nte,R2s,freqs,f_str,dt,use_pinv,t0)
% NSA Number of Signal Averages simulation
%  nsam = NSA(nte,R2s,freqs,f_str,dt,use_pinv,t0)
%     nte  number of echoes
%     R2s  Linewidth R2*                    [Hz]                       (20)
%   freqs  Freqiencies                      [Hz]  ([0 -125 -215 -392 -716])
%   f_str  String for labelling frequencies on plot   
%                                 (['Lac','Pyr-Hyd','Ala','Pyr','Bi-Carb'})
%      dt  Time points                      [s]           (1d-4:1d-5:16d-3)
%use_pinv  Use pinv (instead of inv)                                (false)
%      t0  Initial delay (TE)               [s]                    (0.5d-3)
%    nsam  #signal averages (#freqs,#dt)
%
% Literature: jmri26_1145_2007_reeder.pdf (adding Gaussian T2* decay)
%
% 1/2022  Rolf Schulte
if (nargin<1), help(mfilename); return; end


%% input parameters
if abs(nte-floor(nte))>1d-10, error('nte not even'); end
if ~exist('R2s','var'),   R2s = []; end
if isempty(R2s),          R2s = 20; end
if ~exist('freqs','var'), freqs = []; end
if isempty(freqs)
    freqs = [0 -125 -215 -392 -716];
    f_str = {'Lac','Pyr-Hyd','Ala','Pyr','Bi-Carb'};
end
if ~exist('f_str','var'), f_str = []; end
if isempty(f_str)
    for l=1:length(freqs), f_str{l} = [num2str(freqs(l)) ' Hz']; end
end
if ~exist('dt','var'),    dt = []; end
if isempty(dt),           dt = (1d-4:1d-5:16d-3); end
if ~exist('use_pinv','var'), use_pinv = []; end
if isempty(use_pinv),     use_pinv = false; end
if ~exist('t0','var'),    t0 = []; end
if isempty(t0),           t0 = 0.5d-3; end

freqs = freqs(:).';
if length(R2s)==1, R2s = R2s*ones(1,length(freqs)); end
R2s = R2s(:).';


%% actual NSA calculation
nsam = zeros(length(freqs),length(dt));
cnd = zeros(1,length(dt));
for ll=1:length(dt)
    tn = (0:dt(ll):(nte-1)*dt(ll)).'+t0;

    % Setting up encoding matrix A
    offres = exp(-(pi*tn*R2s).^2/(4*log(2)));
    A = exp(1i*2*pi*tn*freqs).*offres;
    
    if use_pinv
        nsam(:,ll) = real(1./diag(pinv(A'*A)));
    else
        nsam(:,ll) = real(1./diag(inv(A'*A)));
    end
    cnd(1,ll) = cond(A'*A);
end


%% plotting
clf; set(gca,'FontSize',14)
plot(dt*1e3,nsam.','LineWidth',1.5);
grid on
xlabel('TE [ms]'), ylabel(['NSA with ',num2str(nte),' echoes']),
axis([0,max(dt)*1e3,0,nte])
legend(f_str)


end      % NSA.m
