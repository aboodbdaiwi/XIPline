function NSA(nte,df0,freqs,f_str,dt,use_pinv)
% NSA Number of Signal Averages simulation
%  NSA(nte,df0,freqs,f_str,dt)
%     nte  number of echoes
%     df0  B0 induced offresonance          [Hz]                       (20)
%   freqs  Freqiencies                      [Hz]  ([0 -125 -215 -392 -716])
%   f_str  String for labelling frequencies on plot   
%                                 (['Lac','Pyr-Hyd','Ala','Pyr','Bi-Carb'})
%      dt  Time points                      [s]           (1e-4:1e-5:16e-3)
%use_pinv  Use pinv (instead of inv)                                (false)
%
% SR Reeder, JH Brittain, TM Grist, YF Yen
% Least-Square Chemical Shift Separation for 13C Metabolic Imaging
% Additionally dB0-inhomogeneity was included into the model according to:
% Anm = Sum_k { w_k exp[ i*2*pi*( freq_n-dfreq_k )*TE_m ] }
% with:
% n     ... index for nth Metabolite
% m     ... index for mth Echo time
% k     ... summation index for linewidth dispersion
% w_k   ... normalized linewidth function/weight [i.e. Sum_k(w_k)=1]
%
% 2/2009  Florian Wiesinger
% 2/2021  Rolf Schulte
if (nargin<1), help(mfilename); return; end


%% input parameters
if abs(nte-floor(nte))>1d-10, error('nte not even'); end
if ~exist('df0','var'), df0 = []; end
if isempty(df0),        df0 = 20; end
if ~exist('freqs','var'), freqs = []; end
if isempty(freqs)
    freqs = [0 -125 -215 -392 -716];
    f_str = {'Lac','Pyr-Hyd','Ala','Pyr','Bi-Carb'};
end
if ~exist('f_str','var'), f_str = []; end
if isempty(f_str)
    for l=1:length(freqs), f_str{l} = [num2str(freqs(l)) ' Hz']; end
end
if ~exist('dt','var'), dt = []; end
if isempty(dt), dt = (1e-4:1e-5:16e-3); end
if ~exist('use_pinv','var'), use_pinv = []; end
if isempty(use_pinv),        use_pinv = false; end


%% parameters
offf = linspace(-1000,1000,10000);
sigma = sqrt(-(df0/2)^2/(2*log(0.5)));
w = exp(-offf.^2/(2*sigma^2));
w = w/(sum(w));
nsam = zeros(length(freqs),length(dt));


%% actual NSA calculation
for i1=1:length(dt)
    tn = (0:dt(i1):(nte-1)*dt(i1));
    % Setting up encoding matrix A
    [fm,tm] = meshgrid(freqs,tn);
    offres = meshgrid(sum(meshgrid(w,tn).*exp(1i*tn.'*offf*2*pi),2),freqs).';
    A = exp(1i*2*pi*fm.*tm).*offres;
    if use_pinv
        nsam(:,i1) = real(1./diag(pinv(A'*A)));
    else
        nsam(:,i1) = real(1./diag(inv(A'*A)));
    end
end


%% plotting
clf; set(gca,'FontSize',14)
plot(dt*1e3,nsam.','LineWidth',1.5);
grid on
xlabel('TE [ms]'), ylabel(['NSA with ',num2str(nte),' echoes']),
axis([0,max(dt)*1e3,0,nte])
legend(f_str)


end      % NSA.m
