function hyper_flip(TR,n,T1,phi_test_deg)
%HYPER_FLIP  Calculate SNR-optimal flip angle for a single resonance
%  acquired with an FID sequence with non-recoverable magnetisation
%
% hyper_flip(TR,n,T1,phi_test)
%       TR   Repetition time                               [s]       (1)
%        n   Number of excitations                                   (64)
%       T1   T1 relaxation time                            [s]       (25)
% phi_test   Flip angle for testing (optional)             [deg]     ([])
%
% Mathematics:
% The signal for the lth excitation (l=[1:n]) is:
% signal_n = sin(phi)*x^(n-1), with x = cos(phi)*exp(-TR/T1)
% sum(Signal_n) = sin(phi)*(1-x^n)/(1-x)
%
% The SNR is calculated according to:
% SNR = sum(signal_n) / sqrt(n)
% hence relative to the SNR for a single 90deg excitation
%
% 12.3.2009  Florian Wiesinger
% 11/2020  Rolf Schulte
if (nargin<1), help(mfilename); return; end


%% default input parameters
if ~exist('TR','var'), TR = []; end
if isempty(TR), TR = 1; end
if ~exist('n','var'),  n = []; end
if isempty(n), n = 64; end
if ~exist('T1','var'), T1 = []; end
if isempty(T1), T1 = 25; end
if ~exist('phi_test_deg','var'), phi_test_deg = []; end


%% calculations
phi = linspace(0,90,1000)*pi/180;      % flip angles [rad]
x = cos(phi)*exp(-TR/T1);
SNR = sin(phi).*(1-x.^n)./(1-x)/sqrt(n);

SNR_opt = max(SNR);
phi_opt = phi(SNR==SNR_opt);
signal_opt = sin(phi_opt)*(cos(phi_opt)*exp(-TR/T1)).^(0:n-1);


%% print + plot info
clf, set(gcf,'DefaultAxesFontSize',14);
subplot(2,1,1),
plot(phi*180/pi,SNR,'LineWidth',2),
line(phi_opt*180/pi*[1 1],[0 1],'color','r');
xlabel('\phi  [deg]'); ylabel('SNR'); grid on,
title(sprintf('TR=%g s, n=%g, T_1=%g s',TR,n,T1));
axis([0,max(phi)*180/pi,0,ceil(max(SNR)*10)/10]);

subplot(2,1,2),
plot(0:(n-1),signal_opt,'LineWidth',2);
xlabel('excitation'); ylabel('signal'); grid on,
ma = max(signal_opt);
strg = sprintf('\\phi_{opt}=%.3g^{\\circ}, SNR=%.3g %%',...
    phi_opt*180/pi,SNR_opt*100);

fprintf('Optimal Flip Angle = %.3g [deg]; relative SNR = %.3g %%\n',...
    phi_opt*180/pi,SNR_opt*100);


%% info about test flip angle 
if ~isempty(phi_test_deg)
    phi_test = phi_test_deg*pi/180;
    signal_test = sin(phi_test)*(cos(phi_test)*exp(-TR/T1)).^(0:n-1);
    x = cos(phi_test)*exp(-TR/T1);
    SNR_test = sin(phi_test).*(1-x.^n)./(1-x)/sqrt(n);
    hold on, plot(0:(n-1),signal_test,'--m','linewidth',2);
    fprintf('Test Flip Angle = %g [deg]; relative SNR = %.3g %%\n',...
        phi_test_deg,SNR_test*1e2);
    legend('\phi_{opt}','\phi_{test}');
    strg = [strg sprintf('; \\phi_{test}=%.3g^{\\circ}, SNR=%.3g %%',...
        phi_test_deg,SNR_test*100)];
    ma = max([ma signal_test]);
    hold off
end
ma = ceil(ma*10)/10;
if n>1, axis([0,n-1,0,ma]); end
title(strg);


end      % hyper_flip.m
