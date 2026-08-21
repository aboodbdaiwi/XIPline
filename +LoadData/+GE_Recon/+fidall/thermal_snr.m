function thermal_snr(TR,phi)
%THERMAL_SNR  Calculate SNR per unit time for different TRs and flip angles
%   for (spoiled) pulse-and-acquire and thermal polarisation
% thermal_snr(TR,phi)
%
%    TR  Relative repetition times: normalised->TR/T1         (0.01:0.01:5)
%   phi  Excitation flip angle                        [deg]          (1:90)
%
% Literature: Handbook of MRI Pulse Seq. Bernstein; p587, Eq14.8
%
%  9/2022  Rolf Schulte

% if (nargin<1), help(mfilename); return; end
plt_sat = false;


%% default input parameters
t_unit = 1;                       % unit time [s]
T1 = 1;                           % normalise TR to T1
if ~exist('TR','var'),  TR = []; end
if isempty(TR),         TR = (0.01:0.01:5); end
if ~exist('phi','var'), phi = []; end
if isempty(phi),        phi = (1:90); end


%% signal, snr and saturation
phi = phi(:);                     % excitation flip angle [deg]
TR = TR(:).';                     % repetition time [s]
phi_rad = phi*pi/180;             % flip angle [rad]
S_initial = sin(phi_rad);         % signal of first acquisition (full M0)
S_spoil = (sin(phi_rad)*(1-exp(-TR/T1)))./...
    (1-cos(phi_rad)*exp(-TR/T1));                % steady-state signal
snr = bsxfun(@times,S_spoil,sqrt(t_unit./TR));   % snr
sat = bsxfun(@rdivide,S_initial,S_spoil);        % saturation 


%% for Ernst angle
ernst_rad = acos(exp(-TR/T1));                   % Ernst angle [rad]
ernst = ernst_rad*180/pi;                        % Ernst angle [deg]

S_initial_ernst = sin(ernst_rad);           % 1st signal with Ernst angle
S_spoil_ernst = (sin(ernst_rad).*(1-exp(-TR/T1)))./...
    (1-cos(ernst_rad).*exp(-TR/T1));        % steady-state signal
snr_ernst = S_spoil_ernst.*sqrt(T1./TR);    % snr
snr_ernst = snr_ernst/max(snr_ernst);       % normalise to max
sat_ernst = S_initial_ernst./S_spoil_ernst; % saturation
sar_ernst = ernst_angle(TR,T1).^2./TR*T1;   % power deposition
sar_ernst = sar_ernst/max(sar_ernst);       % normalise to max



%% plotting
% SNR and saturation over all TR and phi
figure(10); clf;
set(gcf,'DefaultAxesFontSize',12);
if plt_sat, subplot(2,1,1); end
imagesc(TR,phi,snr);
title('SNR per unit time');
xlabel('TR/T1');
ylabel('phi [deg]');
hold on; plot(TR,ernst,'w','LineWidth',2); hold off
colorbar;

if plt_sat
    subplot(2,1,2);
    imagesc(TR,phi,sat,[1 3]);
    line([1 1],[200 200],'color','w')
    title('Saturation = S_{initial}/S_{spoil}');
    xlabel('TR/T1');
    ylabel('phi [deg]');
    hold on; plot(TR,ernst,'w'); hold off
    colorbar
end
drawnow


%% Ernst angle with SNR, saturation
figure(11); clf
set(gcf,'DefaultAxesFontSize',12);
subplot(2+plt_sat,1,1);
plot(TR,ernst,'b','LineWidth',1.5);
title('Ernst angle');
xlabel('TR/T1');
ylabel('phi [deg]');
grid on;

subplot(2+plt_sat,1,2);
plot(TR,snr_ernst*100,'b',TR,sar_ernst*100,'r','LineWidth',1.5);
title('SNR and SAR per unit time at Ernst angle');
xlabel('TR/T1');
ylabel('[%]');
legend({'SNR','SAR'});
grid on;

if plt_sat
    subplot(3,1,3);
    plot(TR,sat_ernst);
    title('Saturation');
    xlabel('TR/T1');
    ylabel('S_{initial}/S_{spoil}');
    grid on;
end


%% print info
tr_t1 = [0.2 0.5 1 1.5 2 5];
snr0 = snr_ernst(1);
sar0 = ernst_angle(0.01,1).^2./0.01;
% plot(t,ernst_angle(t,1).^2./t);grid on

for l=1:length(tr_t1)
    [~,ii] = min(abs(TR/T1-tr_t1(l)));
    fprintf('TR/T1 = %3g: SNR = %g [%%]; SAR = %g [%%]\n',...
        tr_t1(l), round(snr_ernst(ii)/snr0*100),...
        round(ernst_angle(tr_t1(l),1).^2/tr_t1(l)/sar0*100));
end


end      % thermal_snr.m
