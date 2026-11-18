function [disfitObj, qc] = fit_xe_dissolved_robust(disData1_avg, t, useVoigt)
%FIT_XE_DISSOLVED_ROBUST Robust three-component hyperpolarized 129Xe FID fit.
%
% Fits the original complex FID (real and imaginary simultaneously) using
% variable projection:
%   1. Scale only the optimization target to order unity.
%   2. Nonlinearly optimize resonant frequencies and linewidths.
%   3. At each iteration, solve the three COMPLEX component amplitudes
%      exactly by linear least squares. No zero-phase assumption is needed.
%   4. Restore the original amplitude units and populate the existing
%      Calibration.NMR_TimeFit or NMR_TimeFit_v object.
%
% Requires MATLAB Optimization Toolbox (lsqnonlin) and XIPline Calibration.
%
% Usage:
%   [disfitObj, qc] = fit_xe_dissolved_robust(disData1_avg, t, false);
%   [disfitObj, qc] = fit_xe_dissolved_robust(disData1_avg, t, true);
%
% Component order: RBC (near 0 Hz), Barrier (~-720 Hz), Gas (~-7700 Hz).
% Initial guesses and frequency bounds below are intended for the
% particular 3T acquisition shown; adjust if center frequency changes.
%
% IMPORTANT: Does NOT remove early FID points or apply line broadening.
% Inspect the complex residual before accepting physiological parameters.

if nargin < 3
    useVoigt = false;
end
useVoigt = logical(useVoigt);

fid = double(disData1_avg(:));
t = double(t(:));

assert(numel(fid) == numel(t) && numel(fid) >= 20, ...
    'FID and time vectors must have the same length (at least 20 points).');
assert(all(isfinite(real(fid))) && all(isfinite(imag(fid))) && ...
       all(isfinite(t)) && all(diff(t) > 0), ...
    'FID must be finite, and time must be finite and strictly increasing.');

% Time must be in seconds, as required by Calibration.NMR_Mix.
dt = median(diff(t));
if dt > 1e-3 || dt <= 0
    warning('Sample interval is %.4g s. Verify time is in SECONDS.', dt);
end

ampScale = max(abs(fid));
assert(isfinite(ampScale) && ampScale > 0, 'FID has no valid signal.');
y = fid / ampScale;

% Parameters: [f_RBC f_Barrier f_Gas  LW_RBC LW_Barrier LW_Gas (GW_Barrier)]
p0 = [0, -720, -7700,   240, 240, 40];
lb = [-450, -1500, -8800,  10, 10, 5];
ub = [ 450,  -150, -6500, 3000, 3000, 750];

if useVoigt
    % XIPline's NMR_Mix_v uses a Gaussian factor only in component 2
    % when RBC and gas Gaussian widths are fixed at zero here.
    p0 = [p0, 150];
    lb = [lb, 0];
    ub = [ub, 1200];
end

options = optimoptions('lsqnonlin', ...
    'Algorithm', 'trust-region-reflective', ...
    'Display', 'off', ...
    'MaxIterations', 1000, ...
    'MaxFunctionEvaluations', 30000, ...
    'FunctionTolerance', 1e-11, ...
    'StepTolerance', 1e-11, ...
    'OptimalityTolerance', 1e-10);

[p,~,residual,exitflag,output] = lsqnonlin( ...
    @(par) complex_residual(par,t,y,useVoigt), p0, lb, ub, options);

% Recover optimal complex amplitudes (relative to the scaled input).
X = spectrum_basis(p,t,useVoigt);
a = X \ y;
fittedFID = ampScale * (X*a);

% Reconstruct the standard Calibration fit object in original signal units.
area  = ampScale * abs(a).';
freq  = p(1:3);
fwhm  = p(4:6);
phase = rad2deg(angle(a)).';
fwhmG = [0,0,0];
if useVoigt
    fwhmG(2) = p(7);
    disfitObj = Calibration.NMR_TimeFit_v(fid, t, area, freq, fwhm, ...
        fwhmG, phase, 0, numel(t));
else
    disfitObj = Calibration.NMR_TimeFit(fid, t, area, freq, ...
        fwhm, phase, 0, numel(t));
end

qc = struct();
qc.relativeComplexFIDError = norm(fid - fittedFID)/norm(fid);
qc.relativeMagnitudeFIDError = norm(abs(fid)-abs(fittedFID))/norm(abs(fid));
qc.exitflag = exitflag;
qc.output = output;
qc.frequency_Hz = freq;
qc.lorentzianFWHM_Hz = fwhm;
qc.gaussianFWHM_Hz = fwhmG;
qc.area = area;
qc.phase_deg = phase;
qc.originalSignal = fid;
qc.fittedSignal = fittedFID;

fprintf('\n--- 129Xe three-component complex FID fit ---\n');
fprintf('Model: %s\n', ternary(useVoigt,'Lorentzian + Voigt barrier','Lorentzian'));
fprintf('Optimizer exitflag: %d\n', exitflag);
fprintf('Relative complex FID error: %.4f\n',qc.relativeComplexFIDError);
fprintf('Relative magnitude FID error: %.4f\n',qc.relativeMagnitudeFIDError);
fprintf('        Component      Area        Freq(Hz)   FWHM(Hz)    Phase(deg)\n');
names = {'RBC','Barrier','Gas'};
for k = 1:3
    fprintf('%12s   %10.3g  %+10.2f  %10.2f  %+10.2f\n', ...
        names{k}, area(k), freq(k), fwhm(k), phase(k));
end
if useVoigt
    fprintf('Barrier Gaussian width parameter: %.2f Hz\n',p(7));
end
if exitflag <= 0
    warning('Optimizer did not converge; examine the fitted parameters and residual.');
end

% Plot measured versus fitted using the SAME time and FFT grid.
N = numel(t);
freqAxis = (-floor(N/2):ceil(N/2)-1).'/(N*dt);
measSpec = dt*fftshift(fft(fid));
fitSpec  = dt*fftshift(fft(fittedFID));

fig = figure('Color','w','Name','129Xe Complex FID Fit QC', ...
    'Position',[100 100 1200 850]);
tiledlayout(3,1,'TileSpacing','compact','Padding','compact');

nexttile;
plot(freqAxis,abs(measSpec),'b.','DisplayName','Measured'); hold on;
plot(freqAxis,abs(fitSpec),'r-','LineWidth',1.5,'DisplayName','Fitted');
xlim([-10000 3000]); grid on; legend('Location','best');
ylabel('Magnitude');
title(sprintf('Dissolved-phase spectrum | Relative complex FID error = %.3f', ...
    qc.relativeComplexFIDError));

nexttile;
plot(t*1e3,abs(fid),'b.','DisplayName','Measured'); hold on;
plot(t*1e3,abs(fittedFID),'r-','LineWidth',1.5,'DisplayName','Fitted');
grid on; legend('Location','best');
xlabel('Time (ms)'); ylabel('|FID|');
title('Time-domain magnitude');

nexttile;
plot(t*1e3,real(fid-fittedFID),'k-', 'DisplayName','Real residual'); hold on;
plot(t*1e3,imag(fid-fittedFID),'m-', 'DisplayName','Imaginary residual');
grid on; legend('Location','best');
xlabel('Time (ms)'); ylabel('Residual (a.u.)');
title('Complex time-domain residual');

% Save Fitting Figure
outputFolder = fullfile(pwd, 'GasExchange_Analysis');

if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

% Save high-resolution PNG
exportgraphics(fig, ...
    fullfile(outputFolder, 'Xe129_Dissolved_Fit.png'), ...
    'Resolution', 300);

% Save editable MATLAB figure
savefig(fig, ...
    fullfile(outputFolder, 'Xe129_Dissolved_Fit.fig'));

fprintf('Fit figures saved to: %s\n', outputFolder);

end

function r = complex_residual(p,t,y,useVoigt)
X = spectrum_basis(p,t,useVoigt);
a = X\y;
err = X*a - y;
r = [real(err); imag(err)];
end

function X = spectrum_basis(p,t,useVoigt)
f = p(1:3);
lw = p(4:6);
X = zeros(numel(t),3);
for j = 1:3
    X(:,j) = exp(-pi*lw(j)*t + 1i*2*pi*f(j)*t);
end
if useVoigt
    % Match the exact Gaussian damping in XIPline's NMR_Mix_v:
    X(:,2) = X(:,2).*exp(-4*log(2)*(p(7)^2)*(t.^2));
end
end

function s = ternary(c,yes,no)
if c, s=yes; else, s=no; end
end
