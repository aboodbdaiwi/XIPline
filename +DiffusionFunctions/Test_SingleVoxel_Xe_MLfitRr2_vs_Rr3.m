%% SINGLE-VOXEL 129Xe MORPHOMETRY: MLfitRr2 vs MLfitRr3
% Put MLfitRr2.m and MLfitRr3.m on the MATLAB path (or in +DiffusionFunctions).
% Uses the same data, noise variance, initial guess, bounds and constraints.
% No external image files required. Requires Optimization Toolbox.
clearvars; close all; clc;
rng(42,'twister');

% ===== EDIT THESE TO MATCH YOUR ACQUISITION =====
bvalues = [0 2 4 6 8 10 12];  % s/cm^2: REPLACE with your actual b-values
D0 = 0.14;                    % cm^2/s
Delta_ms = 3.5;               % milliseconds: REPLACE with actual Delta
deltaT = Delta_ms / 1000;     % seconds
R_true_um = 330;              % outer radius, micrometers
r_true_um = 140;              % inner radius, micrometers
S0_true = 1.0;
SNR_at_b0 = 30;                % S0 / Gaussian-channel noise sigma
nNoiseRealizations = 50;      % set 1 for a single noisy realization
% ===============================================

% Original fitting settings (human 129Xe)
lb_Rr = [0.0250 0.0075];
ub_Rr = [0.0450 0.0350];
x0_Rr = [0.0300 0.0140];

% Resolve either package functions or plain files.
% Check both functions inside the DiffusionFunctions package
assert(exist('DiffusionFunctions.MLfitRr2','file') == 2, ...
    'Cannot find DiffusionFunctions.MLfitRr2');

assert(exist('DiffusionFunctions.MLfitRr3','file') == 2, ...
    'Cannot find DiffusionFunctions.MLfitRr3');

bvalues = double(bvalues(:).');
assert(numel(bvalues)>=3 && all(isfinite(bvalues)) && all(bvalues>=0), 'Invalid b-values.');
assert(r_true_um < R_true_um && SNR_at_b0>0 && nNoiseRealizations>=1, 'Invalid simulation parameters.');

trueR = R_true_um * 1e-4; true_r = r_true_um * 1e-4;
trueSignal = xeSignal(trueR,true_r,S0_true,bvalues,D0,deltaT);
sigma = S0_true / SNR_at_b0;
s2 = sigma^2;
assert(all(isfinite(trueSignal)), 'True parameters generate an invalid 129Xe signal.');

fprintf('\n129Xe SINGLE VOXEL TEST\n');
fprintf('Truth: R=%.1f um, r=%.1f um, h=%.1f um, S0=%.3f\n', R_true_um,r_true_um,R_true_um-r_true_um,S0_true);
fprintf('D0=%.4f cm^2/s, Delta=%.3f ms, SNR(b0)=%.1f\n',D0,Delta_ms,SNR_at_b0);
fprintf('b-values: %s\n',mat2str(bvalues));

% Two controlled scenarios: exact signal, then one noisy realization.
noisySignal = abs(trueSignal + sigma*randn(size(trueSignal)) + 1i*sigma*randn(size(trueSignal)));
cases = {trueSignal, noisySignal};
caseNames = {'NOISE-FREE','RICIAN NOISY'};
fitResults = nan(2,2,3);
fitObjectives = nan(2,2);

for k = 1:2
    y = cases{k};
    % Use original caller's S0 estimation, but ensure the same x0/bounds for both functions.
    s0 = y(1);  % original CM_MorphometryFit uses M(j,1)
    lb = [lb_Rr, 0.5*s0];
    ub = [ub_Rr, 1.5*s0];
    x0 = [x0_Rr, s0];
    fprintf('\n--- %s ---\n',caseNames{k});
    try
        [xo,fo] = fitOld(y,bvalues,s2,x0,lb,ub);
        fitResults(k,1,:) = xo;
        fitObjectives(k,1) = fo;
        printFit('MLfitRr2',xo,fo,R_true_um,r_true_um);
    catch ME
        fprintf('MLfitRr2 FAILED: %s\n',ME.message);
    end
    try
        [xn,fn] = fitNew(y,bvalues,s2,x0,lb,ub);
        fitResults(k,2,:) = xn;
        fitObjectives(k,2) = fn;
        printFit('MLfitRr3',xn,fn,R_true_um,r_true_um);
    catch ME
        fprintf('MLfitRr3 FAILED: %s\n',ME.message);
    end
end

% Repeat noisy measurements for the SAME single voxel to test systematic bias.
Rfits = nan(nNoiseRealizations,2); rfits = nan(nNoiseRealizations,2);
valid = false(nNoiseRealizations,2);
for rep = 1:nNoiseRealizations
    y = abs(trueSignal + sigma*randn(size(trueSignal)) + 1i*sigma*randn(size(trueSignal)));
    s0 = y(1);
    lb = [lb_Rr 0.5*s0]; ub = [ub_Rr 1.5*s0]; x0 = [x0_Rr s0];
    for method = 1:2
        try
            if method == 1
                [x,f] = fitOld(y,bvalues,s2,x0,lb,ub);
            else
                [x,f] = fitNew(y,bvalues,s2,x0,lb,ub);
            end
            if all(isfinite(x)) && isfinite(f) && x(1)>x(2) && x(2)>0
                Rfits(rep,method) = x(1)*1e4;
                rfits(rep,method) = x(2)*1e4;
                valid(rep,method) = true;
            end
        catch
            % Failed fit counted below.
        end
    end
end
fprintf('\n--- %d RICIAN NOISE REALIZATIONS OF ONE VOXEL ---\n',nNoiseRealizations);
for method = 1:2
    label = 'MLfitRr2'; if method==2, label='MLfitRr3'; end
    good = valid(:,method);
    if any(good)
        fprintf('%s: success %d/%d; R = %.1f +/- %.1f um (bias %+0.1f); r = %.1f +/- %.1f um (bias %+0.1f)\n', ...
            label,nnz(good),nNoiseRealizations,mean(Rfits(good,method)),std(Rfits(good,method)), ...
            mean(Rfits(good,method))-R_true_um,mean(rfits(good,method)),std(rfits(good,method)), ...
            mean(rfits(good,method))-r_true_um);
    else
        fprintf('%s: NO SUCCESSFUL FITS\n',label);
    end
end

% Show one noisy signal and its predicted fits.
figure('Color','w','Name','129Xe single-voxel fitting comparison');
plot(bvalues,noisySignal,'ko','MarkerFaceColor','k','DisplayName','Rician noisy data'); hold on;
plot(bvalues,trueSignal,'k-','LineWidth',1.8,'DisplayName','Ground truth');
colors = {'b','r'}; labels = {'MLfitRr2','MLfitRr3'};
bfine = linspace(min(bvalues),max(bvalues),250);
for method=1:2
    x = squeeze(fitResults(2,method,:));
    if all(isfinite(x))
        yp = xeSignal(x(1),x(2),x(3),bfine,D0,deltaT);
        if all(isfinite(yp))
            plot(bfine,yp,colors{method},'LineWidth',1.5,'DisplayName',labels{method});
        end
    end
end
xlabel('b (s/cm^2)'); ylabel('Signal (normalized)');
title(sprintf('Single voxel: true R=%g um, r=%g um, SNR=%g',R_true_um,r_true_um,SNR_at_b0));
legend('Location','northeast'); grid on; hold off;

% Plot distributions for the same underlying voxel across noise realizations.
figure('Color','w','Name','129Xe fitting bias comparison');
tiledlayout(1,2);
nexttile; boxplot(Rfits,'Labels',{'Rr2','Rr3'}); hold on;
yline(R_true_um,'k--','True R'); ylabel('R (um)'); title('External radius'); grid on;
nexttile; boxplot(rfits,'Labels',{'Rr2','Rr3'}); hold on;
yline(r_true_um,'k--','True r'); ylabel('r (um)'); title('Internal radius'); grid on;

% Save numerical results for subsequent analysis.
results = struct('bvalues',bvalues,'D0',D0,'deltaT',deltaT,'R_true_um',R_true_um, ...
    'r_true_um',r_true_um,'S0_true',S0_true,'SNR_at_b0',SNR_at_b0, ...
    'trueSignal',trueSignal,'noisySignal',noisySignal,'fitResults',fitResults, ...
    'fitObjectives',fitObjectives,'Rfits_um',Rfits,'rfits_um',rfits,'valid',valid);
save('SingleVoxel_Xe_FitComparison.mat','results');
fprintf('\nSaved SingleVoxel_Xe_FitComparison.mat\n');

function printFit(name,x,f,Rtrue,rtrue)
    fprintf('%s: R=%.2f um (%+.2f), r=%.2f um (%+.2f), h=%.2f um, S0=%.4f, fval=%.4g\n', ...
        name,x(1)*1e4,x(1)*1e4-Rtrue,x(2)*1e4,x(2)*1e4-rtrue, ...
        (x(1)-x(2))*1e4,x(3),f);
end

function signal = xeSignal(R,r,S0,b,D0,deltaT)
    % Published 129Xe HUMAN equations used in both supplied fitting functions.
    L1 = sqrt(2*D0*deltaT); L2 = sqrt(4*D0*deltaT);
    q = 1-r/R;
    DL0 = D0*exp(-2.81*q^1.76);
    betaL = 21.2*(R/L1)^1.5*exp(-3.65*q^(-0.5));
    c1 = 1.13*(R/L2)-1.40*(R/L2)^2;
    c2 = 3.55-11.27*(R/L2)+7.44*(R/L2)^2;
    u = c1*q+c2*q^2;
    DT0 = D0*exp(-0.74*(L2/R)^1.47)*(1+u);
    ct0 = 1.89-5.82*(R/L2)+3.6*(R/L2)^2;
    ct1 = -4.21+14.1*(R/L2)-9.52*(R/L2)^2;
    ct2 = 2.14-7.35*(R/L2)+5.03*(R/L2)^2;
    betaT = ct0+ct1*(r/R)+ct2*(r/R)^2;
    DT = DT0*(1+betaT.*b*DT0);
    DL = DL0*(1-betaL.*b*DL0);
    if any(~isfinite(DT)) || any(~isfinite(DL)) || any(DL<=0) || any(DT<0) || any(DL<DT)
        signal = nan(size(b)); return;
    end
    z = b.*(DL-DT);
    orient = ones(size(z));
    ii = z>1e-8;
    orient(ii) = sqrt(pi)./(2*sqrt(z(ii))).*erf(sqrt(z(ii)));
    orient(~ii) = 1-z(~ii)/3+z(~ii).^2/10;
    signal = S0*exp(-b.*DT).*orient;
end
