function [fit_area, fit_freq, fit_fwhm, fit_phase] = calcTimeDomainSignalFit(obj)
% Fits exponentially decaying components to the given time
% domain signal and provides the results without saving them to
% this object. To save the results to the NMR_Mix, use
% fitTimeDomainSignal

%    fitoptions = optimoptions(@lsqcurvefit);
%                         fitoptions.Display = 'iter-detailed';
%                         fitoptions.Display = 'final-detailed';
fitoptions.Display = 'off';
fitoptions.MaxIter = 10000;
fitoptions.TolFun=1E-900;
fitoptions.TolX = 1E-15;
fitoptions.FinDiffType = 'central';
fitoptions.Algorithm = 'trust-region-reflective';
fitoptions.MaxFunEvals = 5000;

% Put all components into a single matrix
guess = [obj.area.*exp(1i*pi*obj.phase/180); 1i*2*pi*obj.freq-pi*obj.fwhm];

fit_params = lsqcurvefit(@calcUnconstrainedTimeSig,guess,obj.t,...
    obj.timeDomainSignal,[],[],fitoptions);


% Separate out the components from the matrix
fit_area = abs(fit_params(1,:));
fit_freq = imag(fit_params(2,:))/(2*pi);
fit_fwhm = -real(fit_params(2,:))/pi;
fit_phase = angle(fit_params(1,:))*180/pi;
end

function complexSig = calcUnconstrainedTimeSig(nmr_params,t)
% A function used in fitting to allow constraints for complex
% fitting. 
nComp = size(nmr_params,2);
complexSig = zeros(size(t));
for iComp = 1:nComp
    complexSig = complexSig + nmr_params(1,iComp).*exp(nmr_params(2,iComp)*t);
end
end