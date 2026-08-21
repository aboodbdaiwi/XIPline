
function spectralDomainSignal = calcSpectralDomainSignal(obj,f)
% Calculates the spectral domain signal from the mix of NMR
% components at the given spectral frequencies (f is in Hz).
% Note, this function returns the net spectrum from the sum of
% all components in the NMR_Mix. The individual component
% spectrums can be obtained with calcComponentSpectralDomainSignal.
nFreqSamples = length(f);
nComp = length(obj.area);
spectralDomainSignal = zeros(nFreqSamples, 1);
for iComp=1:nComp
    spectralDomainSignal = spectralDomainSignal + ...
        obj.area(iComp)*exp(1i*pi/180*obj.phase(iComp))./...
        (1i*2*pi*(f-obj.freq(iComp))+pi*obj.fwhm(iComp));
end
end
