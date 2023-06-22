%% NMR_MIX
% A class that represents a series of exponentially decaying components. 
% The class knows how to calculate time or spectral domain signals. The 
% key assumption is that all components experience perfectly exponential 
% decay. This can be modeled as:
%
% s(t) = area*exp(i*phase)*exp(-pi*fwhm*t)*exp(i*2*pi*freq*t)
% S(f) = (area*exp(i*phase))/(i*2*pi*(freq-f)-pi*fwhm)
%
% Author: Scott Haile Robertson (2015)
% Website: www.ScottHaileRobertson.com
%
classdef NMR_Mix_v < handle
    
    properties
        area;  % Intensity (area under curve) of each component (arbs)
        freq;  % Frequency of each component in Hz
        fwhm; % FWHM of each component in Hz
        fwhmG; % FWHM of Gaussian
        phase; % Phase of each component in degrees
    end
    
    methods
        function obj = NMR_Mix_v(area, freq, fwhm, fwhmG, phase)
            % Constructs an NMR_Mix object with the given time domain
            % signal  and optional knowledge of component intensities
            % (arbs), frequencies (in Hz), fwhms (in Hz), and phases (in
            % radians). 
            obj = obj.addComponents(area, freq, fwhm, fwhmG, phase);
        end
        
        function obj = resetComponents(obj,area, frequency, fwhm, fwhmG, phase, varargin)
            % Resets component(s) to NMR_Mix object, then resorts to keep
            % frequencies in order
            obj.area = area(:)';
            obj.freq = frequency(:)';
            obj.fwhm = fwhm(:)';
            obj.fwhmG = fwhmG(:)';
            obj.phase = phase(:)';
            
            % Resort object by frequencies
            obj = obj.sortByFreq();
        end
        
        function obj = addComponents(obj,area, frequency, fwhm, fwhmG, phase)
            % Adds component(s) to NMR_Mix object, then resorts to keep
            % frequencies in order
            obj.area = [obj.area area(:)'];
            obj.freq = [obj.freq frequency(:)'];
            obj.fwhm = [obj.fwhm fwhm(:)'];
            obj.fwhmG = [obj.fwhmG fwhmG(:)'];
            obj.phase = [obj.phase phase(:)'];
            
            % Resort object by frequencies
            obj = obj.sortByFreq();
        end
        
        function obj = removeComponents(obj, componentNumbers)
            % Removes the component(s) given their index number
            obj.area(componentNumbers) = [];
            obj.freq(componentNumbers) = [];
            obj.fwhm(componentNumbers) = [];
            obj.fwhmG(componentNumbers) = [];
            obj.phase(componentNumbers) = [];
        end
        
        function obj = sortByFreq(obj)
            % Sorts all components by thier frequencies in decending order
            [sortFreq sortIdx] = sort(obj.freq, 'descend');
            
            % Sort fits
            obj.area = obj.area(sortIdx);
            obj.freq = sortFreq;
            obj.fwhm = obj.fwhm(sortIdx);
            obj.fwhmG = obj.fwhmG(sortIdx);
            obj.phase = obj.phase(sortIdx);
        end
        
        function obj = applyGlobalPhaseShift(obj, phaseShiftInDegrees)
            % Phase to RBC peak at first TE
            obj.phase = obj.phase+phaseShiftInDegrees;
        end
        
        function componentTimeDomainSignal = calcComponentTimeDomainSignal(obj,t)
            % Calculates the time domain signal from the individual components
            % of the NMR mix at the given time points (t is in sec). Note, 
            % this function returns the time signal for each individual
            % component. The overall "mix" signal can be obtained with 
            % calcTimeDomainSignal
            nTimePoints = length(t);
            nComp = length(obj.area);
            componentTimeDomainSignal = zeros(nTimePoints,nComp);           
            for iComp = 1
                componentTimeDomainSignal(:,iComp) = obj.area(iComp)*...
                    exp(1i*pi/180*obj.phase(iComp) + 1i*2*pi*t*obj.freq(iComp))...
                    .* exp(-pi*t*obj.fwhm(iComp));
            end
            
            for iComp = 2:nComp
                componentTimeDomainSignal(:,iComp) = obj.area(iComp)* ...
                    exp(1i*pi/180*obj.phase(iComp) + 1i*2*pi*t*obj.freq(iComp))...
                    .* exp(-t.^2*4*log(2)*(obj.fwhmG(iComp))^2) .* exp(-pi*t*obj.fwhm(iComp));
            end
        end
        
        function timeDomainSignal = calcTimeDomainSignal(obj,t)
            % Calculates the net time domain signal from the mix of NMR
            % components at the given time points (t is in sec). Note, this
            % function returns the net time signal from the sum of all 
            % components in the NMR_Mix. The individual component signals
            % can be obtained with calcComponentTimeDomainSignal
            nTimePoints = length(t);
            nComp = length(obj.area);
            timeDomainSignal = zeros(nTimePoints, 1);           
            for iComp = 1
                timeDomainSignal = timeDomainSignal + ...
                    obj.area(iComp)* exp(1i*pi/180*obj.phase(iComp) + 1i*2*pi*t*obj.freq(iComp))...
                    .* exp(-pi*t*obj.fwhm(iComp));
            end
            
            for iComp = 2:nComp
                timeDomainSignal = timeDomainSignal + ...
                    obj.area(iComp)* exp(1i*pi/180*obj.phase(iComp) + 1i*2*pi*t*obj.freq(iComp))...
                    .* exp(-t.^2*4*log(2)*(obj.fwhmG(iComp))^2) .* exp(-pi*t*obj.fwhm(iComp));
            end
        end
        
%         function componentSpectralDomainSignal = calcComponentSpectralDomainSignal(obj,f)
%             % Calculates the spectral domain signal from the individual 
%             % components of the NMR_mix at the given spectral frequencies 
%             % (f is in Hz). Note, this function returns the spectrum for each individual
%             % component. The overall "mix" spectrum can be obtained with 
%             % calcSpectralDomainSignal.
%             nFreqSamples = length(f);
%             nComp = length(obj.area);
%             componentSpectralDomainSignal = zeros(nFreqSamples,nComp);
%             for iComp=1:nComp
%                 componentSpectralDomainSignal(:,iComp) = ...
%                     obj.area(iComp)*exp(1i*pi/180*obj.phase(iComp))./...
%                     (1i*2*pi*(f-obj.freq(iComp))+pi*obj.fwhm(iComp));
%             end
%         end
%         
%         function spectralDomainSignal = calcSpectralDomainSignal(obj,f)
%             % Calculates the spectral domain signal from the mix of NMR 
%             % components at the given spectral frequencies (f is in Hz). 
%             % Note, this function returns the net spectrum from the sum of
%             % all components in the NMR_Mix. The individual component 
%             % spectrums can be obtained with calcComponentSpectralDomainSignal.
%             nFreqSamples = length(f);
%             nComp = length(obj.area);
%             spectralDomainSignal = zeros(nFreqSamples, 1);
%             for iComp=1:nComp
%                 spectralDomainSignal = spectralDomainSignal + ...
%                     obj.area(iComp)*exp(1i*pi/180*obj.phase(iComp))./...
%                     (1i*2*pi*(f-obj.freq(iComp))+pi*obj.fwhm(iComp));
%             end
%         end
%         
%         function timeDomainComponentLorentzianSignal = calcTimeDomainComponentLorentzianSignal(obj,t)
%             % Calculates time domain signal of unphase lorenzian signals 
%             % of each NMR component at the given spectral frequencies 
%             % (f is in Hz). The lorentzians are just the real part of the 
%             % spectrum.
%             % Note, this function returns a curve for each individual
%             % component. The overall "mix" of time domain lorentzian curves 
%             % can be obtained with calcTimeDomainLorentzianSignal.           
%             nTimePoints = length(t);
%             nComp = length(obj.area);
%             timeDomainComponentLorentzianSignal = zeros(nTimePoints,nComp);
%             for iComp=1:nComp
%                 timeDomainComponentLorentzianSignal(:,iComp) = obj.area(iComp)* ...
%                     exp(-pi*t*obj.fwhm(iComp) + 1i*2*pi*t*obj.freq(iComp));
%             end
%         end
%         
        function componentLorentzianCurves = calcComponentLorentzianCurves(obj,f)
            % Calculates unphased lorenzian curves for each NMR component 
            % at the given spectral frequencies (f is in Hz). The 
            % lorentzians are just the real part of the spectrum.
            % Note, this function returns a curve for each individual
            % component. The overall "mix" of lorentzian curves can be
            % obtained with calcLorentzianCurves.
            
            t = linspace(0,1,length(f)+1)/(f(2)-f(1));
            t = t(1:end-1);
            nTimePoints = length(f);
            nComp = length(obj.area);            
            componentTimeDomainSignal = zeros(nTimePoints,nComp);
            
            for iComp=1:nComp
                componentTimeDomainSignal(:,iComp) = ...
                    obj.area(iComp)* exp(1i*2*pi*t*obj.freq(iComp))...
                    .* exp(-t.^2*4*log(2)*(obj.fwhmG(iComp))^2) .* exp(-pi*t*abs(obj.fwhm(iComp)));
            end
            
            componentLorentzianCurves = fftshift(fft(componentTimeDomainSignal),1);
            
        end
%         
%         function lorentzianCurves = calcLorentzianCurves(obj,f)
%             % Calculates unphased lorenzian curves for the net NMR_mix
%             % at the given spectral frequencies (f is in Hz). The 
%             % lorentzians are just the real part of the spectrum.
%             % Note, this function returns a curve for the net signal from 
%             % all components. The individual lorentzian curves can be
%             % obtained with calcComponentLorentzianCurves.
%             nFreqSamples = length(f);
%             nComp = length(obj.area);
%             lorentzianCurves = zeros(nFreqSamples,1);
%             for iComp=1:nComp
%                 lorentzianCurves = lorentzianCurves + ...
%                     (obj.area(iComp).*pi*obj.fwhm(iComp))./...
%                     ((2*pi*(obj.freq(iComp)-f)).^2+(pi*obj.fwhm(iComp))^2);
%             end
%         end
%         
%         function componentDispersiveCurves = calcComponentDispersiveCurves(obj,f)
%             % Calculates unphased dispersive curves for each NMR component 
%             % at the given spectral frequencies (f is in Hz). The 
%             % dispersive functions are just the imaginary part of the spectrum.
%             % Note, this function returns a curve for each individual
%             % component. The overall "mix" of dispersive functions can be
%             % obtained with calcDispersiveCurves.
%             nFreqSamples = length(f);
%             nComp = length(obj.area);
%             componentDispersiveCurves = zeros(nFreqSamples,nComp);
%             for iComp=1:nComp
%                 componentDispersiveCurves(:,iComp) = ...
%                     (obj.area(iComp).*(2*pi*(obj.freq(iComp)-f)))./...
%                     ((2*pi*(obj.freq(iComp)-f)).^2+(pi*obj.fwhm(iComp))^2);
%             end
%         end
%         
%         function dispersiveCurves = calcDispersiveCurves(obj,f)
%             % Calculates unphased dispersive functions for the net NMR_mix
%             % at the given spectral frequencies (f is in Hz). The 
%             % dispersive functions are just the imaginary part of the spectrum.
%             % Note, this function returns a curve for the net signal from 
%             % all components. The individual dispersive functions can be
%             % obtained with calcComponentDispersiveCurves.
%             nFreqSamples = length(f);
%             nComp = length(obj.area);
%             dispersiveCurves = zeros(nFreqSamples,1);
%             for iComp=1:nComp
%                 dispersiveCurves = dispersiveCurves + ...
%                     (obj.area(iComp).*(2*pi*(obj.freq(iComp)-f)))./...
%                     ((2*pi*(obj.freq(iComp)-f)).^2+(pi*obj.fwhm(iComp))^2);
%             end
%         end
        
        function estSNRdB = estimateSNR(obj)
            fitSignal = obj.calcTimeDomainSignal(obj.t);
            fitNoise = fitSignal - obj.timeDomainSignal;
            estSNRdB = snr(fitSignal, fitNoise);
        end
    end
end