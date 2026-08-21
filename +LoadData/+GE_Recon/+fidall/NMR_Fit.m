function obj = NMR_Fit(time_domain_signal, t, ...
    area, freq, fwhm, phase, line_broadening, zeroPadSize)
% Initialize structure for fitting.

obj.area = area(:)';
obj.freq = freq(:)';
obj.fwhm = fwhm(:)';
obj.phase = phase(:)';

obj.t = t(:); % make a column vector
nSamples = length(obj.t);

% Default line broadening is zero
obj.lineBroadening = line_broadening;
if(isempty(obj.lineBroadening))
    obj.lineBroadening = 0;
end

% Apply linebroadening
time_domain_signal = time_domain_signal(:); % make a column vector
obj.timeDomainSignal = time_domain_signal.*exp(-pi*obj.lineBroadening*obj.t); % Note that pi*fwhm = 1/t2star

% Default zero padding is none
obj.zeroPadSize = zeroPadSize;
if(isempty(obj.zeroPadSize))
    obj.zeroPadSize = nSamples;
end

% Calculate spectrum
dwell_time = (obj.t(2)-obj.t(1));
obj.spectralDomainSignal = dwell_time...
    *fftshift(fft(obj.timeDomainSignal,obj.zeroPadSize));

% Calculate freq samples
obj.f = linspace(-0.5,0.5,obj.zeroPadSize+1)/dwell_time;
obj.f = obj.f(1:(end-1)); % Take off last sample to have nSamples
obj.f = obj.f(:); % Make a column vector

end