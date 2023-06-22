% This demo illustrates how the NMR_Mix and NMR_Fit objects work
%
% Author: Scott Haile Robertson
%

% Create a model of a 3 component system, providing the area (intensity),
% frequency, fwhm, and phase of each component
area = [10 10 10];      %arbs
freq = [400 30 -4E3];   %Hz
fwhm = [300 200 100];   %Hz
phase = [75 -40 0];     %deg
nmrMix = NMR_Mix(area, freq, fwhm, phase);

% Simulate the noiseless time domain signal
dwell_time = 62E-6; %sec
npts = 512;
t = (dwell_time*((1:npts) - 1))'; %sec
signal = nmrMix.calcTimeDomainSignal(t);

% Add some noise (note this should probably be white gaussian, but hey -
% its just a demo)
noiseMag = 2;
noisySignal = signal + noiseMag*([rand(size(signal)) + 1i*rand(size(signal))]-[0.5+0.5i]);

% Select line broadening and zeropadding values
zeroPadSize =npts*2^2; % This has no effect on fitting, only plotting/display
linebroadening = 0; %Hz

%% Fitting Option 1 (easiest, but least robust)
% Fit the noisy FID data in time domain with no guesses
% Note that this method often "misses" narrow peaks
disp('Performing automatic peak fitting for 3 peaks...')
nmrFit = NMR_TimeFit(noisySignal,t,[],[],[],[],linebroadening,zeroPadSize);
nmrFit = nmrFit.autoAddComponents(3);

% Describe and show the fit
disp('Fits without any guessses:');
nmrFit.describe();
figure();
%%nmrFit.plotFit();
nmrFit.plotSpectralFit
figure();
nmrFit.plotTimeFit();

%% Fitting Option 2 (works well when you have little prior knowledge)
% Use a interactive fitting tool to guide the fitting
disp('Fitting peaks with fitTool...')
nmrFit2 = NMR_TimeFit(noisySignal,t,[],[],[],[],linebroadening,zeroPadSize);
nmrFit2 = nmrFit2.fitTool();

% Describe and show the fit
disp('Fits using fitTool:');
nmrFit2.describe();
figure();
nmrFit2.plotFit();
figure();
nmrFit2.plotTimeFit();

%% Fitting Option 3 (works really well, but requires prior knowledge)
% Now lets assume we had decent guesses of the frequencies and fwhms, lets 
% re-fit with some decent starting guesses
disp('Using automated fitting with decent starting guesses...')
area_guess = [1 1 1]; 
freq_guess = [390 32 -4100];
fwhm_guesses = [320 193 104];
phase_guesses = [0 0 0];
nmrFit3 = NMR_TimeFit(noisySignal,t,...
    area_guess,freq_guess,fwhm_guesses,phase_guesses,...
    linebroadening,zeroPadSize);
nmrFit3 = nmrFit3.fitTimeDomainSignal();

% Describe and show the fit
disp('Fits with reasonable guessses:');
nmrFit3.describe();
figure();
nmrFit3.plotSpectralFit();
figure();
nmrFit3.plotTimeFit();

%% Fitting Option 4 (best, but requires lots of prior knowledge)
% We can also constrain the fits to only allow parameters to vary within
% a certain range. If the algorithm keeps missing peaks, you can constrain
% it and often get better results
disp('Using automated fitting with decent starting guesses...')

% Provide intensity (area under the curve) guesses, and positivity
% constraint
area_guess = [1 1 1]; 
area_lowerBounds = [0 0 0]; % Always positive
area_upperBounds = 1E10*[1 1 1]; % Make it finite, but huge

% Provide freq guesses and constraints
freq_guess = [390 32 -4100];
freq_lowerBounds = [250 -100 -4500];
freq_upperBounds = [500 100 -3000];

% Provide fwhm guesses and constraints
fwhm_guesses = [320 193 104];
fwhm_lowerBounds = [0 0 0];
fwhm_upperBounds = inf*[1 1 1];

% I tend to not constrain the phase (allow -inf to inf)
phase_guesses = [1 1 1];
phase_lowerBounds = -inf*[1 1 1];
phase_upperBounds = inf*[1 1 1];

% Create Fit object, then get a rough fit
nmrFit4 = NMR_TimeFit(noisySignal,t,...
    area_guess,freq_guess,fwhm_guesses,phase_guesses,...
    linebroadening,zeroPadSize);
nmrFit4 = nmrFit4.fitTimeDomainSignal();

% Set constraint bounds, then re-fit
nmrFit4 = nmrFit4.setBounds(area_lowerBounds,area_upperBounds,...
    freq_lowerBounds,freq_upperBounds,...
    fwhm_lowerBounds,fwhm_upperBounds,...
    phase_lowerBounds,phase_upperBounds);
nmrFit4 = nmrFit4.fitTimeDomainSignal();

% Describe and show the fit
disp('Fits with reasonable guessses and constraints:');
nmrFit4.describe();
figure()
nmrFit4.plotFit();
figure();
nmrFit4.plotTimeFit();
% Note that the these methods often agree well with high SNR data, but be 
% aware that the automatic fitting sometimes has trouble "finding" a peak