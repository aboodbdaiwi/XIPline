% simple script to fit an desired number of peaks
% edit the fit guess section if you have trouble fitting.
% Skip as many as you like, and then average as many as you like
% 5/24/17 mess around to take 3-peak fit to 4 peaks
% 6/30/17 fix TE90 calculation
% 7/21/17 include Voigt fitting option
% 1/26/18 enable rmoving of initial points.

clc; clear all; close all;
nSkip = 2001; % number of fids to skip at beginning (1601 to get dyn. gas, 2001 for Dixon bonus)
nAvg = 1;    % for files with multiple fids, allow averaging
nPeak = 1;   % how many peaks to fit; only applies if dissolved flag == 0
dissolved = 1; % if set to 1, fit multi-peak and do multiple dissolved iterations with constraints, override nPeak
Voigt = 1; % if set to 1, use voigt profile instead of lorentzian
duke_spectro = 1 ;  % set to 1 only for basic duke spectroscopy (Soher) sequence to calculate n_fids properly. 
points_off = 0; % low BW spectra have some spurious high signal in first few FIDs

%% first read in the file and create the twix object
[file, path] = uigetfile('*.*', 'Select file');  % starts in current dir 
file_with_path = strcat(path, file);  % join path and filename to open
twix_obj = mapVBVD(file_with_path); 
filename = file(15:end-4); % return filename back to something short
filename = strrep(filename, '_', '-'); % replace underscores to avoid subscript problems

%% dealing with fields that may or may not be in header for all sequences
if(~isfield(twix_obj.hdr.Meas, 'flTransRefAmpl'))
      header.imagetwix_obj.hdr.Meas.flTransRefAmpl = 0;
end

if(~isfield(twix_obj.hdr.MeasYaps.sTXSPEC.aRFPULSE{1,1}, 'flAmplitude'))
      twix_obj.hdr.MeasYaps.sTXSPEC.aRFPULSE{1,1}.flAmplitude = 0;
end

if(~isfield(twix_obj.hdr.Meas, 'GradDelayTime'))
      twix_obj.hdr.Meas.GradDelayTime = 0;
end

if(~isfield(twix_obj.hdr.Config, 'NCha'))
      twix_obj.hdr.Config.NCha = 0;
end

if(~isfield(twix_obj.hdr.Config, 'NSlc'))
      twix_obj.hdr.Config.NSlc = 0;
end

%% read key header variables
Seq_name = twix_obj.hdr.Config.SequenceDescription;
weight = twix_obj.hdr.Dicom.flUsedPatientWeight;  %
ref_amp = twix_obj.hdr.Spice.TransmitterReferenceAmplitude; % Reference amplitude entered
freq = twix_obj.hdr.Dicom.lFrequency;  % seems to work in all sequences
TE = twix_obj.hdr.Meas.alTE(1);   
TR = twix_obj.hdr.Config.TR;
dwell_time = twix_obj.hdr.MeasYaps.sRXSPEC.alDwellTime{1,1};  % appears most robust and generic
rf_amp1 = twix_obj.hdr.MeasYaps.sTXSPEC.aRFPULSE{1,1}.flAmplitude; 
freq = twix_obj.hdr.Dicom.lFrequency;  % seems to work in all sequences
%fft_factor = twix_obj.hdr.MeasYaps.asCoilSelectMeas{1,1}.aFFT_SCALE{1,1}.flFactor;  % fft scale factor
grad_delay= twix_obj.hdr.Meas.GradDelayTime;  % a factor in radial images that seems to vary
num_coils = twix_obj.hdr.Config.NCha; % number of coil channels
num_slices = twix_obj.hdr.Config.NSlc  % number of slices
num_pts = twix_obj.image.sqzSize(1);  % num_pts is always first in array. Second can be coils or FIDs
%num_fids = twix_obj.image.sqzSize(end); % get info from image instead of hdr when possible. works for multi-coil, not multi-slice
num_fids = twix_obj.hdr.Config.NoOfFourierLines; % seems to be generic to both GRE and radial
if duke_spectro == 1
    num_fids = twix_obj.hdr.Config.NRepMeas; % this works for basic duke spectroscopy sequence. The above contains number of points for spectro
end % if
% for multislice GRE can use Config.ImageLines, NImagLins, PhaseEncodingLines, RawLn, NoOfFourierLines 
% for radial can use Config.NLinMeas, NProj, NoOfFourierLines, 

%% print out key header variable values
fprintf('\n\n');
fprintf('\tName = %s\n',filename');
fprintf('\tName = %s\n',Seq_name');
fprintf('\tWeight = %0.1f kg \n',weight);
fprintf('\tTE=%0.0f us\n',TE);
fprintf('\tTR=%0.0f us\n',TR);
fprintf('\tFrequency=%0.0f Hz\n',freq);
fprintf('\tvoltage=%0.1f V\n',rf_amp1);
%fprintf('\tgradient_delay=%0.3f\n',grad_delay);
fprintf('\tdwell_time=%0.0f ns\n',dwell_time);
fprintf('\tnum_coils=%0.0f\n',num_coils);
fprintf('\tnum_slices=%0.0f\n',num_slices);
fprintf('\tnPts=%0.0f\n',num_pts);
fprintf('\tnFrames=%0.0f\n',num_fids);

%% Adjust variables as needed
dwell_s = dwell_time*1E-9;  % convert from nanoseconds 
if nAvg> num_fids
    nAvg = num_fids
    a=sprintf('Not enough fids. Setting nAvg to %g.\n',num_fids);
    disp(a);
end % if

%% create the actual FID data and deal with dimensions and slices
twix_obj.image.flagIgnoreSeg = true; %Add this flag to ignore the extra radial dimensions. Doesn't hurt other files.
%theFID = squeeze(double(twix_obj.image()));  % gets rid of dimensions of length 1
theFID = squeeze(double(twix_obj.image.unsorted())); % this version takes away the slice dimension and just tacks them end to end
fprintf('\nAnalysis of %0.0f averaged FIDs, skipping first %0.0f frames\n',nAvg,nSkip);
Data = theFID(:,nSkip+1:nSkip+1+nAvg-1);  % skip initial and final desired frames
MeanData = mean(Data,2);
if points_off~=0
    fprintf('\nRemoving first %0.0f points from averaged FID\n',points_off);
    MeanData(1:points_off) = []; % remove points
    MeanData(num_pts-points_off+1:points_off)=0; % just replace end of FID with zeroes
end %if

%% Fit Mean Data
t = double((0:(length(MeanData)-1))*dwell_s');

if dissolved == 0   % simple spectrum
    guess = zeros(1,nPeak); % populate area, frequency, FWHM and phase guess with right numbers of zeroes
%    meanfitObj = NMR_TimeFit(MeanData,t,guess,guess,guess,guess,0,10000); % populate with zero guesses
%     meanfitObj = NMR_TimeFit(MeanData,t,[1 1 ],[-50 -100],[20 20],[0 0 0],0,10000); % provide some guesses
    meanfitObj = NMR_TimeFit(MeanData,t,2e-7,100,5,0,0,10000); % provide some guesses

    meanfitObj.fitTimeDomainSignal();
    figure();
    meanfitObj.plotTimeAndSpectralFit;
    meanfitObj.describe()

    %% Calculate and report final gas phase FID spectral SNR
    timeFit=meanfitObj.calcTimeDomainSignal(meanfitObj.t);  % calculate fitted data in time domain
    timeRes=timeFit - meanfitObj.timeDomainSignal;
    n25pct = round(length(timeRes)/4);
    std25 = std(timeRes(end-n25pct:end)); % take standard deviation of tail of residuals
    SNR=meanfitObj.area/std25; % will be a matrix for multiple peak fit
    
else % do the whole nine yards
    %meanfitObj = NMR_TimeFit(MeanData,t,1,0,30,0,0,10000); % y,x,area, freq, FWHM, phase, lb, zeropadsize
    %meanfitObj = NMR_TimeFit(MeanData,t,1e-6,15000,250,0,0,10000); % try this for single peak. Specify frequency
    %meanfitObj = NMR_TimeFit(MeanData,t,[1e-6 2e-4],[11000 0],[200 40],[0 0],0,10000); % try fitting to 2 peaks
    %meanfitObj = NMR_TimeFit(MeanData,t,[0 0],[0 0],[0 0],[0 0],0,10000); % try no guesses
    %meanfitObj = NMR_TimeFit(MeanData,t,[2e-4 3e-4 2e-4 1e-5],[28 500 653 7400],[280 300 350 44],[-125 147 100 79],0,10000); % try fitting to 2 peaks
    %meanfitObj = NMR_TimeFit(MeanData,t,[1 1 1],[0 -700 -7400],[240 240 40],[0 0 0],0,10000); % standard calibration guesses
    %meanfitObj = NMR_TimeFit(MeanData,t,[1 1 ],[0 -340],[20 30],[0 0 0],0,10000); % cyclohexane water
    if Voigt == 1
        fprintf('Using Voigt Lineshape (Unconstrained) to fit\n');
%         meanfitObj = NMR_TimeFit_v(MeanData,t,[1 1 1],[7400 6700 0],[250 200 30],[0 200 0],[0 0 0],0,length(t)); % For spectra where gas is on-resonance somehow
        meanfitObj = NMR_TimeFit_v(MeanData,t,[1 1 1],[0 -700  -7400],[250 200 30],[0 200 0],[0 0 0],0,10000); % zero fill a bit. Shouldn't change fits?
        meanfitObj.fitTimeDomainSignal();
        guessV=[meanfitObj.area' meanfitObj.freq' meanfitObj.fwhmL' meanfitObj.fwhmG' meanfitObj.phase'];  % modified for Voigt
        figure();
        meanfitObj.plotTimeAndSpectralFit;
        meanfitObj.describe()
        RbcBarr=meanfitObj.area(1)/meanfitObj.area(2);
        deltaPhase = meanfitObj.phase(2)-meanfitObj.phase(1); % RBC-barrier phase difference
        deltaPhase = mod(abs(deltaPhase),180); % deal with wrap around, but also negative phase
        fprintf('RBC:barrier = %4.3f, delta_phase = %4.1f\n',RbcBarr,deltaPhase);
        fprintf('\nConstrained Voigt Fit.\n');
%         f = warndlg('Update GuessV. The Click OK', 'User Input');  %      % PAUSE HERE TO UPDATE FIT GUESSES and Constrain
%         drawnow     % Necessary to print the message
%         waitfor(f);
        meanfitObj = NMR_TimeFit_v(MeanData,t,guessV(:,1)',guessV(:,2)',guessV(:,3)',guessV(:,4)', guessV(:,5)',0,10000); % populate with revised guess
        meanfitObj.setBounds(0.33*guessV(:,1)',3*guessV(:,1),...
            guessV(:,2)'-1,guessV(:,2)'+1,...
            0.9*guessV(:,3)',1.1*guessV(:,3),0.9*guessV(:,4)'-1,1.1*guessV(:,4)+1,...
            [-inf -inf -inf],[inf inf inf]);
        meanfitObj.fitTimeDomainSignal(); %this is the actual fitting.
        figure();
        meanfitObj.plotTimeAndSpectralFit;
        meanfitObj.describe()
        RbcBarr=meanfitObj.area(1)/meanfitObj.area(2);
        deltaPhase = meanfitObj.phase(2)-meanfitObj.phase(1); % RBC-barrier phase difference
        deltaPhase = mod(abs(deltaPhase),180); % deal with wrap around, but also negative phase
        fprintf('RBC:barrier = %4.3f, delta_phase = %4.1f\n',RbcBarr,deltaPhase);
    end %if
end % for % comment this if you want to do the lorentzian stuff

%     %% now just return to regular old Lorentzian fitting the rest of the way
%     fprintf('\nReturnn to Using Lorentzian Lineshapes \n');
%     meanfitObj = NMR_TimeFit(MeanData,t,[1 1 1],[0 -700  -7400],[240 240 40],[0 0 0],0,length(t));
%     meanfitObj.fitTimeDomainSignal();
%     figure();
%     meanfitObj.plotTimeAndSpectralFit;
%     meanfitObj.describe()
%     fitparams_array1=[meanfitObj.area' meanfitObj.freq' meanfitObj.fwhm' meanfitObj.phase']; % array of fit params in old matrix format
%     
%     %% re-fit with 2 barrier peaks Lorentzian
%     fprintf('\nFit Four Peaks, Unconstrained\n');
%     % set guesses equal to prior 3-peak fit
%     guess(:,1)=meanfitObj.area;
%     guess(:,2)=meanfitObj.freq;
%     guess(:,3)=meanfitObj.fwhm;
%     guess(:,4)=meanfitObj.phase;
%     guessNew=[guess(1:2,:);guess(2:3,:)]; % duplicate barrier row in new guess matrix
%     guessNew(2,1)=guessNew(2,1)/3; % barrier 1 area is about 1/3 of all barrier
%     guessNew(3,1)=guessNew(3,1)/2; % barrier 2 area is also smaller
%     guessNew(1,2)=guessNew(1,2)+34; % add 1ppm to RBC frequency
%     guessNew(2,2)=guessNew(2,2)+4*34; % add 4ppm for barrier 1 frequency
%     guessNew(3,2)=guessNew(3,2)-34; % subtract 1 ppm for barrier 2 frequency
%     %guessNew(1,3)=guessNew(1,3)*10/8; % increase RBC linewidth by typical
%     %guessNew(2,3)=guessNew(2,3)*14.6/8.6; % increase barrier1 width nearly 2x
%     %guessNew(3,3)=guessNew(3,3)*9.3/8.6; % increase barrier 2 width slightly
%     guessNew(1,4)=guessNew(1,4)+12; % add just a bit to RBC phase
%     guessNew(2,4)=guessNew(2,4)-140; % decrease barrier1 phase massively
%     guessNew(3,4)=guessNew(3,4)+20; % increase barrier 2 phase slightly
%     guessNew(:,4)=mod(guessNew(1:4),360); % fix phases if needed
%     
%     % pass new guesses and other parameters to Object. This is not the fit.
%     meanfitObj = NMR_TimeFit(MeanData,t,guessNew(:,1)',guessNew(:,2)',guessNew(:,3)',guessNew(:,4)',0,10000); % populate with revised guess
%     
%     meanfitObj.fitTimeDomainSignal(); %this is the actual fitting.
%     figure();
%     meanfitObj.plotTimeAndSpectralFit;
%     meanfitObj.describe()
%     %% create old fashioned matrix of fit parameters for first fit for copying
%     fitparams_array2=[meanfitObj.area' meanfitObj.freq' meanfitObj.fwhm' meanfitObj.phase']; 
%     
%     % one more time with guesses and bounds. First repopulate with guesses
%     fprintf('\nFit Four Peaks with Constraints\n');
%     guessNew1(:,1)=meanfitObj.area;
%     guessNew1(:,2)=meanfitObj.freq;
%     guessNew1(:,3)=meanfitObj.fwhm;
%     guessNew1(:,4)=meanfitObj.phase;
%     meanfitObj = NMR_TimeFit(MeanData,t,guessNew1(:,1)',guessNew1(:,2)',guessNew1(:,3)',guessNew1(:,4)',0,10000); % populate with revised guess
%     
%     % then set bounds (after passing guesses) to pass the error checking
%     % meanfitObj.setBounds([1e-6 1e-6 1e-6 1e-6],[inf inf inf inf],...
%     %             [-200 -700 -1000 -10000],[100 100 100 100],...
%     %             [20 20 20 20],[400 400 401 400],...
%     %             [-inf -inf -inf -inf],[inf inf inf inf]);
%     meanfitObj.setBounds(0.33*guessNew1(:,1)',3*guessNew1(:,1),...
%         guessNew1(:,2)'-100,[100 100 100 100],...
%         0.5*guessNew1(:,3)',2*guessNew1(:,3),...
%         [-inf -inf -inf -inf],[inf inf inf inf]);
%     
%     meanfitObj.fitTimeDomainSignal(); %this is the actual fitting.
%     figure();
%     meanfitObj.plotTimeAndSpectralFit;
%     meanfitObj.describe()
%     %% create old fashioned matrix of fit parameters for first fit for copying
%     fitparams_array3=[meanfitObj.area' meanfitObj.freq' meanfitObj.fwhm' meanfitObj.phase'];
%     
%     % now go back to 3 peak fit, but constrain RBC frequency to 4-peak value
%     fprintf('\nFit Three Peaks with Constrained RBC frequency\n');
%     guess(1,2)=meanfitObj.freq(1); %take previous guesses, but fix RBC
%     meanfitObj = NMR_TimeFit(MeanData,t,guess(:,1)',guess(:,2)',guess(:,3)',guess(:,4)',0,10000); % populate with revised guess
% %    meanfitObj = NMR_TimeFit(MeanData,t,[1 1 1],[-1.57 -686.85 -7367.8],[282 266 40],[0 0 0],0,10000); 
%     % Typical fit with RBC frequency constrained to not go negative
%     meanfitObj.setBounds([-inf -inf -inf],[inf inf inf],...
%         [(guess(1,2)'-10) -inf -inf],[(guess(1,2)'+10) inf inf],...
%         (guess(:,3)'-10),(guess(:,3)+10)',...
%         [-inf -inf -inf],[inf inf inf]);
%     % fit with RBC and Barrier phase constrained - to try to connect to Dixon Pipeline
% %     meanfitObj.setBounds([-inf -inf -inf],[inf inf inf],...
% %         [(guess(1,2)'-20) -inf -inf],[(guess(1,2)'+100) inf inf],...
% %         (guess(:,3)'-40),(guess(:,3)+40)',...
% %         [93 3 -inf],[97 7 inf]);
%     % fit with freqs constrained ±1 and widths constrained ±1. Not typical
% %     meanfitObj.setBounds([-inf -inf -inf],[inf inf inf],...
% %                 [guess(1,2)'-1 guess(2,2)'-1 guess(3,2)'-1],[guess(1,2)'+1 guess(2,2)'+1 guess(3,2)'+1],...
% %                 0.99*guess(:,3)',1.01*guess(:,3)',...
% %                 [-inf -inf -inf],[inf inf inf]);
%     
%     meanfitObj.fitTimeDomainSignal(); %this is the actual fitting.
%     figure();
%     meanfitObj.plotTimeAndSpectralFit;
%     meanfitObj.describe()
%     %% create old fashioned matrix of fit parameters for first fit for copying
%     fitparams_array4=[meanfitObj.area' meanfitObj.freq' meanfitObj.fwhm' meanfitObj.phase'];
% 
% end %if
% 
% %% Calculate spectral SNR
% timeFit=meanfitObj.calcTimeDomainSignal(meanfitObj.t);  % calculate fitted data in time domain
% timeRes=timeFit - meanfitObj.timeDomainSignal;
% n25pct = round(length(timeRes)/4);
% std25 = std(timeRes(end-n25pct:end)); % take standard deviation of tail of residuals
% SNR=meanfitObj.area/std25
% 
% % Optional Calculate TE90 (like for fat-water) if multiple peaks
% if nPeak<2
%     fprintf('Not enough peaks for TE90');
% else
%     %% Calculate and report TE90 NEW
%     deltaPhase = meanfitObj.phase(2)-meanfitObj.phase(1); % RBC-barrier phase diff in cal spectrum
%     deltaPhase = mod(deltaPhase,180); % deal with wrap around
%     deltaF = abs(meanfitObj.freq(2)-meanfitObj.freq(1)); % absolute RBC-barrier freq difference
%     deltaTe90= (90-deltaPhase)/(360*deltaF); % how far off are we?
%     te90 =TE + deltaTe90*1e6; % in usec
%     fprintf('TE90=%3.2f ms\n',te90/1000); % report TE90 in ms
%     if abs(deltaPhase)> 90
%         fprintf(2,'WARNING! Phi_cal =%3.0f%c!; Use min TE!\n',deltaPhase,char(176));
%     end %if
% end % if
