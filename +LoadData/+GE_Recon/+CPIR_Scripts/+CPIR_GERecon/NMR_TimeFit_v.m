%% NMR_TimeFit
% A class that fits FID to a series of exponentially decaying components using a
% NMR_Mix object. This class knows how to fit to a time domain signal with
% or without constraints.
%
classdef NMR_TimeFit_v < NMR_Fit_v
    methods
        function obj = NMR_TimeFit_v(time_domain_signal, t, ...
                area, freq, fwhm, fwhmG, phase, line_broadening, zeroPadSize)
            % Construct an NMR_TimeFit
            obj = obj@NMR_Fit_v(time_domain_signal, t, ...
                area, freq, fwhm, fwhmG, phase, line_broadening, zeroPadSize);
        end
        function [fit_area, fit_freq, fit_fwhm, fit_fwhmG, fit_phase, ci_area, ...
                ci_freq, ci_fwhm, ci_fwhmG, ci_phase] = fitComponentToResidual(obj)
            % Fits a single exponentially decaying component to the
            % residual signal. This is useful in peak fitting. Note,
            % however, that this function only fits based on residuals,
            % thus you are fitting the residuals, and not the signal. After
            % fitting the residuals, it is recommended that you "refit" all
            % the components as they will likely change slightly to better
            % accomodate the new component.
            
            %% Step one:
            % Use the magnitude of the residual freq domain signal to
            % make a guess of the next peaks frequency
            
            % If there are already components, fit them first, then
            % subtract them to get the residual spectrum
            if(~isempty(obj.area))
                % Calculate spectrum at starting time
                dwell_time = (obj.t(2)-obj.t(1));
                
                % Calculate freq samples
                residualSpectrum = obj.spectralDomainSignal - ...
                    dwell_time*fftshift(fft(obj.calcTimeDomainSignal(obj.t)));
            else
                residualSpectrum = obj.spectralDomainSignal;
            end
            
            % Find a guess frequency for the next peak
            [maxVal maxIdx] = max(abs(residualSpectrum));
            peakFreq = obj.f(maxIdx);
            
            %% Step two:
            % Fit the residual signal in the time domain, giving reasonable
            % bounds
            freq_halfBound = 100;
            fwhm_UpperBound = 500;
            residualTimeDomainSignal = obj.timeDomainSignal - obj.calcTimeDomainSignal(obj.t);
            largestPeakFit = NMR_TimeFit(residualTimeDomainSignal, obj.t, ...
                1, peakFreq, 100, 0, obj.lineBroadening, obj.zeroPadSize);
            largestPeakFit = largestPeakFit.setBounds(0, inf, ...
                peakFreq-freq_halfBound, peakFreq+freq_halfBound,...
                0, fwhm_UpperBound,-inf,inf);
            [fit_area, fit_freq, fit_fwhm, fit_fwhmG, fit_phase, ci_area, ci_freq, ci_fwhm, ci_fwhmG, ci_phase, v] = ...
                largestPeakFit.calcTimeDomainSignalFit();
        end
        
        function obj = autoAddComponent(obj)
            % This function attemps to automatically fit and add a new
            % exponentially decaying signal to the NMR_mix object
            
            % Fit next component to residual signal, then add it to this
            % NMR_Mix object
            [add_area, add_freq, add_fwhm, add_phase junk_area junk_freq ...
                junk_fwhm junkfwhmG junk_phase] = obj.fitComponentToResidual();
            
            % Add fitted component to NMR_mix
            obj = obj.addComponents(add_area, add_freq, add_fwhm, add_fwhmG, add_phase);
            
            % Refit all components after addition of latest component
            obj = obj.fitTimeDomainSignal();
            
            if(~isempty(obj.ub))
                % Remove bounds in case sorted order changed
                obj.ub = [];
                obj.lb = [];
                warning('Deleting bounds - you need to call setBounds again with the new components bounds');
            end
        end
        
        function obj = autoAddComponents(obj, nComponents)
            % This function attemps to automatically fit and add
            % n (defined by nComponents) new exponentially decaying signals
            % to the NMR_mix object
            
            for iComp = 1:nComponents
                obj = obj.autoAddComponent();
            end
        end
        
        function obj = fitTool(obj)
            % This function attemps to adds exponentially decaying signal
            % components until the user stops the fit.
            
            % Create a figure to show results
            thisFigure = gcf();
            clf;
            
            % Keep fitting until user says stop
            continueFitting = true;
            while(continueFitting)
                % Create temporary NRM_Fix object and add newest component
                tempFit = obj.copyFitObject();
                tempFit = tempFit.autoAddComponent();
                
                % Show fit
                tempFit.plotTimeAndSpectralFit();
                
                % Report new fits
                tempFit.describe();
                
                % If user wants to keep component, add it to the object
                output_str1 = lower(input('Keep fit? [y/n]','s'));
                keepFit = length(findstr(output_str1, 'y')>0);
                if(keepFit) 
                    obj = tempFit.copyFitObject();
                    
                    % If user wants to fit more peaks, continue fitting...
                    output_str2 = lower(input('\tFit more peaks? [y/n]','s'));
                    continueFitting = length(findstr(output_str2, 'y')>0);
                else
                    % Assume user is done since they didnt like the last
                    % fit
                    continueFitting = 0;
                end
            end
        end
        
        function copyObj = copyFitObject(obj)
            % Undo linebroadening
            timeSig = obj.timeDomainSignal.*exp(pi*obj.lineBroadening*obj.t); 
            
            % Create object
            copyObj = NMR_TimeFit_v(timeSig, obj.t, ...
                obj.area, obj.freq, obj.fwhm, obj.fwhmG, obj.phase, obj.lineBroadening , obj.zeroPadSize);

            % Set bounds
            if(~isempty(obj.ub))
                copyObj.ub = obj.ub;
                copyObj.lb = obj.lb;
            end
        end
        
        function [fit_area, fit_freq, fit_fwhm, fit_fwhmG, fit_phase, ci_area, ci_freq, ci_fwhm, ci_fwhmG, ci_phase] = calcTimeDomainSignalFit(obj)
            % Fits exponentially decaying components to the given time
            % domain signal and provides the results without saving them to
            % this object. To save the results to the NMR_Mix, use
            % fitTimeDomainSignal
            
            fitoptions = optimoptions('lsqcurvefit');
            %                        fitoptions.Display = 'iter-detailed';
            %                         fitoptions.Display = 'final-detailed';
            fitoptions.Display = 'off';
            fitoptions.MaxIter = 10000;
            fitoptions.TolFun=1E-900;
            fitoptions.TolX = 1E-10;
            fitoptions.FinDiffType = 'central';
            fitoptions.Algorithm = 'trust-region-reflective'; % Not used 'levenberg-marquardt'
%             fitoptions.Algorithm = 'levenberg-marquardt';
            fitoptions.MaxFunEvals = 30000;
            
            % Put all components into a single matrix
            guess = [obj.area; obj.freq; obj.fwhm; ...
                obj.fwhmG; obj.phase];
            
            [fit_params,resnorm,residual,exitflag,output,lambda,J] = lsqcurvefit(@obj.calcConstrainedTimeSig,guess,obj.t,...
                [real(obj.timeDomainSignal),imag(obj.timeDomainSignal)],...
                obj.lb,obj.ub,fitoptions);

            % Separate out the components from the matrix
            fit_vec = fit_params(1,:).*exp(1i*pi*fit_params(5,:)/180);
            fit_area = abs(fit_vec); % Prevents neg amp
            fit_freq = fit_params(2,:);
            fit_fwhm = fit_params(3,:);
            fit_fwhmG = fit_params(4,:);
            fit_phase = atan2(imag(fit_vec), real(fit_vec))*180/pi; % Handles neg amp
            
            %% Check that no aliased frequencies were detected, if so correct them
            maxBW = max(obj.f)-min(obj.f);
            halfBW = 0.5*maxBW;
            aliasedIdx = find(abs(fit_freq)>halfBW);
            nAliased = length(aliasedIdx);
            while(nAliased > 0)
                for iAliased =1:nAliased
                    idxAliased = aliasedIdx(iAliased);
                    while(fit_freq(idxAliased)<-halfBW)
                        fit_freq(idxAliased) = fit_freq(idxAliased)+ maxBW;
                    end
                    while(fit_freq(idxAliased)>halfBW)
                        fit_freq(idxAliased) = fit_freq(idxAliased)- maxBW;
                    end
                end
                
                % Fit again
                guess = [fit_area; fit_freq; fit_fwhm;  fit_fwhmG; fit_phase];
                [fit_params,resnorm,residual,exitflag,output,lambda,J] = lsqcurvefit(@obj.calcConstrainedTimeSig,guess,obj.t,...
                    [real(obj.timeDomainSignal),imag(obj.timeDomainSignal)],...
                    obj.lb,obj.ub,fitoptions);

                % Separate out the components from the matrix
                fit_vec = fit_params(1,:).*exp(1i*pi*fit_params(5,:)/180);
                fit_area = abs(fit_vec); % Prevents neg amp
                fit_freq = fit_params(2,:);
                fit_fwhm = fit_params(3,:);
                fit_fwhmG = fit_params(4,:);
                fit_phase = atan2(imag(fit_vec), real(fit_vec))*180/pi; % Handles neg amp
                
                aliasedIdx = find(abs(fit_freq)>halfBW);
                nAliased = length(aliasedIdx);
            end

             % Calculate 95% confidence interval using jacobian
%             ci = nlparci(fit_params,residual,'jacobian',J);
            ci = ones(length(fit_area)*5,2);
            ci = reshape(ci,[size(fit_params) 2]);
            ci_vec = squeeze(ci(1,:,:).*exp(1i*pi*ci(5,:,:)/180));
            ci_area = [abs(ci_vec)];
            ci_freq = squeeze(ci(2,:,:));
            ci_fwhm = squeeze(ci(3,:,:));
            ci_fwhmG = squeeze(ci(4,:,:));
            ci_phase = angle(ci_vec)*180/pi;            
        end
        
        function obj = fitTimeDomainSignal(obj)
            % Fits exponentially decaying components to the given time
            % domain signal and saves the results to this NMR_Mix object.
            % To just return the fits and not save the results, use
            % calcTimeDomainSignalFit
            [fit_area, fit_freq, fit_fwhm, fit_fwhmG, fit_phase, ci_area, ci_freq, ci_fwhm, ci_fwhmG, ci_phase] = obj.calcTimeDomainSignalFit();
            
            % Save fits
            obj = obj.resetComponents(fit_area, fit_freq, fit_fwhm, fit_fwhmG, fit_phase, ...
                ci_area, ci_freq, ci_fwhm, ci_fwhmG, ci_phase);
        end
        
        function complexSig = calcUnconstrainedTimeSig(obj,nmr_params,t)
            % A function used in fitting to allow constraints for complex
            % fitting. This is the same as calling calcTimeDomainSignal of
            % the NMR_Mix object, only it separates the real and imaginary
            % parts to allow for constraints to be used in fitting.
            nComp = size(nmr_params,2);
            complexSig = zeros(size(t));
            for iComp = 1:nComp
                complexSig = complexSig + nmr_params(1,iComp).*exp(nmr_params(2,iComp)*t);
            end
        end
        
        function realImagSig = calcConstrainedTimeSig(obj,nmr_params,t)
            % A function used in fitting to allow constraints for complex
            % fitting. This is the same as calling calcTimeDomainSignal of
            % the NMR_Mix object, only it separates the real and imaginary
            % parts to allow for constraints to be used in fitting.
            nComp = numel(nmr_params)/5;
            nmr_params = reshape(nmr_params,[5 nComp]);
            tmpNmrMix = NMR_Mix_v(nmr_params(1,:), nmr_params(2,:), ...
                nmr_params(3,:), nmr_params(4,:), nmr_params(5,:));
            complexSig = tmpNmrMix.calcTimeDomainSignal(t);
            realImagSig = [real(complexSig) imag(complexSig)];
        end
        
        function ax1 = plotTimeAndSpectralFit(obj)
            % Calculate fitted and residual spectrums usign time domain
            % signal so that amplitudes are correct even if signal is
            % truncated at the end
            dwell_time = (obj.t(2)-obj.t(1));
            zeroPaddedTime = min(obj.t(:)) + dwell_time*((1:obj.zeroPadSize)-1)';
            
            % Calculate spectrum
            zeroPaddedFreq = linspace(-0.5,0.5,obj.zeroPadSize+1)/dwell_time;
            zeroPaddedFreq = zeroPaddedFreq(1:(end-1)); % Take off last sample to have nSamples
            zeroPaddedFreq = zeroPaddedFreq(:); % Make a column vector
            
            fittedZeroPaddedSpectrum = dwell_time*fftshift(fft(obj.calcTimeDomainSignal(zeroPaddedTime)));
            fittedSpectrum = dwell_time*fftshift(fft(obj.calcTimeDomainSignal(obj.t)));
            individualSpectrums = dwell_time*fftshift(fft(obj.calcComponentTimeDomainSignal(zeroPaddedTime),[],1),1);
            residualSpectrum = obj.spectralDomainSignal - fittedSpectrum;
            
            % Calculate fitted and residual spectrums
            individualSignals = obj.calcComponentTimeDomainSignal(zeroPaddedTime);
            fittedZeroPaddedSignal = obj.calcTimeDomainSignal(zeroPaddedTime);
            fittedSignal = obj.calcTimeDomainSignal(obj.t);
            residualSignal = obj.timeDomainSignal - fittedSignal;
            
            
            nComponents = length(obj.area);
            fMat = repmat(zeroPaddedFreq,[1 nComponents]);
            tMat = repmat(zeroPaddedTime,[1 nComponents]);
            
            legendStrings = cell(1, nComponents);
            for iComp=1:nComponents
                legendStrings{iComp} = ['C\_' sprintf('%03.0f',iComp)];
            end
            
            % Show results to user
            ax3 = subplot(4,2,3);
            plot(obj.f,abs(obj.spectralDomainSignal),'.k','markersize',16);
            hold on;
            plot(zeroPaddedFreq,abs(fittedZeroPaddedSpectrum),'-g','Linewidth',2);
            plot(obj.f,abs(residualSpectrum),'.r','markersize',8);
            hold off;
            ylabel('Magnitude Intensity');
            set(ax3,'xticklabel',{[]}) ;
            set(ax3,'XDir','reverse');
            
            ax5 = subplot(4,2,5);
            plot(obj.f,real(obj.spectralDomainSignal),'.k','markersize',16);
            hold on;
            plot(zeroPaddedFreq,real(fittedZeroPaddedSpectrum),'-g','Linewidth',2);
            plot(obj.f,real(residualSpectrum),'.r','markersize',8);
            hold off;
            ylabel('Real Intensity');
            set(ax5,'xticklabel',{[]});
            set(ax5,'XDir','reverse');
            
            ax7 = subplot(4,2,7);
            plot(obj.f,imag(obj.spectralDomainSignal),'.k','markersize',16);
            hold on;
            plot(zeroPaddedFreq,imag(fittedZeroPaddedSpectrum),'-g','Linewidth',2);
            plot(obj.f,imag(residualSpectrum),'.r','markersize',8);
            hold off;
            xlabel('Spectral Frequency (Hz)');
            ylabel('Imaginary Intensity');
%             legend('Measured','Fitted','Residual');
            set(ax7,'XDir','reverse');
            
            ax1 = subplot(4,2,1);
            plot(fMat,real(individualSpectrums),'Linewidth',2);
%             legend(legendStrings);
            ylabel('Component Intensity');
            set(ax1,'xticklabel',{[]}) ;
            set(ax1,'XDir','reverse');
            
            % Keep all x axes in sinc
            linkaxes([ax1,ax3,ax5 ax7],'x');   
            ax1.XLim = [-6000,2000];
            
            % Show results to user
            ax4 = subplot(4,2,4);
            plot(obj.t,abs(obj.timeDomainSignal),'.k','markersize',16);
            hold on;
            plot(zeroPaddedTime,abs(fittedZeroPaddedSignal),'-g','Linewidth',2);
            plot(obj.t,abs(residualSignal),'.r','markersize',8);
            hold off;
            ylabel('Magnitude Intensity');
            set(ax4,'xticklabel',{[]}) ;
            
            ax6 = subplot(4,2,6);
            plot(obj.t,real(obj.timeDomainSignal),'.k','markersize',16);
            hold on;
            plot(zeroPaddedTime,real(fittedZeroPaddedSignal),'-g','Linewidth',2);
            plot(obj.t,real(residualSignal),'.r','markersize',8);
            hold off;
            ylabel('Real Intensity');
            set(ax6,'xticklabel',{[]});
            
            ax8 = subplot(4,2,8);
            plot(obj.t,imag(obj.timeDomainSignal),'.k','markersize',16);
            hold on;
            plot(zeroPaddedTime,imag(fittedZeroPaddedSignal),'-g','Linewidth',2);
            plot(obj.t,imag(residualSignal),'.r','markersize',8);
            hold off;
            xlabel('Time');
            ylabel('Imaginary Intensity');
            
            ax2 = subplot(4,2,2);
            plot(tMat,real(individualSignals));
            ylabel('Component Intensity');
            set(ax2,'xticklabel',{[]}) ;
            
            % Keep all x axes in sinc
            linkaxes([ax2,ax4,ax6 ax8],'x');
        end
        
        function ax1 = plotSpectralFit(obj)
            % Calculate fitted and residual spectrums usign time domain
            % signal so that amplitudes are correct even if signal is
            % truncated at the end
            dwell_time = (obj.t(2)-obj.t(1));
            zeroPaddedTime = min(obj.t(:)) + dwell_time*((1:obj.zeroPadSize)-1)';
            
            % Calculate spectrum
            zeroPaddedFreq = linspace(-0.5,0.5,obj.zeroPadSize+1)/dwell_time;
            zeroPaddedFreq = zeroPaddedFreq(1:(end-1)); % Take off last sample to have nSamples
            zeroPaddedFreq = zeroPaddedFreq(:); % Make a column vector
            
            fittedZeroPaddedSpectrum = dwell_time*fftshift(fft(obj.calcTimeDomainSignal(zeroPaddedTime)));
            fittedSpectrum = dwell_time*fftshift(fft(obj.calcTimeDomainSignal(obj.t)));
            individualSpectrums = dwell_time*fftshift(fft(obj.calcComponentTimeDomainSignal(zeroPaddedTime),[],1),1);
            residualSpectrum = obj.spectralDomainSignal - fittedSpectrum;
            
            nComponents = length(obj.area);
            fMat = repmat(zeroPaddedFreq,[1 nComponents]);
            
            legendStrings = cell(1, nComponents);
            for iComp=1:nComponents
                legendStrings{iComp} = ['C\_' sprintf('%03.0f',iComp)];
            end
            
            % Show results to user
            ax2 = subplot(4,1,2);
            plot(obj.f,abs(obj.spectralDomainSignal),'.k','markersize',16);
            hold on;
            plot(zeroPaddedFreq,abs(fittedZeroPaddedSpectrum),'-g','Linewidth',2);
            plot(obj.f,abs(residualSpectrum),'.r','markersize',8);
            hold off;
            ylabel('Magnitude Intensity');
            set(ax2,'xticklabel',{[]}) ;
            set(ax2,'XDir','reverse');
            
            ax4 = subplot(4,1,3);
            plot(obj.f,real(obj.spectralDomainSignal),'.k','markersize',16);
            hold on;
            plot(zeroPaddedFreq,real(fittedZeroPaddedSpectrum),'-g','Linewidth',2);
            plot(obj.f,real(residualSpectrum),'.r','markersize',8);
            hold off;
            ylabel('Real Intensity');
            set(ax4,'xticklabel',{[]});
            set(ax4,'XDir','reverse');
            
            ax5 = subplot(4,1,4);
            plot(obj.f,imag(obj.spectralDomainSignal),'.k','markersize',16);
            hold on;
            plot(zeroPaddedFreq,imag(fittedZeroPaddedSpectrum),'-g','Linewidth',2);
            plot(obj.f,imag(residualSpectrum),'.r','markersize',8);
            hold off;
            xlabel('Spectral Frequency (Hz)');
            ylabel('Imaginary Intensity');
            legend('Measured','Fitted','Residual');
            set(ax5,'XDir','reverse');
            
           
            ax1 = subplot(4,1,1);
            plot(fMat,real(individualSpectrums),'Linewidth',2);
%             legend(legendStrings);
            ylabel('Component Intensity');
            set(ax1,'xticklabel',{[]}) ;
            set(ax1,'XDir','reverse');
            
            % Keep all x axes in sinc
            linkaxes([ax1,ax2,ax4 ax5],'x');   
        end
        
        function ax1 = plotTimeFit(obj)
            % Calculate fitted and residual spectrums
            dwell_time = (obj.t(2)-obj.t(1));
            zeroPaddedTime = min(obj.t(:)) + dwell_time*((1:obj.zeroPadSize)-1)';
            individualSignals = obj.calcComponentTimeDomainSignal(zeroPaddedTime);
            fittedZeroPaddedSignal = obj.calcTimeDomainSignal(zeroPaddedTime);
            fittedSignal = obj.calcTimeDomainSignal(obj.t);
            residualSignal = obj.timeDomainSignal - fittedSignal;
            
            % Calculate lorentzian curves for each component
            nComponents = length(obj.area);
            tMat = repmat(zeroPaddedTime,[1 nComponents]);
            
            legendStrings = cell(1, nComponents);
            for iComp=1:nComponents
                legendStrings{iComp} = ['C\_' sprintf('%03.0f',iComp)];
            end
            
            % Show results to user
            ax2 = subplot(4,1,2);
            plot(obj.t,abs(obj.timeDomainSignal),'.k','markersize',16);
            hold on;
            plot(zeroPaddedTime,abs(fittedZeroPaddedSignal),'-g','Linewidth',2);
            plot(obj.t,abs(residualSignal),'.r','markersize',8);
            hold off;
            ylabel('Magnitude Intensity');
            set(ax2,'xticklabel',{[]}) ;
            
            %             ax3 = subplot(5,1,3);
            %             phaseSig = angle(zeroPaddedTimeimeDomainSignal);
            %             phaseFit = angle(fittedSignal);
            %             phaseDelta = phaseSig - phaseFit;
            %             plot(zeroPaddedTime, phaseSig,'-b');
            %             hold on;
            %             plot(zeroPaddedTime,phaseFit,'-g');
            %             plot(zeroPaddedTime,phaseDelta,'-r');
            %
            %             hold off;
            %             ylabel('Phase (Radians)');
            %             set(ax3,'xticklabel',{[]}) ;
            
            ax4 = subplot(4,1,3);
            plot(obj.t,real(obj.timeDomainSignal),'.k','markersize',16);
            hold on;
            plot(zeroPaddedTime,real(fittedZeroPaddedSignal),'-g','Linewidth',2);
            plot(obj.t,real(residualSignal),'.r','markersize',8);
            hold off;
            ylabel('Real Intensity');
            set(ax4,'xticklabel',{[]});
            
            ax5 = subplot(4,1,4);
            plot(obj.t,imag(obj.timeDomainSignal),'.k','markersize',16);
            hold on;
            plot(zeroPaddedTime,imag(fittedZeroPaddedSignal),'-g','Linewidth',2);
            plot(obj.t,imag(residualSignal),'.r','markersize',8);
            hold off;
            xlabel('Time');
            ylabel('Imaginary Intensity');
            legend('Measured','Fitted','Residual');
            
            %             if(~isempty(obj.fref))
            %                 % Add PPM axis
            %                 ax1ppm = subplot(5,1,1);
            %                 set(ax1ppm,'units','normalized',...
            %                     'XAxisLocation','top','YAxisLocation','right',...
            %                     'YTick',[],'YTickLabel',[],'Color','none');
            %
            %                 ax1 = axes('Position',get(ax1ppm,'Position'));
            %             else
            ax1 = subplot(4,1,1);
            %             end
            
            plot(tMat,real(individualSignals));
            legend(legendStrings);
            ylabel('Component Intensity');
            set(ax1,'xticklabel',{[]}) ;
            %             if(~isempty(obj.fref))
            %                 set(ax1ppm,'XDir','reverse');
            %             end
            
            % Keep all x axes in sinc
            linkaxes([ax1,ax2,ax4 ax5],'x');
            
            %             if(~isempty(obj.fref))
            %                 % Initialize XLim in correct units
            %                 set(ax1ppm,'xlim',NMR_Mix.relFreqToPpm(get(ax2,'XLim'),obj.fref));
            %
            %                 % Keep ppm axiz in sinc with freq axis
            %                 xLimListener = addlistener( ax1, 'XLim', 'PostSet', ...
            %                     @(src,evt) set(ax1ppm,'XLim',...
            %                     NMR_Mix.relFreqToPpm(get(ax1,'XLim'),obj.fref)) );
            %             end
        end
    end
end
