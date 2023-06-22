%% NMR_Fit
% A class that fits FID to a series of exponentially decaying components using a
% NMR_Mix object. This class knows how to fit with
% or without constraints.
%
classdef NMR_Fit < Calibration.NMR_Mix
    
    properties
        ub;
        lb;
        t;
        f;
        timeDomainSignal;
        spectralDomainSignal;
        zeroPadSize;
        lineBroadening;
        
        % 95% confidence interval
        ci_area; % 95% confidence interval of area
        ci_freq; % 95% confidence interval of frequency
        ci_fwhm; % 95% confidence interval of FWHM
        ci_phase; % 95% confidence interval of phase
    end
    
    methods
        function obj = applyGlobalPhaseShift(obj, phaseShiftInDegrees)
            % Call super class to change NMR_Mix object
            obj = applyGlobalPhaseShift@Calibration.NMR_Mix(obj,phaseShiftInDegrees);
            
            % Adjust time domain signal
            obj.timeDomainSignal = exp(pi*1i*phaseShiftInDegrees/180)*obj.timeDomainSignal;
            obj.spectralDomainSignal = exp(pi*1i*phaseShiftInDegrees/180)*obj.spectralDomainSignal;
        end
        
        function obj = NMR_Fit(time_domain_signal, t, ...
                area, freq, fwhm, phase, line_broadening, zeroPadSize)
            % Constructs a fitting object. Note - this does not perform the
            % fitting.
            
            % Construct NMR_Mix object
            obj = obj@Calibration.NMR_Mix(area,freq,fwhm,phase);
            
            obj.t = t(:); % make a column vector
            
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
                obj.zeroPadSize = length(obj.t);
            end
            
            % Calculate spectrum
            dwell_time = (obj.t(2)-obj.t(1));
            obj.spectralDomainSignal = dwell_time...
                *fftshift(fft(obj.timeDomainSignal));
                   
            % Calculate freq samples
            obj.f = linspace(-0.5,0.5,length(obj.spectralDomainSignal)+1)/dwell_time;
            obj.f = obj.f(1:(end-1)); % Take off last sample to have nSamples
            obj.f = obj.f(:); % Make a column vector
            
            %Delete confidence intervals if they exist
            obj.ci_area = [];
            obj.ci_freq = [];
            obj.ci_fwhm = [];
            obj.ci_phase = [];
        end
        
        function obj = resetComponents(obj,area, frequency, fwhm, phase, varargin)
            % Resets component(s) to NMR_Mix object, then resorts to keep
            % frequencies in order
            obj.area = area(:)';
            obj.freq = frequency(:)';
            obj.fwhm = fwhm(:)';
            obj.phase = phase(:)';
            
            % reset confidence intervals
            if(nargin > 5)
                obj.ci_area = varargin{1};
                obj.ci_freq = varargin{2};
                obj.ci_fwhm = varargin{3};
                obj.ci_phase = varargin{4};
            end
            
            % Resort object by frequencies
            obj = obj.sortByFreq();
        end
        
        function obj = removeComponents(obj, componentNumbers)
            % Removes the component(s) given their index number
            obj.area(componentNumbers) = [];
            obj.freq(componentNumbers) = [];
            obj.fwhm(componentNumbers) = [];
            obj.phase(componentNumbers) = [];
            
            %Delete confidence intervals if they exist
            obj.ci_area = [];
            obj.ci_freq = [];
            obj.ci_fwhm = [];
            obj.ci_phase = [];
        end
        
        function obj = sortByFreq(obj)
            % Sorts all components by thier frequencies in decending order
            [sortFreq sortIdx] = sort(obj.freq, 'descend');
            
            % Sort fits
            obj.area = obj.area(sortIdx);
            obj.freq = sortFreq;
            obj.fwhm = obj.fwhm(sortIdx);
            obj.phase = obj.phase(sortIdx);
            
            % Sort confidence intervals
            if(~isempty(obj.ci_area))
                obj.ci_area = obj.ci_area(sortIdx,:);
                obj.ci_freq = obj.ci_freq(sortIdx,:);
                obj.ci_fwhm = obj.ci_fwhm(sortIdx,:);
                obj.ci_phase = obj.ci_phase(sortIdx,:);
            end
        end
        
        function obj = addComponents(obj,area, frequency, fwhm, phase)
            % Adds component(s) to NMR_Mix object, then resorts to keep
            % frequencies in order
            obj.area = [obj.area area(:)'];
            obj.freq = [obj.freq frequency(:)'];
            obj.fwhm = [obj.fwhm fwhm(:)'];
            obj.phase = [obj.phase phase(:)'];
            
            %Delete confidence intervals if they exist
            obj.ci_area = [];
            obj.ci_freq = [];
            obj.ci_fwhm = [];
            obj.ci_phase = [];
            
            % Resort object by frequencies
            obj = obj.sortByFreq();
        end
        
        function obj = setBounds(obj, area_lb, area_ub, freq_lb, freq_ub,...
                fwhm_lb, fwhm_ub, phase_lb, phase_ub)
            % Sets the constraints for fitting this NMR_Fit object.
            % Adding additional components to the NMR_Fit will erase the
            % bounds, so only add bounds immediately before you fit the
            % signal..
            obj.ub = [area_ub(:)'; freq_ub(:)'; fwhm_ub(:)'; phase_ub(:)'];
            obj.lb = [area_lb(:)'; freq_lb(:)'; fwhm_lb(:)'; phase_lb(:)'];
            
            % Check that constraints are possible
            if(any(obj.ub < obj.lb))
                error('Impossible constraint set.');
            end
            % Force guesses to be within constraints
            if(~isempty(obj.area))
                if(any(obj.area < obj.lb(1,:)) | any(obj.area > obj.ub(1,:)))
                    lbbad = obj.area < obj.lb(1,:)
                    obj.area(lbbad) = obj.lb(1,lbbad);
                    ubbad = obj.area > obj.ub(1,:)
                    obj.area(ubbad) = obj.ub(1,ubbad);
                    warning('Bad area bound');
                end
                if(any(obj.freq < obj.lb(2,:)) | any(obj.freq > obj.ub(2,:)))
                    lbbad = obj.freq < obj.lb(2,:)
                    obj.freq(lbbad) = obj.lb(2,lbbad);
                    ubbad = obj.freq > obj.ub(2,:)
                    obj.freq(ubbad) = obj.ub(2,ubbad);
                    warning('Bad freq bound');
                end
                if(any(obj.fwhm < obj.lb(3,:)) | any(obj.fwhm > obj.ub(3,:)))
                    lbbad = obj.fwhm < obj.lb(3,:)
                    obj.fwhm(lbbad) = obj.lb(3,lbbad);
                    
                    ubbad = obj.fwhm > obj.ub(3,:)
                    obj.fwhm(ubbad) = obj.ub(3,ubbad);
                    warning('Bad fwhm bound');
                end
                if(any(obj.phase < obj.lb(4,:)) | any(obj.phase > obj.ub(4,:)))
                    lbbad = obj.phase < obj.lb(4,:)
                    obj.phase(lbbad) = obj.lb(4,lbbad);
                    ubbad = obj.phase > obj.ub(4,:)
                    obj.phase(ubbad) = obj.ub(4,ubbad);
                    warning('Bad phase bound');
                end
            end
        end
        
        function describe(obj)
            % This function will report on the fitted component
            % intensities, frequencies, linewidths, and phases. Note that
            % if line broadening is used, this method subtracts off the
            % line broadening from the fitted values, thus reporting them
            % as if there was no line broadening applied
%             if(~isempty(obj.ci_area))
%                 disp('95% confidence intervals are shown in parentheses');
%                 disp('area (arbs)                               Freq (Hz)                      Linewidth(Hz)                   Phase(degrees)');
%                 for iComp = 1:length(obj.area)
%                     disp([sprintf('%8.3e',obj.area(iComp)) ' (' sprintf('%8.3e',min(obj.ci_area(iComp,:))) ':' sprintf('%8.3e',max(obj.ci_area(iComp,:))) ')   ' ...
%                         sprintf('%+8.2f',obj.freq(iComp)) ' (' sprintf('%+8.2f',min(obj.ci_freq(iComp,:))) ':' sprintf('%+8.2f',max(obj.ci_freq(iComp,:))) ')   ' ...
%                         sprintf('%8.2f',obj.fwhm(iComp)-obj.lineBroadening) ' (' sprintf('%8.2f',min(obj.ci_fwhm(iComp,:))-obj.lineBroadening) ':' sprintf('%8.2f',max(obj.ci_fwhm(iComp,:))-obj.lineBroadening) ')   ' ...
%                         sprintf('%+9.2f',obj.phase(iComp))  ' (' sprintf('%+9.2f',min(obj.ci_phase(iComp,:))) ':' sprintf('%+9.2f',max(obj.ci_phase(iComp,:))) ')' ...
%                         ]);
%                 end
%             else
                disp('area (arbs)  Freq (Hz)  Linewidth(Hz)  Phase(degrees)');
                for iComp = 1:length(obj.area)
                    disp([sprintf('%8.3e',obj.area(iComp)) ' ' ...
                        sprintf('%+8.2f',obj.freq(iComp))  '  ' ...
                        sprintf('%8.2f',obj.fwhm(iComp)-obj.lineBroadening) '  ' ...
                        sprintf('%+9.2f',obj.phase(iComp))]);
                end
%             end
        end
    end
end
