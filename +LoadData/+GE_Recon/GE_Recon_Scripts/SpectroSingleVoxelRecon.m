function [ plottedSpectrum, quantitationResults ] = SpectroSingleVoxelRecon( pfileFullPath )
%% SpectroSingleVoxelRecon
%
% Copyright (c) 2019, 2025, GE HealthCare
% All Rights Reserved
% 
% This material is proprietary to GE HealthCare. The methods and
% techniques described herein are considered trade secrets and/or 
% confidential. Reproduction or distribution, in whole or in part, is
% forbidden except by express written permission of GE HealthCare.
% GE is a trademark of General Electric Company used under trademark license.
% Resulting outputs are not for diagnostic purposes.
%
%   Reconstruct single voxel spectroscopy data using the GE Single Voxel
%   Reconstruction algorithm. The reconstructed spectrum is plotted into a
%   square image matrix. The spectrum is quantitated if quantitation
%   is enabled. The plotted ppm range can be adjusted using rhuser38 and
%   rhuser39 as described below.
%
%   plottedSpectrum = SpectroSingleVoxelRecon(pfileFullPath);


    %% Load pfile
    % Load pfile and hold on to the pfile directory. The pfile directory
    % will be used to save the final plotted spectrum dicom image.
    [pfilePath] = fileparts(pfileFullPath);
    pfileHandle = GERecon('Pfile.Load', pfileFullPath);
    header = GERecon('Pfile.Header');
    
    %% Extract Parameters of Interest
    % Extract parameters describing the single voxel acquisition and
    % reconstruction from the currently loaded pfile. These parameters
    % include product  processing options for the water-suppressed signal 
    % and non-water-suppressed reference FIDs.
    % Initialize space for channel combined data to zero.
    params = GERecon('Spectro.SingleVoxel.Parameters');  
    
    channelCombinedReference = zeros(params.acquiredFidLength, 1);
    channelCombinedSignal = zeros(params.acquiredFidLength, 1);
    channelImagesFigureHandle = figure;
    accumulatedChannelWeight = 0.0;
    for channel=1:pfileHandle.channels     
        
        %% Extract raw data
        % Single voxel scans do not use the slice or echo dimension; thus,
        % set both the slice and echo indices to one. Extract all views 
        % for the current channel. Each view is one acquired FID and is 
        % either a water-suppressed signal FID or a non-water-suppressed 
        % reference FID. The signal and reference FIDs are extracted 
        % separately below. The non-water-suppressed reference FIDs are 
        % always ordered in the pfile before the water-suppressed signal 
        % FIDs for product scans. The first dimension of the raw data is
        % the time index and the second dimension is the view (or FID)
        % index.
        % The data in the pfile has been accumulated based on how many
        % NEX's were prescribed. However, the data has not been scaled back
        % down based on the number of NEX's. For Single Voxel Spectroscopy
        % scans, the number of NEX's is indicated by header.RawHeader.navs;
        % thus, use this value to scale data according to the number of
        % NEX's.
        % Note that the GERecon utility can also extract raw data in the
        % Spectro Prescan files using the Pfile.ViewData function.
        slice = 1;
        echo = 1;
        allChannelData = GERecon('Pfile.KSpace', slice, echo, channel);        
        allChannelData = allChannelData / header.RawHeader.navs;
                        
        %% Process signal FIDs
        % Process each water-suppressed signal FID that was acquired for 
        % this channel and combine them into a single water-suppressed 
        % signal FID. If rf chopping was enabled then every other FID is
        % negated.
        % Note that water-suppressed-signal FIDs are always ordered in the
        % pfile after the non-water-suppressed reference FIDs.
        % The processing steps applied to each water-suppressed signal FID
        % are specified in the params structure. A description of each
        % processing step can be found in the ProcessFid function below.
        signalFids = allChannelData(:, (params.numNonWaterSuppressedReferenceFids+1):end);
        combinedSignalFid = zeros(params.acquiredFidLength, 1);
        for fidIndex=1:params.numWaterSuppressedSignalFids
            processedSignalFid = ProcessFid(signalFids(:,fidIndex), params.signalFrameProcessingOptions);    

            if(params.chop && (mod(fidIndex,2) ~= 1))
                processedSignalFid = processedSignalFid * -1.0;
            end
            
            combinedSignalFid = combinedSignalFid + (processedSignalFid ./ params.numWaterSuppressedSignalFids);
        end
                
        %% Process reference FIDs
        % If non-water-suppressed reference FIDs were acquired for this
        % channel, process them and combine them into a single
        % non-water-suppressed reference FID. The combined non-water
        % suppressed reference FID is processed one more time and the phase
        % correction vector applied to the combined non-water-suppressed
        % reference FID is also applied to the combined water-suppressed
        % signal FID. The phase correction vector contains the same phase
        % shifts applied to the combinedReferenceFid. Thus, when multiplied
        % with the signal FID, these phase shifts are also applied to the
        % signal FID.
        if(params.numNonWaterSuppressedReferenceFids > 0)
            
            referenceFids = allChannelData(:, 1:params.numNonWaterSuppressedReferenceFids);
            
            combinedReferenceFid = zeros(params.acquiredFidLength, 1);
            for fidIndex=1:params.numNonWaterSuppressedReferenceFids
                processedReferenceFid = ProcessFid(referenceFids(:,fidIndex), params.referenceFrameProcessingOptions);

                if(params.chop && (mod(fidIndex,2) ~= 1))
                    processedReferenceFid = processedReferenceFid * -1.0;
                end

                combinedReferenceFid = combinedReferenceFid + (processedReferenceFid ./ params.numNonWaterSuppressedReferenceFids);
            end                      
                        
            % Compute magnitude squared weight
            weight = max(abs(combinedReferenceFid)) * max(abs(combinedReferenceFid));
            
            [combinedReferenceFid, referenceCorrectionVector] = ProcessFid(combinedReferenceFid, params.combinedReferenceFrameProcessingOptions);
            
            combinedSignalFid = combinedSignalFid .* referenceCorrectionVector;            
        else            
            % Compute magnitude squared weight
            weight = max(abs(combinedSignalFid)) * max(abs(combinedSignalFid));
            combinedReferenceFid = 0;
        end
               
        %% Display Intermediate Channel Images
        % Perform water subtraction (if enabled) and transform the current
        % channel's data. Plot and display the transformed data.
        % Intermediate channel data plotting can be disabled by setting
        % displayChannelImages = 0.        
        displayChannelImages = 1;
        if(displayChannelImages == 1)                              
            if(params.subtractWater)
                currentChannelCombinedSignalFid = GERecon('Spectro.SingleVoxel.SubtractWater', combinedSignalFid, combinedReferenceFid);
            else
                currentChannelCombinedSignalFid = combinedSignalFid;
            end

            transformedChannelData = GERecon('Spectro.SingleVoxel.Transform', currentChannelCombinedSignalFid);           
            
            % Compute starting/ending indices to plot based on the
            % parameters for this acquisition/reconstruction (see
            % documentation for ComputeStartingEndingIndicesToPlot below
            % for more details).
            [startingIndex, endingIndex] = ComputeStartingEndingIndicesToPlot(transformedChannelData, params.lowestPpmToPlot, params.ppmRangeToPlot, params.centerFrequencyMHz, params.temperatureAdjustedCenterFrequencyPpm, params.spectralBandwidthHz);

            figure(channelImagesFigureHandle);            
            subplot(2,1,1,'replace');
            totalXAxisPoints = endingIndex-startingIndex;
            plot(real(transformedChannelData(startingIndex:endingIndex)));xlim([0,(totalXAxisPoints-1)]);colorbar;title(['Channel ' num2str(channel)]);
            subplot(2,1,2,'replace');
            plottedSpectrumImage = GERecon('Spectro.SingleVoxel.Plot', transformedChannelData, params.lowestPpmToPlot, params.ppmRangeToPlot);
            imagesc(plottedSpectrumImage);colormap(gray);colorbar;title(['Channel ' num2str(channel) ' Image Plot']);
        end
        
        %% Accumulate Multi-Channel Data
        % Individual channel data is scaled by the maximum magnitude
        % squared of the non-water-suppressed reference FID (if 
        % non-water-suppressed reference FIDs ARE acquired) or by the 
        % maximum magnitude squared of the water-suppressed signal FID 
        % (if non-water-suppressed reference FIDs are NOT acquired).
        % The cumulative channel weights are saved to scale the channel
        % combined FIDs back down after the channel combine.
        accumulatedChannelWeight = accumulatedChannelWeight + weight;
        channelCombinedReference = channelCombinedReference + sqrt(weight) * combinedReferenceFid;
        channelCombinedSignal = channelCombinedSignal + sqrt(weight) * combinedSignalFid;
    end 
    
    %% Water Subtraction
    % Water subtraction uses a non-water-suppressed reference FID to
    % estimate the residual water contribution to the water-suppressed
    % FID's signal. The residual water signal in the water-suppressed
    % signal FID is subtracted. Additional parameters required for water
    % subtraction are taken from the currently active pfile's header.
    % The non-water-suppressed reference FID is plotted along with the
    % water-suppressed signal FID before and after sater subtraction.
        
    % Subtract water if water subtraction is enabled.
    if(params.subtractWater)
        figure;
        subplot(3,1,1,'replace');
        [startingIndex, endingIndex] = ComputeStartingEndingIndicesToPlot(transformedChannelData, params.lowestPpmToPlot, params.ppmRangeToPlot, params.centerFrequencyMHz, params.temperatureAdjustedCenterFrequencyPpm, params.spectralBandwidthHz);
        transformedChannelData = GERecon('Spectro.SingleVoxel.Transform', channelCombinedSignal); 
        plot(real(transformedChannelData(startingIndex:endingIndex)));xlim([0 (endingIndex-startingIndex)]);title(['Channel ' num2str(channel) ' - Channel Combined Signal']);    

        subplot(3,1,2,'replace');
        transformedChannelData = GERecon('Spectro.SingleVoxel.Transform', channelCombinedReference);       
        plot(real(transformedChannelData(startingIndex:endingIndex)));xlim([0 (endingIndex-startingIndex)]);title(['Channel ' num2str(channel) ' - Channel Combined Reference']);                
        
        channelCombinedSignal = GERecon('Spectro.SingleVoxel.SubtractWater', channelCombinedSignal, channelCombinedReference);
        
        subplot(3,1,3,'replace');
        transformedChannelData = GERecon('Spectro.SingleVoxel.Transform', channelCombinedSignal);       
        plot(real(transformedChannelData(startingIndex:endingIndex)));xlim([0 (endingIndex-startingIndex)]);title(['Channel ' num2str(channel) ' - Channel Combined Water Subtracted Signal']);                    
    end   
    
    %% Multi-Channel Scaling
    % Scale down the combined channel data based on the cumulative 
    % weight accumulated in the channel loop above.
    multiChannelWeight = sqrt(accumulatedChannelWeight);
    channelCombinedReference = channelCombinedReference / multiChannelWeight;
    channelCombinedSignal = channelCombinedSignal / multiChannelWeight;
    
    %% Quantitation
    % Quantitate the spectrum using the GE Single Voxel Spectroscopy
    % quantitation algorithm. The algorithm parameters will be taken from
    % the currently active pfile in the GERecon utility. This algorithm
    % will filter and perform a Fourier transform internally. Thus, the
    % input data must be in the time domain.
    %
    % This algorithm computes the following quantitated results and ratios:
    %
    % * Snr
    % * NAA
    % * Cr
    % * Ch
    % * mI
    % * H2O
    % * Noise
    % * NAA/Cr
    % * Cr/Cr
    % * Ch/Cr
    % * mI/Cr
    % * H20/Cr
    if(params.quantitateSpectrum)
        quantitationResults = GERecon('Spectro.SingleVoxel.Quantitate', channelCombinedSignal, channelCombinedReference);
    else
        % These default values are interpreted by the GE Image Viewer to 
        % correctly not display quantitation information if it is not
        % computed, but still display single voxel location information.
        quantitationResults = struct('Snr', -6.0);
        quantitationResults.('NAA') = -6.0;        
        quantitationResults.('Cr') = -6.0;        
        quantitationResults.('Ch') = -6.0;        
        quantitationResults.('mI') = -6.0;
        quantitationResults.('H2O') = -6.0;
        quantitationResults.('Noise') = -6.0;
        quantitationResults.('NAA_Cr') = -6.0;
        quantitationResults.('Cr_Cr') = -6.0;
        quantitationResults.('Ch_Cr') = -6.0;
        quantitationResults.('mI_Cr') = -6.0;
        quantitationResults.('H2O_Cr') = -6.0;
    end
   
    %% Spectral Transform
    % Transform to the frequency domain for plotting. A spectral filter is
    % applied prior to transforming. If needed, the filter can be obtained
    % using the Spectro.SingleVoxel.SpectralFilter command (as is shown
    % below).
    transformedSignal = GERecon('Spectro.SingleVoxel.Transform', channelCombinedSignal);
    filter = GERecon('Spectro.SingleVoxel.SpectralFilter');
    figure;
    subplot(2,1,1);
    plot(abs(filter));title('Filter Applied With Transform Command');
    subplot(2,1,2);
    [startingIndex, endingIndex] = ComputeStartingEndingIndicesToPlot(transformedSignal, params.lowestPpmToPlot, params.ppmRangeToPlot, params.centerFrequencyMHz, params.temperatureAdjustedCenterFrequencyPpm, params.spectralBandwidthHz);
    plot(real(transformedSignal(startingIndex:endingIndex)));title('Final Transformed Spectrum');
    
    %% Plot data to an image and generate DICOM
    % This function will extract a region of interest to plot. The region
    % of interest is based on the user CVs 38 and 39 (rhuser38 and
    % rhuser39) which specify a lowest ppm to plot and a ppm range to
    % plot. (i.e. plot from lowestPpmToPlot to 
    % (lowestPpmToPlot + ppmRangeToPlot)). These values correspond to the 
    % lowestPpmToPlot and ppmRangeToPlot fields in the params structure. 
    % This command plots the real components of the complex input spectrum.    
    plottedSpectrum = GERecon('Spectro.SingleVoxel.Plot', transformedSignal, params.lowestPpmToPlot, params.ppmRangeToPlot);
    figure;imagesc(plottedSpectrum);colormap('gray');colorbar;title('Final Image');
    
    imageNumber = 0;
    matlabDicomPath = fullfile(pfilePath, 'matlabDicom', filesep);
    orientation = GERecon('Pfile.Orientation',1);
    corners = GERecon('Pfile.Corners',1);
    
    % Create Spectroscopy Quantitation Values Dicom Tag
    quantitationValuesDicomTag.Group = hex2dec('0043');
    quantitationValuesDicomTag.Element = hex2dec('1093');
    quantitationValuesDicomTag.VRType = 'DS';   
    quantitationValuesDicomTag.Value = [num2str(quantitationResults.('NAA')) '\' num2str(quantitationResults.('Cr')) '\' num2str(quantitationResults.('Ch')) '\' num2str(quantitationResults.('mI')) '\' num2str(quantitationResults.('H2O'))];

    % Create Spectroscopy Quantitation Ratios Tag
    quantitationRatiosDicomTag.Group = hex2dec('0043');
    quantitationRatiosDicomTag.Element = hex2dec('1094');
    quantitationRatiosDicomTag.VRType = 'DS';   
    quantitationRatiosDicomTag.Value = [num2str(quantitationResults.('NAA_Cr')) '\' num2str(quantitationResults.('Cr_Cr')) '\' num2str(quantitationResults.('Ch_Cr')) '\' num2str(quantitationResults.('mI_Cr')) '\' num2str(quantitationResults.('H2O_Cr'))];    

    % Create Spectroscopy Parameters Tag (contains SNR and Noise
    % Quantitation Results)
    nucleus = header.RawHeader.user6;
    spectroParametersDicomTag.Group = hex2dec('0043');
    spectroParametersDicomTag.Element = hex2dec('108f');
    spectroParametersDicomTag.VRType = 'DS';   
    spectroParametersDicomTag.Value = [num2str(nucleus) '\' num2str(quantitationResults.('Snr')) '\' num2str(quantitationResults.('Noise'))];        

    GERecon('Dicom.Write', [matlabDicomPath 'Image_0000.dcm'], plottedSpectrum, imageNumber, orientation, corners, (header.SeriesData.se_no * 100), header.SeriesData.se_desc, quantitationValuesDicomTag, quantitationRatiosDicomTag, spectroParametersDicomTag);
end

function [processedFid, cumulativePhaseCorrectionVector] = ProcessFid(fidToProcess, processingParameters)
%% ProcessFids
% Process the given FIDs using the included single voxel algorithms based on
% the values in the processingParameters structure.

    processedFid = fidToProcess;
    cumulativePhaseCorrectionVector = ones(size(fidToProcess));
    appliedScalar = 1.0;
    
    %% Normalize
    % Divide each complex value in the FID by the maximum magnitude of the
    % FID. The scaled FID and  maximum magnitude scalar are returned. The 
    % maximum magnitude scalar may be used to scale the data back up at a 
    % later processing step if desired.
    if(processingParameters.normalize)
        [processedFid, appliedScalar] = GERecon('Spectro.SingleVoxel.Normalize', fidToProcess);
    end
    
    %% Low Pass Filter
    % Filter the given FID with a low pass FIR filter. The FID is passed 
    % through the FIR filter once in the forward direction followed 
    % by once in the reverse direction. This results in zero phase 
    % distortion in the filtered FID and is analogous to MATLAB's filtfilt
    % function.
    if(processingParameters.lowPassFilter)
        processedFid = GERecon('Spectro.SingleVoxel.LowPassFilter', processedFid);
    end
    
    %% Shift Max Frequency to DC
    % Applies a linear phase ramp to the input data to shift the frequency
    % with maximum magnitude to the DC location in the frequency domain.  
    % This is done by transforming the input data, finding the maximum 
    % magnitude point in frequency space, and constructing a linear phase 
    % shift correction vector that will bring the maximum magnitude point 
    % in the frequency domain to the DC frequency.  The correction vector 
    % is applied to the input data in the time domain.  The complex 
    % correction vector is also returned.
    if(processingParameters.shiftMaximumMagnitudeFrequencyToDC)
        [processedFid, cumulativePhaseCorrectionVector] = GERecon('Spectro.SingleVoxel.ShiftMaximumMagnitudeFrequencyToDC', processedFid);
    end

    %% Subtract First Point Phase
    % Compute the phase of the first point in the input FID and subtract 
    % that phase angle from all points in the FID. This results in the 
    % first point in the FID having a phase angle of zero. A phase 
    % correction vector with a constant complex value at each point is also 
    % returned. The constant complex value is the phase that was subtracted 
    % from each point in the input FID:
    %  $e^{(-i*phaseOfFirstPoint)}$
    if(processingParameters.subtractPhaseOfFirstPointInFid)
        [processedFid, phaseCorrectionVector] = GERecon('Spectro.SingleVoxel.SubtractPhaseOfFirstPointInFid', processedFid);        
        cumulativePhaseCorrectionVector = cumulativePhaseCorrectionVector .* phaseCorrectionVector;
    end
    
    %% Remove Linear Phase Trend
    % Estimates a linear phase trend by computing and unwrapping the phase
    % of the input FID. A linear phase ramp extending from the phase of the 
    % first point in the FID to the phase of the last point in the FID is 
    % constructed and subtracted from the input FID. The input FID with the 
    % linear phase trend removed and the phase correction vector are both 
    % returned. The correction vector is in complex format representing:
    %  $e^{-i*linearPhaseEstimateAtGivenPoint}$
    if(processingParameters.removeLinearPhaseTrend)
        [processedFid, phaseCorrectionVector] = GERecon('Spectro.SingleVoxel.RemoveLinearPhaseTrend', processedFid);        
        cumulativePhaseCorrectionVector = cumulativePhaseCorrectionVector .* phaseCorrectionVector;
    end

    %% Smoothed Phase Removal
    % Compute and unwrap the phase of the input vector. Smooth the phase
    % that exists across the FID using a spline smoothing function. The 
    % smoothed phase estimate is subtracted from the input vector. The 
    % input vector with the smoothed phase subtracted is returned. The 
    % smoothed phase vector that was subtracted from the input data is also 
    % returned. The phase correction vector is in complex format
    % representing:
    %  $e^{-i*smoothedPhaseEstimateAtGivenPoint}$
    if(processingParameters.smoothedPhaseRemoval)
        [processedFid, phaseCorrectionVector] = GERecon('Spectro.SingleVoxel.SmoothedPhaseRemoval', processedFid);        
        cumulativePhaseCorrectionVector = cumulativePhaseCorrectionVector .* phaseCorrectionVector;
    end    
    
    %% Scale Data Up
    % If the FID was normalized above, scale the data back up by the same
    % amount the data was previously scaled down.
    if(processingParameters.normalize)
        processedFid = processedFid * appliedScalar;
    end    
end

function [startingIndex, endingIndex] = ComputeStartingEndingIndicesToPlot(spectrumVector, lowestPpmToPlot, ppmRangeToPlot, centerFrequencyMHz, centerFrequencyPpmValue, bandwidthHz)
%% ComputeStartingEndingIndicesToPlot
% Given the ppm value of the center frequency, the acquisition bandwidth, 
% and the center frequency in MHz, compute the starting/ending vector 
% indices that contain the ppm range defined by lowestPpmToPlot and 
% ppmRangeToPlot. These ppm inputs map to rhuser38 and rhuser39 in product 
% single voxel spectroscopy sequences.
%
%   Example: To plot from -0.6ppm to 4.3ppm 
%            lowestPpmToPlot = -0.6ppm 
%            ppmRangeToPlot = 4.9ppm 
%          
%   [startingIndex, endingIndex] = ComputeStartingEndingIndicesToPlot(spectrumVector, -0.6, 4.9, centerFrequencyMHz, centerFrequencyPpmValue, bandwidthHz);
%

%% Compute points to plot
%
% $ppmDeltaFromCenterFreq=\frac{(frequency\:\delta\:wrt\:Center\:Freq\:Hz)}{(Center\:Frequency\:MHz)}$
%
% $pointsPerHz=\frac{Num\:Points\:In\:Vector}{Bandwidth\:Hz}$
%
% $(frequency\:\delta\:wrt\:Center\:Freq\:Hz)*pointsPerHz\:=\:(points\:from\:max\:frequency\:to\:point\:of\:interest)$
%
    highestPpmValueToPlot = lowestPpmToPlot + ppmRangeToPlot;
    ppmDeltaFromCenterFrequencyPpmToHighestPpmValueToPlot = centerFrequencyPpmValue - highestPpmValueToPlot;
    frequencyDeltaInHzFromCenterFreqToHighestPpmToPlot = centerFrequencyMHz * ppmDeltaFromCenterFrequencyPpmToHighestPpmValueToPlot;

    ppmDeltaFromCenterFrequencyPpmToLowestPpmValueToPlot = centerFrequencyPpmValue - lowestPpmToPlot;
    frequencyDeltaInHzFromCenterFreqToLowestPpmToPlot = centerFrequencyMHz * ppmDeltaFromCenterFrequencyPpmToLowestPpmValueToPlot;            

    maxMagnitude = max(abs(spectrumVector));
    tolerance = 0.0001;
    waterPeakIndex = find( abs(abs(spectrumVector) - maxMagnitude) < tolerance );                    
    pointsPerHz = size(spectrumVector,1) / bandwidthHz;
    startingIndex = ceil(waterPeakIndex + frequencyDeltaInHzFromCenterFreqToHighestPpmToPlot * pointsPerHz);
    endingIndex = ceil(waterPeakIndex + frequencyDeltaInHzFromCenterFreqToLowestPpmToPlot * pointsPerHz); 
    
    if(startingIndex < 1) || (startingIndex > size(spectrumVector,1))
       startingIndex = 1;
    end

    if(endingIndex > size(spectrumVector,1))
        endingIndex = size(spectrumVector,1);
    end
end
