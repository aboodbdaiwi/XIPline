function [ phaseCorrectionHandle, nextDabPacketIndex ] = EpiReferenceScanRecon( archiveHandle, rowFlipHandle, plotAppliedCoefficients )
%% EpiReferenceScanRecon - Compute EPI Phase Correction coefficients from ScanArchive
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
%   [ phaseCorrectionHandle, nextDabPacketIndex ] = EpiReferenceScanRecon(archiveHandle, rowFlipHandle, plotAppliedCoefficients)
%   Computes EPI phase correction coefficients from the given scan archive.
%   A phase correction handle is returned which contains the computed
%   coefficients. The next dab packet index in the scan archive is also
%   returned. This is useful for archives which contain an integrated
%   reference acquisition. For these archives, reconstruction can continue
%   from the dab packet index that the reference processing stopped at.
%   An archiveHandle containing the data to compute coefficients from is 
%   required. An archiveHandle is required rather than a ScanArchive path
%   to support the case of integrated reference scan in which the archive
%   is opened by a parent Epi recon script which calls this function to
%   compute coefficients from a portion of the data in the archive.
%   A row flip handle is required to apply the EPI row flip prior to
%   computing coefficients from the reference data. Similar to the
%   archiveHandle, a handle is required here because this allows a parent
%   script to open the row flip handle and share it with this script.
%   An optional plotAppliedCoefficients parameter may be supplied to 
%   control whether the computed coefficients are applied to the original 
%   row-flipped reference data and plotted. Applying the coefficients to 
%   the reference data is not required for recon, but is useful to 
%   demonstrate the effect that the coefficients have.


    %% Load
    % Initialize an empty phase correction
    % reference handle. The phase correction reference handle will be filled
    % with coefficients and returned to the caller. The filled handle can be
    % used to apply phase correction using the computed coefficients.
    phaseCorrectionHandle = GERecon('Epi.LoadReference');
    
    %% Account for reference views
    % Some multi-phase EPI scans acquire additional non-phase-encoded views
    % at the top or bottom of the EPI echo train. These views are not
    % needed for reference scan processing. This code determines the number
    % of additional reference views acquired with the reference scan data
    % so they can be discarded when the data is processed.    
    refViewsTop = archiveHandle.DownloadData.rdb_hdr_rec.rdb_hdr_extra_frames_top;
    refViewsBottom = archiveHandle.DownloadData.rdb_hdr_rec.rdb_hdr_extra_frames_bot;    
    
    %% Initialize figures
    % Initialize figures used to display coefficients
    coefficientsFigureHandle = figure;        
    if(nargin < 3)
        plotAppliedCoefficients = 1;
    end    
    if(plotAppliedCoefficients)
        applyCoefficientsFigureHandle = figure;
    end

    %% Loop and compute coefficients
    % Compute linear and constant phase correction coefficients from
    % reference scan data. This loop will continue until one phase of the 
    % scan has been reconstructed. The reference scan is always a single
    % phase of the scan.
    % Counters are initialized here to keep track of when one phase of the
    % scan has completed. Each packet will either be a scan control 
    % packet (opcode = 0), hyper dab packet (opcode = 6), or a diffusion
    % hyper dab packet (opcode = 14). Only the hyper dab or diffusion hyper
    % dab packets have data associated with them.
    acquisitionPassCounter = 0;
    numAcquisitionsPerRepetition = archiveHandle.DownloadData.rdb_hdr_rec.rdb_hdr_npasses / archiveHandle.DownloadData.rdb_hdr_rec.rdb_hdr_reps;
    for dabPacketIndex = 1:archiveHandle.ControlCount          
        % Obtain the next control packet from the scan archive and
        % determine if this packet has associated acquired data to compute
        % coefficients from
        currentControl = GERecon('Archive.Next', archiveHandle);        
        packetHasAcquisitionData = 0;
        if(currentControl.opcode == 0)
            packetHasAcquisitionData = 0;
            acquisitionPassCounter = acquisitionPassCounter + 1;
            if(acquisitionPassCounter == numAcquisitionsPerRepetition)
                % One repetition of the scan has been processed. The 
                % integrated reference scan data is always in the first 
                % phase of the scan. Thus, now that the first phase of the
                % scan is complete, the ref scan processing is complete.
                % Set the output nextDabPacketIndex to the next dab packet
                % index to process. The calling script can pickup the
                % reconstruction from this point.
                nextDabPacketIndex = dabPacketIndex + 1;
                return;
            end
        elseif(currentControl.opcode == 6)
            packetHasAcquisitionData = 1;
        elseif(currentControl.opcode == 14)
            packetHasAcquisitionData = 1;
        end
       
        if(packetHasAcquisitionData) 
            % Obtain slice information for the slice we're about to reconstruct
            acquisitionPassCounterOneBased = acquisitionPassCounter + 1;
            oneBasedPassSliceNum = currentControl.sliceNum + 1;
            sliceInfo = GERecon('Archive.Info', archiveHandle, acquisitionPassCounterOneBased, oneBasedPassSliceNum);
            geometricSliceNumber = sliceInfo.Number;
            
            % Extract acquired data and handle bottom-up kSpace acquisitions
            kSpace = currentControl.Data;
                        
            if(currentControl.viewSkip < 0)
                % This is a bottom up EPI scan. Thus, the kSpace must be
                % flipped along the phase encoding direction prior to
                % processing. This is because the raw data is held in
                % acquisition time order. Since the views are acquired
                % bottom-up. The first view in the raw data matrix should
                % actually be the last/bottom view of the kSpace matrix
                kSpace = flip(kSpace,3);
            end            
            
            % Account for additional reference views that may be acquired with multi-phase scans
            kSpaceTotalNumViews = size(kSpace,3);                
            kSpaceWithoutRefViews = kSpace(:,:,(refViewsTop+1):(kSpaceTotalNumViews-refViewsBottom));
            
            numChannels = size(kSpaceWithoutRefViews,2);            
            for channel = 1:numChannels
                % Apply rowflip
                rowFlippedKSpace = GERecon('Epi.ApplyImageRowFlip', rowFlipHandle, squeeze(kSpaceWithoutRefViews(:,channel,:))); 
                
                % Compute phase correction coefficients, store in reference handle, and plot coefficients                
                computedCoefficients = GERecon('Epi.ComputeCoefficients', phaseCorrectionHandle, rowFlippedKSpace, geometricSliceNumber, channel);
                figure(coefficientsFigureHandle);
                subplot(2,1,1,'replace');
                plot(computedCoefficients(:,1));title(['Constant Coefficients - Slice: ' num2str(geometricSliceNumber) ', Channel: ' num2str(channel)]);
                subplot(2,1,2,'replace');
                plot(computedCoefficients(:,2));title(['Linear Coefficients - Slice: ' num2str(geometricSliceNumber) ', Channel: ' num2str(channel)]);            
                drawnow;
            end
            
            % If Phase alignment bit is set in the pool header
            % phase alignment is performed on the computed coefficients
            GERecon('Epi.PhaseAlignment',phaseCorrectionHandle, geometricSliceNumber);
            
            for channel = 1:numChannels
                if(plotAppliedCoefficients)
                    % Apply phase correction
                    phaseCorrectedKSpace = GERecon('Epi.ApplyPhaseCorrection', phaseCorrectionHandle, rowFlippedKSpace, geometricSliceNumber, channel);

                    % Plot data before and after shifiting data
                    figure(applyCoefficientsFigureHandle);
                    subplot(2,1,1,'replace');
                    imagesc(abs(rowFlippedKSpace)');colormap(gray);colorbar;axis off;title(['Raw Slice ' num2str(geometricSliceNumber) ', Channel ' num2str(channel)]);
                    subplot(2,1,2,'replace');
                    imagesc(abs(phaseCorrectedKSpace)');colormap(gray);colorbar;axis off;title(['Corrected Slice: ' num2str(geometricSliceNumber) ', Channel: ' num2str(channel)]);
                    drawnow;  
                end
            end
        end
    end            
end

