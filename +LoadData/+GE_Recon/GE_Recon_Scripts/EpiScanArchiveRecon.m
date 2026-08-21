function [ ] = EpiScanArchiveRecon( epiScanArchiveFilePath )
%% EpiScanArchiveRecon - Reconstruct data from an EPI scan archive
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
%   EpiScanArchiveRecon(epiScanArchiveFilePath)
%   will reconstruct the single-shot EPI data in the given scan archive. 
%   All secondary input files are obtained from the files contained in the 
%   scan archive itself. These secondary input files are extracted from the
%   scan archive when it is loaded.
%   This script currently supports a basic single-shot EPI reconstruction 
%   with the following steps/algorithms:
%       1. Extract Raw Data
%       2. Row flip based on rowflip.param file in scan archive
%       3. Static EPI Phase Correction (based on either the ref.h5 file in the
%          scan archive or the coefficients computed from an integrated
%          reference pass of a diffusion acquisition)
%       4. Ramp sampling based on vrgf.dat in scan archive
%       5. 2D Transform
%       6. Channel Combination (Sum of Squares or ASSET)
%       7. Gradwarp
%       8. DICOM Generation
%
%   Support for other EPI specific reconstruction steps (such as
%   magnitude or complex averaging, realtime field adjustment, etc...) may
%   be added by following the examples in the existing Epi pfile based 
%   rehearsal recon scripts

    %% Load the given Epi ScanArchive
    % Load the given scan archive and extract all secondary input files
    [archiveFilePath archiveFileName archiveFileExt] = fileparts(epiScanArchiveFilePath);
    archive = GERecon('Archive.Load', epiScanArchiveFilePath);    
    
    %% Load rowflip.param
    % rowflip.param is used to control which EPI readouts are flipped along
    % the readout axis. For a single shot scan, this would be every other
    % readout.
    rowFlipHandle = GERecon('Epi.LoadRowFlip', 'rowflip.param');    
    
    %% Load reference data
    % Based on the download data in the archive, determine if this scan has
    % an integrated reference scan or if the reference coefficients should
    % be loaded from the ref.h5 file on disk. 
    isIntegratedStaticPC = archive.DownloadData.rdb_hdr_rec.rdb_hdr_ref == 5;
    isIntegratedRefBitSet = bitget(archive.DownloadData.rdb_hdr_rec.rdb_hdr_pcctrl, 5);   
    isIntegratedRef = isIntegratedStaticPC || isIntegratedRefBitSet;    
    firstImageAcquisitionPhase = 0;
    firstImageDabPacketIndex = 1;
    if(isIntegratedRef)
        [phaseCorrectionHandle, nextDabPacketIndex] = EpiReferenceScanRecon(archive, rowFlipHandle);
        firstImageDabPacketIndex = nextDabPacketIndex;
        firstImageAcquisitionPhase = 1;
    else       
        % Read phase correction coefficients from ref.h5
        phaseCorrectionHandle = GERecon('Epi.LoadReference', 'ref.h5');
    end   
    
    %% Load vrgf.dat
    % Look in the scan archive directory for a vrgf.dat file. Certain fMRI
    % scans may use an alternative vrgf_kernels.dat file. Support may be
    % added for vrgf_kernels.dat by following the example in
    % EpiMultiPhaseRecon.m
    vrgfHandle = GERecon('RampSampling.Load', 'vrgf.dat');    
    kernelMatrix = GERecon('RampSampling.KernelMatrix', vrgfHandle);
    figure;
    subplot(2,1,1);
    imagesc(kernelMatrix);colormap(gray);colorbar;axis off;title('VRGF Interpolation Kernel Matrix');
    subplot(2,1,2);
    numRows = size(kernelMatrix, 1);
    plot(kernelMatrix(numRows/4, :));title(['Sinc Interpolation Kernel for Interpolated Point: ' num2str(numRows/4) ' of ' num2str(numRows)]);    

    %% Load ASSET Calibration (if it exists)
    % Check the exam data path for any ASSET calibration files that were
    % included in the archive. If an ASSET calibration file is included in
    % the archive then load the calibration here.
    examDataDirectory = num2str(archive.DownloadData.rdb_hdr_exam.ex_no);
    useAsset = 0;
    if(exist(num2str(archive.DownloadData.rdb_hdr_exam.ex_no), 'dir'))
       % If the exam data directory exists, check if it contains an ASSET calibration file 
       fileList = dir(examDataDirectory);
       for i=1:size(fileList,1)
          findResult = strfind(fileList(i).name, 'Asset');
          if(size(findResult) > 0)
             assetCalibrationFile = fullfile(num2str(archive.DownloadData.rdb_hdr_exam.ex_no), fileList(i).name);
             disp(['Loading Calibration File: ' assetCalibrationFile]); 
             useAsset = 1;
             GERecon('Asset.LoadCalibration', assetCalibrationFile);    
          end
       end
    end
    
    %% Parameter and variable initialization
    % Extract information from the archive's DownloadData object to set
    % parameters used in this reconstruction.
    % The refViewsTop and refViewsBottom variables are used to discard
    % readouts acquired for dynamic phase correction. The dynamic phase
    % correction algorithm is not currently supported by this script.
    refViewsTop = archive.DownloadData.rdb_hdr_rec.rdb_hdr_extra_frames_top;
    refViewsBottom = archive.DownloadData.rdb_hdr_rec.rdb_hdr_extra_frames_bot;     
    imageSize = archive.DownloadData.rdb_hdr_rec.rdb_hdr_im_size;    
    % Determine if the phase encoding direction was reversed for this scan.
    % Reverse phase encoding direction is supported for DW-EPI scans or
    % other single-echo EPI scans. When reverse phase encoding is used, the
    % image must be reversed along the phase-encoding direction after the
    % 2D transform.
    isDiffusion = bitget(archive.DownloadData.rdb_hdr_rec.rdb_hdr_data_collect_type1, 22);
    numEchoes = archive.DownloadData.rdb_hdr_rec.rdb_hdr_nechoes;
    flipBitsSet = bitget(archive.DownloadData.rdb_hdr_rec.rdb_hdr_dacq_ctrl,2) || bitget(archive.DownloadData.rdb_hdr_rec.rdb_hdr_dacq_ctrl,3);
    alternatePhaseEncodingDirection = (flipBitsSet && isDiffusion) || (flipBitsSet && numEchoes == 1);    
    
    %% Loop and Recon
    % Loop over all control packets in this scan archive. Each control
    % packet corresponds to either an Epi Readout (opcode 6 or 14) or an
    % end of pass/end of scan(opcode 0). PassDone and ScanDone control
    % packets do not have data associated with them; thus, only data
    % associated with opcode 6 (HyperScanOpcode) or 14 (DiffusionHyperScanOpcode)
    % is reconstructed.
    % If this scan contained an integrated reference pass, the
    % EpiReferenceScanRecon script already processed some control packets
    % and returned the next control packet to process. Thus, the loop
    % starts from the [firstImageDabPacketIndex] variable and continues until
    % all control packets are processed.
    numAcquisitionsPerRepetition = archive.DownloadData.rdb_hdr_rec.rdb_hdr_npasses / archive.DownloadData.rdb_hdr_rec.rdb_hdr_reps;
    acquisitionPassCounter = 0;
    phaseIndex = firstImageAcquisitionPhase;
    imageNumber = 0;
    figure;    
    for dabPacketIndex = firstImageDabPacketIndex:archive.ControlCount        
        % Extract the next control packet and associated data from the scan
        % archive. Inspect the opcode in the control packet to determine if
        % the packet contains data to reconstruct.
        currentControl = GERecon('Archive.Next', archive);        
        reconstructDataForThisPacket = 0;
        if(currentControl.opcode == 0)
            reconstructDataForThisPacket = 0;
            % For diffusion scans, the same volume may be acquired multiple times for the various diffusion passes (T2, B-Value/Diffusion Direction). 
            % Thus, when this counter reaches the number of acquisitions per repetition, reset the counter back to 0.
            acquisitionPassCounter = acquisitionPassCounter + 1;
            if(acquisitionPassCounter == numAcquisitionsPerRepetition)
                acquisitionPassCounter = 0;
                phaseIndex = phaseIndex + 1;
            end
        elseif(currentControl.opcode == 6)
            reconstructDataForThisPacket = 1;
        elseif(currentControl.opcode == 14)
            reconstructDataForThisPacket = 1;
        end
       
        if(reconstructDataForThisPacket)       
            % Obtain slice information for the slice we're about to reconstruct
            acquisitionPassCounterOneBased = acquisitionPassCounter + 1;
            oneBasedPassSliceNum = currentControl.sliceNum + 1;
            sliceInfo = GERecon('Archive.Info', archive, acquisitionPassCounterOneBased, oneBasedPassSliceNum);
            geometricSliceNumber = sliceInfo.Number;
            
            % Extract raw image kSpace. Additional reference views used
            % with dynamic phase correction are not supported in this
            % rehearsal recon; thus, they are discarded here. Also, handle 
            % the bottom up case in which view-increment is negative prior 
            % to any other processing.
            kSpace = currentControl.Data;
            
            if(currentControl.viewSkip < 0)
                % This is a bottom up EPI scan. Thus, the kSpace must be
                % flipped along the phase encoding direction prior to
                % processing. This is because the raw data is held in
                % acquisition time order. Since the views are acquired
                % bottom-up. The first view in the raw data matrix should
                % actually be the last view of the kSpace matrix.
                kSpace = flip(kSpace,3);
            end
            
            kSpaceTotalNumViews = size(kSpace,3);                
            kSpaceWithoutRefViews = kSpace(:,:,(refViewsTop+1):(kSpaceTotalNumViews-refViewsBottom));
            
            % Reconstruct each channel
            numChannels = size(kSpace,2);
            channelImages = single(zeros(imageSize, imageSize, numChannels));
            for channel=1:numChannels
                % Apply rowflip
                rowFlippedKSpace = GERecon('Epi.ApplyImageRowFlip', rowFlipHandle, squeeze(kSpaceWithoutRefViews(:,channel,:)));                                
                
                % Apply phase correction
                phaseCorrectedKSpace = GERecon('Epi.ApplyPhaseCorrection', phaseCorrectionHandle, rowFlippedKSpace, geometricSliceNumber, channel);

                % Interpolate ramp sampled data
                interpolatedData = GERecon('RampSampling.Interpolate', vrgfHandle, phaseCorrectedKSpace);

                % Filter, 2D-IFFT
                image = GERecon('Transform', interpolatedData);

                % If ASSET and Homodyne are both enabled for this scan then 
                % ASSET is run on both the high pass filtered and low pass 
                % filtered images generated by the Homodyne algorithm. To 
                % enable this use case, the Transform command returns the high 
                % pass filtered and low pass filtered images in indices one and 
                % two of the third dimension of the return image. For this 
                % case, the channel images array must have space to store the
                % additional high pass filtered / low pass filtered images. 
                % Resize the channelImages matrix here to enable this use case.
                % Note that this resize happens only once, the first time 
                % through this loop.
                if( (size(image,3) > 1) && (size(channelImages,4) == 1) )
                    channelImages = single(zeros(size(channelImages,1), size(channelImages,2), size(channelImages,3), 2));                    
                end

                channelImages(:,:,channel, :) = image;
            end
            
            % As described above, if a reverse phase encoding direction was
            % used for this scan, the data must be flipped along the phase
            % encoding direction after the 2D transform. This processing is
            % done here.
            if(alternatePhaseEncodingDirection)
                channelImages = fliplr(channelImages);
            end
            
            % Channel combine processing - either Sum of Squares or ASSET
            if(useAsset)
                channelCombinedImage = GERecon('Asset.Unalias', channelImages, sliceInfo.Corners);
            else                
                channelCombinedImage = GERecon('SumOfSquares', channelImages);
            end            

            % Zero out kissoff views and apply gradwarp
            kissoffViews = archive.DownloadData.rdb_hdr_rec.rdb_hdr_kissoff_views;
            channelCombinedImage(:,1:kissoffViews) = 0;
            channelCombinedImage(:,(end-kissoffViews+1):end) = 0;
            gradwarpImage = GERecon('Gradwarp', abs(channelCombinedImage), sliceInfo.Corners, 'XRMW');

            if(archive.DownloadData.rdb_hdr_rec.rdb_hdr_hnover > 0)             
                reconXRes = archive.DownloadData.rdb_hdr_rec.rdb_hdr_rc_xres;
                reconYRes = archive.DownloadData.rdb_hdr_rec.rdb_hdr_rc_yres;
                gradwarpImage = gradwarpImage * (256 / (reconXRes * reconYRes));
            end
            
            % Rotate/Transpose
            rotatedTransposedSlice = GERecon('Orient', gradwarpImage, sliceInfo.Orientation);            

            % Clip to range of shorts (match product functionality)
            rotatedTransposedSlice(rotatedTransposedSlice < 0) = 0;
            rotatedTransposedSlice(rotatedTransposedSlice > 32767) = 32767;
            finalImage = int16(rotatedTransposedSlice);

            % Display final image
            imagesc(finalImage);colormap(gray);colorbar;axis off;title(['Slice: ' num2str(geometricSliceNumber) 'Phase: ' num2str(phaseIndex)]);
            drawnow;

            % Compute image number and create dicom image
            matlabDicomPath = fullfile(archiveFilePath, 'matlabDicoms', filesep);
            seriesNumber = archive.DownloadData.rdb_hdr_series.se_no * 100;
            GERecon('Dicom.Write', [matlabDicomPath 'Image_' num2str(imageNumber, '%03d') '.dcm'], finalImage, imageNumber, sliceInfo.Orientation, sliceInfo.Corners, seriesNumber);
            imageNumber = imageNumber + 1;            
        end
    end

    GERecon('Archive.Close', archive);
end

