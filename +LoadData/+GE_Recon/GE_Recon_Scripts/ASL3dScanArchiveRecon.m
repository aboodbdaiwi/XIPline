function ASL3dScanArchiveRecon(ASL3dScanArchiveFilePath)
%% 3DASLScanArchiveRecon - Reconstruct data from a 3DASL scan archive
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

% Multi-echo is not considered in this application.
    
    %% Load the given Spiral ScanArchive 
    % All secondary input files that are in this ScanArchive are extracted 
    % to the directory that this scan archive exists in 
    archive = GERecon('Archive.Load', ASL3dScanArchiveFilePath);
    archiveFilePath = fileparts(ASL3dScanArchiveFilePath);

    %% Parameter and variable initialization 
    % Extract information from the archive's DownloadData object to set 
    % parameters used in this reconstruction.
    numPointsPerArm = archive.DownloadData.rdb_hdr_rec.rdb_hdr_da_xres;
    numArms = archive.DownloadData.rdb_hdr_rec.rdb_hdr_da_yres-1;
    imageSize = archive.DownloadData.rdb_hdr_rec.rdb_hdr_im_size;
    numChannels = archive.DownloadData.rdb_hdr_rec.rdb_hdr_dab(1).stop_rcv ...
        - archive.DownloadData.rdb_hdr_rec.rdb_hdr_dab(1).start_rcv + 1;
    numSlices = archive.Slices;
    kissoffs = archive.DownloadData.rdb_hdr_rec.rdb_hdr_slblank;
    numProcessedSlices = numSlices-2*kissoffs;
    zStart = kissoffs + 1;
    zEnd = numSlices - kissoffs;
    isChopDataInZ = ~bitget(archive.DownloadData.rdb_hdr_rec.rdb_hdr_data_format, 3);
    perf_scaling_factor = archive.DownloadData.rdb_hdr_rec.rdb_hdr_asl_perf_weighted_scale;

    rawKSpace = zeros(numPointsPerArm, numArms, numChannels, numSlices);
    processedKSpace = zeros(numPointsPerArm, numArms, numChannels, numProcessedSlices);
    channelImages = zeros(imageSize, imageSize, numChannels);
    armStackData = zeros(numPointsPerArm, numChannels, numSlices);

    %% Loop and Recon
    % Loop over all control packets in this scan archive. Each control
    % packet corresponds to either a Spiral2D Readout (opcode 1) or an
    % end of pass/end of scan(opcode 0). PassDone and ScanDone control
    % packets do not have data associated with them; thus, only data
    % associated with opcode 1 (ProgrammableOpcode) is reconstructed.
    passCounter = 1;
    acquisitionCounter = 1;
    nexCounter = 0;
    sliceCounter = 0;
    armCounter = 0;
    for dabPacketIndex = 1:archive.ControlCount
        currentControl = GERecon('Archive.Next', archive);

        if(currentControl.opcode == 1)
            armIndex = currentControl.viewNum;
            
            oneBasedPassSliceNum = currentControl.sliceNum + 1;
            sliceInfo = GERecon('Archive.Info', archive, passCounter, oneBasedPassSliceNum);
            geometricSliceNumber = sliceInfo.Number;
            sliceInfoTable(geometricSliceNumber) = sliceInfo;
            
            % Combine a Frame into an esisting Sources.
            % Operation = 0: Result = Frame
            %             1: Result = Frame + Source
            %             2: Result = Frame - Source
            %             3: Result = Source - Frame
            if(0 == currentControl.operation)
                rawKSpace(:, armIndex, :, geometricSliceNumber) = double(currentControl.Data);
                % Counts slice number
                if(1 == armIndex)
                    sliceCounter = sliceCounter + 1;
                end
                % Counts arm number
                if(1 == geometricSliceNumber)
                    armCounter = armCounter + 1;
                end
            elseif(1 == currentControl.operation)
                rawKSpace(:, armIndex, :, geometricSliceNumber) = ...
                    double(currentControl.Data) + ...
                    squeeze(rawKSpace(:, armIndex, :, geometricSliceNumber));
            elseif(2 == currentControl.operation)
                rawKSpace(:, armIndex, :, geometricSliceNumber) = ...
                    double(currentControl.Data) - ...
                    squeeze(rawKSpace(:, armIndex, :, geometricSliceNumber));
            elseif(3 == currentControl.operation)
                rawKSpace(:, armIndex, :, geometricSliceNumber) = ...
                    squeeze(rawKSpace(:, armIndex, :, geometricSliceNumber))...
                    - double(currentControl.Data);
            end
            
            % Counts nex number
            if((1 == armIndex) && (1 == geometricSliceNumber) && (1 == acquisitionCounter))
                nexCounter = nexCounter + 1;
            end
        elseif(currentControl.opcode == 0)
            passCounter = passCounter + 1;
            if((numSlices == sliceCounter) && (numArms == armCounter))
                rawKSpace = rawKSpace / nexCounter;
                if( acquisitionCounter == 1 )
                    rawKSpace = rawKSpace * perf_scaling_factor;
                end
                for armIndex = 1:numArms
                    armStackData(:,:,:) = squeeze(rawKSpace(:, armIndex, :, :));
                    if(isChopDataInZ)
                        armStackData(:, :, 2:2:end) = -armStackData(:, :, 2:2:end);
                    end
                    
                    transformedData = GERecon('ZTransform', armStackData);
                    processedKSpace(:, armIndex, :, :) = transformedData(:, :, zStart:zEnd);
                end

                for sliceIndex = 1:numProcessedSlices
                    % Convert spiral K-Space to gridded K-space
                    gridKSpace = GERecon('Spiral.Regrid', squeeze(processedKSpace(:,:,:,sliceIndex)));

                    for channelIndex = 1:numChannels
                        % Transform K-Space
                        channelImages(:,:,channelIndex) = ifft2(gridKSpace(:,:,channelIndex));
                    end
                    
                    % Get corners and orientation for this slice location
                    geometricSliceNumber = sliceInfoTable(sliceIndex).Number;
                    corners = sliceInfoTable(sliceIndex).Corners;
                    orientation = sliceInfoTable(sliceIndex).Orientation;

                    % Apply Complex Channel Combination
                    combinedImage = GERecon('PhaseCorrectAndCombine', channelImages);

                    % Create Magnitude Image
                    magnitudeImage = abs(combinedImage);

                    % Apply Gradwarp
                    gradwarpedImage = GERecon('Gradwarp', magnitudeImage, corners);

                    % Clip the image
                    clippedImage = gradwarpedImage;
                    clippedImage(clippedImage < 0) = 0;

                    % Orient the image
                    finalImage = GERecon('Orient', clippedImage, orientation);

                    maskedFinalImage = GERecon('Spiral.ApplyCircularMask', finalImage);

                    % Display
                    imagesc(maskedFinalImage);colormap(gray);colorbar;axis off;title(['Slice: ' num2str(sliceIndex)]);
                    title(['Slice: ' num2str(sliceIndex)]);

                    % Save DICOMs
                    matlabDicomPath = fullfile(archiveFilePath, 'matlabDicoms', filesep);
                    filename = [matlabDicomPath 'Image_' num2str(geometricSliceNumber+numProcessedSlices*(acquisitionCounter-1)-1, '%05d') '.dcm'];
                    GERecon('Dicom.Write', filename, maskedFinalImage, geometricSliceNumber, orientation, corners);
                end % end of sliceIndex
                
                acquisitionCounter = acquisitionCounter + 1;
                passCounter = 1;
                sliceCounter = 0;
                armCounter = 0;
            end % 3D ASL scan consists of one PD and one PW
        end % opcode
    end % loop of all control packets
    
    GERecon('Archive.Close', archive);
end
