function Spiral2DScanArchiveRecon(spiral2DScanArchiveFilePath)
%% Spiral2DScanArchiveRecon - Reconstruct data from a spiral2D scan archive
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
    
    %% Load the given Spiral ScanArchive 
    % All secondary input files that are in this ScanArchive are extracted 
    % to the directory that this scan archive exists in 
    archive = GERecon('Archive.Load', spiral2DScanArchiveFilePath);  

    %% Parameter and variable initialization 
    % Extract information from the archive's DownloadData object to set 
    % parameters used in this reconstruction.
    numPointsPerArm = archive.DownloadData.rdb_hdr_rec.rdb_hdr_da_xres;
    numArms = archive.DownloadData.rdb_hdr_rec.rdb_hdr_da_yres-1;
    imageSize = archive.DownloadData.rdb_hdr_rec.rdb_hdr_im_size;
    numChannels = archive.DownloadData.rdb_hdr_rec.rdb_hdr_dab(1).stop_rcv ...
        - archive.DownloadData.rdb_hdr_rec.rdb_hdr_dab(1).start_rcv + 1;
    numSlices = archive.Slices;
    kSpace = zeros(numPointsPerArm, numArms, numChannels, numSlices);
    channelImages = zeros(imageSize, imageSize, numChannels);

    %% Loop and Recon
    % Loop over all control packets in this scan archive. Each control
    % packet corresponds to either an Spiral2D Readout (opcode 1) or an
    % end of pass/end of scan(opcode 0). PassDone and ScanDone control
    % packets do not have data associated with them; thus, only data
    % associated with opcode 1 (ProgrammableOpcode) is reconstructed.
    passCounter = 1;
    nexCounter = 0;
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
                kSpace(:, armIndex, :, geometricSliceNumber) = double(currentControl.Data);
            elseif(1 == currentControl.operation)
                kSpace(:, armIndex, :, geometricSliceNumber) = ...
                    double(currentControl.Data) + kSpace(:, armIndex, :, geometricSliceNumber);
            elseif(2 == currentControl.operation)
                kSpace(:, armIndex, :, geometricSliceNumber) = ...
                    double(currentControl.Data) - kSpace(:, armIndex, :, geometricSliceNumber);
            elseif(3 == currentControl.operation)
                kSpace(:, armIndex, :, geometricSliceNumber) = ...
                    kSpace(:, armIndex, :, geometricSliceNumber) - double(currentControl.Data);
            end
            
            % Counts nex number
            if((1 == armIndex) && (1 == geometricSliceNumber))
                nexCounter = nexCounter + 1;
            end
        elseif(currentControl.opcode == 0)
            passCounter = passCounter + 1;
        end
    end

    kSpace = kSpace / nexCounter;
    for sliceIndex = 1:numSlices
        % Convert spiral K-Space to gridded K-space
        gridKSpace = GERecon('Spiral.Regrid', squeeze(kSpace(:,:,:,sliceIndex)));
        
        for channelIndex = 1:numChannels
            % Transform K-Space
            channelImages(:,:,channelIndex) = ifft2(gridKSpace(:,:,channelIndex));
        end
        
        % Get corners and orientation for this slice location
        corners = sliceInfoTable(sliceIndex).Corners;
        orientation = sliceInfoTable(sliceIndex).Orientation;
        
        % Apply Channel Combination
        combinedImage = GERecon('SumOfSquares', channelImages);
        
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
        imagesc(maskedFinalImage);
        title(['Slice: ' num2str(sliceIndex)]);
        
        % Save DICOMs
        filename = ['DICOMs/ScanArchive/image' num2str(sliceIndex) '.dcm'];
        GERecon('Dicom.Write', filename, maskedFinalImage, sliceIndex, orientation, corners);
        
        pause(0.1);
    end
    
    GERecon('Archive.Close', archive);
end
    
