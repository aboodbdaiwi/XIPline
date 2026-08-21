function FlexRecon(pfilePath, varargin)
%% FlexRecon - Reconstruct a 3D Flex Pfile with ARC Parallel Imaging
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
% FlexRecon(pfilePath, varargin)
% will reconstruct the Flex data in the given pfile. The sampling pattern
% file (kacqXXX.txt) can be specified as the second argument to describe 
% the acquisition of the K-Space.
%
% Limitations: intensity correction

    % Load Pfile
    pfile = GERecon('Pfile.Load', pfilePath);
    header = GERecon('Pfile.Header', pfile);

    if pfile.isArcEnabled
        if( LoadKacqWithInputString(varargin) )
            % kacq file specified via varargin loaded
        elseif( LoadKacqBasedOnPfileLocation(pfilePath) )
            % kacq file found based on raw data location and header info
        else
            % No kacq specified or found
            fprintf('No kacq file provided or found, ARC processing will rely on raw data to infer sampling pattern\n');
        end
    end

    % Load KSpace. Since 3D Arc Pfiles contain space for the zipped
    % slices (even though the data is irrelevant), only pull out
    % the true acquired K-Space. Z-transform will zip the slices
    % out to the expected extent.
    acquiredSlices = pfile.slicesPerPass;
    
    % 3D Scaling Factor
    scaleFactor = pfile.scaleFactor3d;
    
    pass = 1; % Single pass example
    
    kSpace = zeros(pfile.xRes, pfile.yRes, pfile.reconstructedSlicesPerPass, pfile.echoes, pfile.channels, 'single');
    
    for echo = 1:pfile.echoes
            kSpaceSingleEcho = zeros(pfile.xRes, pfile.yRes, acquiredSlices, pfile.channels, 'single');
            for slice = 1:acquiredSlices
                    sliceInfo.pass = pass;
                    sliceInfo.sliceInPass = slice;

                    % Load K-Space
                    for channel = 1:pfile.channels
                        if (pfile.isArcEnabled || pfile.isZEncoded)
                            kSpaceSingleEcho(:,:,slice,channel) = GERecon('Pfile.KSpace', sliceInfo, echo, channel);
                        else
                            kSpace(:,:,slice,echo,channel) = GERecon('Pfile.KSpace', sliceInfo, echo, channel);
                        end
                    end
            end
            if pfile.isArcEnabled
                % Synthesize KSpace (Apply ARC)
                kSpaceSingleEcho = GERecon('Arc.Synthesize', kSpaceSingleEcho);
            end

            if (pfile.isArcEnabled || pfile.isZEncoded)
                % Transform Across Slices
                kSpace(:,:,:,echo,:) = ifft(kSpaceSingleEcho, pfile.reconstructedSlicesPerPass, 3);
            end
    end
    clear kSpaceSingleEcho;

    if (pfile.isArcEnabled || pfile.isZEncoded)
        % Scale
        kSpace = kSpace * scaleFactor;
    end
    
    % Loop for each slice/channel to create complex images for 2-pt Dixon
    % processing, and also generate in-phase and out-of-phase DICOM images
    for slice = 1:pfile.reconstructedSlicesPerPass
        for echo = 1:pfile.echoes
            for channel = 1:pfile.channels

                % Transform K-Space
                channelImage = GERecon('Transform', kSpace(:,:,slice,echo,channel));

                channelImages(:,:,channel) = channelImage;
            end

            echoChannelImages(:,:,echo,:) = channelImages;

            % Get slice information (corners and orientation) for this slice location
            sliceInfo.pass = pass;
            sliceInfo.sliceInPass = slice;

            info = GERecon('Pfile.Info', sliceInfo);

            % Apply Channel Combination
            combinedImage = GERecon('SumOfSquares', channelImages);

            % Create Magnitude Image
            magnitudeImage = abs(combinedImage);

            % Apply Gradwarp
            gradwarpedImage = GERecon('Gradwarp', magnitudeImage, info.Corners);

            % Orient the image
            finalImage = GERecon('Orient', gradwarpedImage, info.Orientation);

            % Display
            imagesc(finalImage);

            % Display
            title(['Pass: ' num2str(pass) ' Slice: ' num2str(slice) ' Echo: ' num2str(echo)]);

            % Save DICOMs
            seriesNumber = header.SeriesData.se_no * 100 + echo;
            seriesDescription = ['Echo' num2str(echo) ': ' header.SeriesData.se_desc];
            filename = ['DICOMs' num2str(seriesNumber) '/image' num2str(info.Number) '.dcm'];
            GERecon('Dicom.Write', filename, finalImage, info.Number, info.Orientation, info.Corners, seriesNumber, seriesDescription);
            pause(0.05);

        end

        inputFieldmapSlice = GERecon('Flex.GenerateFieldMap', echoChannelImages);

        inputFieldmap(:,:,slice) = inputFieldmapSlice;
    end

    % Flex Processing
    fieldmap = GERecon('Flex.ProcessFieldMap', inputFieldmap);

    water = zeros(size(fieldmap), 'single');
    fat = zeros(size(fieldmap), 'single');
    for slice = 1:pfile.reconstructedSlicesPerPass
        for echo = 1:pfile.echoes
            for channel = 1:pfile.channels

                % Transform K-Space
                channelImage = GERecon('Transform', kSpace(:,:,slice,echo,channel));

                echoChannelImages(:,:,echo,channel) = channelImage;
            end
        end
        [water(:,:,slice), fat(:,:,slice)] = GERecon('Flex.Separate', echoChannelImages, fieldmap(:,:,slice));
    end
    
    [water, fat] = GERecon('Flex.Identify', water, fat);
    
    for slice = 1:pfile.reconstructedSlicesPerPass    
        sliceInfo.sliceInPass = slice;
        info = GERecon('Pfile.Info', sliceInfo);

        % Apply Gradwarp, Orientation, Display and Save
        gradwarpedImage = GERecon('Gradwarp', water(:,:,slice), info.Corners);
        finalImage = GERecon('Orient', gradwarpedImage, info.Orientation);
        seriesNumber = header.SeriesData.se_no;
        seriesDescription = ['Water: ' header.SeriesData.se_desc];
        filename = ['DICOMs' num2str(seriesNumber) '/image' num2str(info.Number) '.dcm'];
        GERecon('Dicom.Write', filename, finalImage, info.Number, info.Orientation, info.Corners, seriesNumber, seriesDescription);

        % Apply Gradwarp, Orientation, Display and Save
        gradwarpedImage = GERecon('Gradwarp', fat(:,:,slice), info.Corners);
        finalImage = GERecon('Orient', gradwarpedImage, info.Orientation);
        seriesNumber = header.SeriesData.se_no * 100;
        seriesDescription = ['Fat: ' header.SeriesData.se_desc];
        filename = ['DICOMs' num2str(seriesNumber) '/image' num2str(info.Number) '.dcm'];
        GERecon('Dicom.Write', filename, finalImage, info.Number, info.Orientation, info.Corners, seriesNumber, seriesDescription);

        % Apply Gradwarp, Orientation, Display and Save
        orientedComplexImageReal = GERecon('Orient', real(fieldmap(:,:,slice)), info.Orientation);
        orientedComplexImageImag = GERecon('Orient', imag(fieldmap(:,:,slice)), info.Orientation);
        angleImage = uint16(65536 * (0.5 + angle(complex(orientedComplexImageReal, orientedComplexImageImag))./(2*pi)));
        finalImage = uint8(255 *ind2rgb(angleImage, hsv(65536)));
        imagesc(finalImage);
        seriesNumber = header.SeriesData.se_no * 100 + 3;
        seriesDescription = ['Fieldmap: ' header.SeriesData.se_desc];
        filename = ['DICOMs' num2str(seriesNumber) '/image' num2str(info.Number) '.dcm'];
        GERecon('Dicom.Write', filename, finalImage, info.Number, info.Orientation, info.Corners, seriesNumber, seriesDescription);

        pause(0.05);
    end
end

function loaded = LoadKacqWithInputString(kacqCell)
    loaded = 0;
    
    if isempty(kacqCell)
        fprintf('No kacq file specified\n');
    else
        % Load Arc Sampling Pattern (kacq_yz.txt)
        kacqFileString = char(kacqCell{1});
        if exist(kacqFileString, 'file')
            fprintf('Using specified kacq file : %s\n', kacqFileString);
            GERecon('Arc.LoadKacq', kacqFileString);
            loaded = 1;
        else
            fprintf('Could not find specified kacq file : %s\n', kacqFileString);
        end
    end
end

function loaded = LoadKacqBasedOnPfileLocation(pfilePath)
    loaded = 0;
    pfile = GERecon('Pfile.Load', pfilePath);
    header = GERecon('Pfile.Header', pfile);
    
    % No kacq file provided, check for kacq next to raw data file
    pathString = fileparts(pfilePath);
    potentialKacqFile = fullfile(pathString, ['kacq_yz.txt.' sprintf('%d', header.RawHeader.kacq_uid)]);
    fprintf('Attempting to load kacq based on %s\n', pfilePath);
    if exist(potentialKacqFile, 'file')
        fprintf('%s exists, and will be used for ARC processing\n', potentialKacqFile);
        GERecon('Arc.LoadKacq', potentialKacqFile);
        loaded = 1;
    else
        fprintf('%s was not found\n', potentialKacqFile);
    end
end
