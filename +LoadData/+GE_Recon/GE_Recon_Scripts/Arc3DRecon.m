function Arc3DRecon(pfilePath, varargin)
%% Arc3DRecon - Reconstruct a 3D ARC Pfile
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
% Arc3DRecon(pfilePath, varargin)
% will reconstruct the 3D ARC data in the given pfile. The sampling pattern
% file (kacqXXX.txt) can be specified as the second argument to describe 
% the acquisition of the K-Space.
%
% Limitations: intensity correction

    % Load Pfile
    pfile = GERecon('Pfile.Load', pfilePath);
    header = GERecon('Pfile.Header', pfile);

    if length(varargin) == 1
        % Load Arc Sampling Pattern (kacq_yz.txt)
        GERecon('Arc.LoadKacq', varargin{1});
    end

    % Load KSpace. Since 3D Arc Pfiles contain space for the zipped
    % slices (even though the data is irrelevant), only pull out
    % the true acquired K-Space. Z-transform will zip the slices
    % out to the expected extent.
    acquiredSlices = pfile.slicesPerPass / header.RawHeader.zip_factor;
    
    % 3D Scaling Factor
    scaleFactor = header.RawHeader.user0;
    if header.RawHeader.a3dscale > 0
        scaleFactor = scaleFactor * header.RawHeader.a3dscale;
    end
    
    scaleFactor = pfile.slicesPerPass / scaleFactor;
    
    for pass = 1:pfile.passes
        for echo = 1:pfile.echoes
    
            kSpace = zeros(pfile.xRes, pfile.yRes, acquiredSlices, pfile.channels);

            for slice = 1:acquiredSlices
                
                sliceInfo.pass = pass;
                sliceInfo.sliceInPass = slice;
                
                for channel = 1:pfile.channels

                    % Load K-Space
                    kSpace(:,:,slice,channel) = GERecon('Pfile.KSpace', sliceInfo, echo, channel);

                end
            end

            % Synthesize KSpace
            kSpace = GERecon('Arc.Synthesize', kSpace);
            
            % Transform Across Slices
            kSpace = ifft(kSpace, pfile.slicesPerPass, 3);
            
            % Scale
            kSpace = kSpace * scaleFactor;

            % Loop for each slice/channel to create a magnitude image
            for slice = 1:pfile.slicesPerPass
                for channel = 1:pfile.channels

                    % Transform K-Space
                    channelImage = GERecon('Transform', kSpace(:,:,slice,channel));

                    channelImages(:,:,channel) = channelImage;
                end

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
                imageNumber = ImageNumber(pass, info.Number, echo, pfile);
                filename = ['DICOMs/image' num2str(imageNumber) '.dcm'];
                GERecon('Dicom.Write', filename, finalImage, imageNumber, info.Orientation, info.Corners);

                pause(0.05);
            end
        end
    end
end

function number = ImageNumber(pass, slice, echo, pfile)
% Image numbering scheme (P = Phase; S = Slice; E = Echo):
% P0S0E0, P0S0E1, ... P0S0En, P0S1E0, P0S1E1, ... P0S1En, ... P0SnEn, ...
% P1S0E0, P1S0E1, ... PnSnEn

    % Need to map the legacy "pass" number to a phase number
    numPassesPerPhase = fix(pfile.passes / pfile.phases);
    phase = fix(pass / numPassesPerPhase);

    slicesPerPhase = pfile.slicesPerPass * numPassesPerPhase * pfile.echoes;
    number = (phase-1) * slicesPerPhase + (slice-1) * pfile.echoes + (echo-1) + 1;
end
