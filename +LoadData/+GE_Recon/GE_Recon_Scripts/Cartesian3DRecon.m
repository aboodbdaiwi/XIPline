function Cartesian3DRecon(pfilePath)
%% CartesianRecon - Reconstruct 3D Cartesian K-Space
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
% Cartesian3DRecon(pfilePath)
% will reconstruct the 3D Cartesian K-Space in the given pfile. This 
% excludes pfiles with ARC enabled.
%
% Limitations: Parallel imaging, intensity correction

    % Load Pfile
    pfile = GERecon('Pfile.Load', pfilePath);
    header = GERecon('Pfile.Header', pfile);
    
    acquiredSlices = pfile.slicesPerPass;
    outputSlices = pfile.reconstructedSlicesPerPass;
    scaleFactor = outputSlices / pfile.scaleFactor3d;
   
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
            
            % Transform Across Slices
            kSpace = ifft(kSpace, outputSlices, 3);
            
            % Scale
            kSpace = kSpace * scaleFactor;
            
            % Loop for each slice/channel to create a magnitude image
            for slice = 1:outputSlices
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
                gradwarpedImage = GERecon('Gradwarp', magnitudeImage, info.Corners, 'XRMW');

                % Orient the image
                finalImage = GERecon('Orient', gradwarpedImage, info.Orientation);

                % Display
                imagesc(finalImage);

                % Display
                title(['Pass: ' num2str(pass) ' Slice: ' num2str(slice) ' Echo: ' num2str(echo)]);

                % Save DICOMs
                imageNumber = (info.Number-1) * pfile.echoes + echo;
                filename = ['DICOMs/image' num2str(imageNumber) '.dcm'];
                GERecon('Dicom.Write', filename, finalImage, imageNumber, info.Orientation, info.Corners);

                pause(0.05);
            end
        end
    end
end

