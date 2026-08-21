function Spiral2DRecon(pfilePath)
%% Spiral2DRecon - Reconstruct a 2D Spiral Pfile
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
% Spiral2DRecon(pfilePath)
% will reconstruct the 2D Spiral data in the given pfile.

    % Load Pfile
    pfile = GERecon('Pfile.Load', pfilePath);
    header = GERecon('Pfile.Header');
    [filePath] = fileparts(pfilePath);

    % Prep
    kSpace = zeros(pfile.xRes, pfile.yRes, pfile.channels);
    channelImages = zeros(pfile.imageSize, pfile.imageSize, pfile.channels);

    for s = 1:pfile.slices
        for c = 1:pfile.channels
            % Load K-Space
            kSpace(:,:,c) = GERecon('Pfile.KSpace', s, 1, c);
        end
        
        % Convert spiral K-Space to gridded K-space
        gridKSpace = GERecon('Spiral.Regrid', kSpace);
           
        for c = 1:pfile.channels
            % Transform K-Space
            channelImages(:,:,c) = ifft2(gridKSpace(:,:,c));
        end

        % Get corners and orientation for this slice location
        corners = GERecon('Pfile.Corners', s);
        orientation = GERecon('Pfile.Orientation', s);
            
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
        imagesc(maskedFinalImage);colormap(gray);colorbar;axis off;title(['Slice: ' num2str(s)]);
        drawnow;
        
        % Compute image number and create dicom image
        matlabDicomPath = fullfile(filePath, 'matlabDicom', filesep);
        fileName = [matlabDicomPath 'Image_' num2str(s-1, '%05d') '.dcm'];
        GERecon('Dicom.Write', fileName, maskedFinalImage, s, orientation, corners, header.SeriesData.se_no * 100);
    end
end
