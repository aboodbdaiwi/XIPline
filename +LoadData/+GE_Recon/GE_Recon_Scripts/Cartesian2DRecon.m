function Cartesian2DRecon(pfilePath)
%% CartesianRecon - Reconstruct 2D Cartesian K-Space
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
% Cartesian2DRecon(pfilePath)
% will reconstruct the 2D Cartesian K-Space in the given pfile.
%
% Limitations: Parallel imaging, intensity correction

    % Load Pfile
    pfile = GERecon('Pfile.Load', pfilePath);
    header = GERecon('Pfile.Header', pfile);
  
    % Reference views may be acquired for select 2D FSE scans. Reference views 
    % are used to correct for phase errors between readouts of an FSE echo train.
    % This example does not support FSE phase correction. Thus, the number of
    % reference views is computed here, and the reference views are discarded below.
    % Reference views are only acquired when the retrospective phase correction
    % algorithm is enabled. 
    fsePhaseCorrectionBitMask = 131072;
    numRefViews = 0;
    if( bitand(header.RawHeader.data_collect_type1, fsePhaseCorrectionBitMask) )
        numRefViews = header.RawHeader.etl;
    end 

    for p = 1:pfile.phases
        for s = 1:pfile.slices

            % Get corners and orientation for this slice location
            corners = GERecon('Pfile.Corners', s);
            orientation = GERecon('Pfile.Orientation', s);

            for e = 1:pfile.echoes
                for c = 1:pfile.channels

                    % Load K-Space
                    kspace = GERecon('Pfile.KSpace', s, e, c);
                    
                    % Discard reference views
                    fullYRes = size(kspace, 2);
                    yRes = fullYRes - numRefViews;
                    kspace = kspace(:, 1:yRes);

                    % Transform K-Space
                    channelImage = GERecon('Transform', kspace);

                    channelImages(:,:,c) = channelImage;
                end

                % Apply Channel Combination
                combinedImage = GERecon('SumOfSquares', channelImages);

                % Create Magnitude Image
                magnitudeImage = abs(combinedImage);

                % Apply Gradwarp
                gradwarpedImage = GERecon('Gradwarp', magnitudeImage, corners,'XRMW');

                % Orient the image
                finalImage = GERecon('Orient', gradwarpedImage, orientation);

                % Display
                imagesc(finalImage);
                title(['Phase: ' num2str(p) 'Slice: ' num2str(s) 'Echo: ' num2str(e)]);

                % Save DICOMs
                imageNumber = ImageNumber(s, e, p, pfile);
                filename = ['DICOMs/image' num2str(imageNumber) '.dcm'];
                GERecon('Dicom.Write', filename, finalImage, imageNumber, orientation, corners);

                pause(0.1);
            end
        end
    end
end

function number = ImageNumber(slice, echo, phase, pfile)
% Image numbering scheme:
% P0S0E0, P0S0E1, ... P0S0En, P0S1E0, P0S1E1, ... P0S1En, ... P0SnEn, ...
% P1S0E0, P1S0E1, ... PnSnEn
    slicesPerPhase = pfile.slices * pfile.echoes;
    number = (phase-1) * slicesPerPhase + (slice-1) * pfile.echoes + (echo-1) + 1;
end
