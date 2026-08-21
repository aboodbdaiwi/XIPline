function [ postageStamps ] = SpectroPerChannelMultiVoxel(pfileFullPath)
%% SpectroMultiChannelMultiVoxel
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
% postageStamps = SpectroPerChannelMultiVoxel(pfileFullPath);
% 
% Reconstruct multi-voxel spectroscopy data on a per channel basis. 
% This script does not implement a channel combination. This script is
% equivalent to the GE Product 1.5T Prose Single Channel Reconstruction.


    %% Load pfile
    % Load pfile and hold on to the pfile directory. The pfile directory
    % will be used to save postage stamp dicom outputs.
    [pfilePath] = fileparts(pfileFullPath);
    pfileHandle = GERecon('Pfile.Load', pfileFullPath);
    header = GERecon('Pfile.Header', pfileHandle);
    
    % Extract spectroscopy parameters from the pfile. These parameters
    % include the excited voxel center and widths in RAS patient
    % coordinates.
    params = GERecon('Spectro.PerChannelMultiVoxel.Parameters');   
    
    %% Channel Loop
    % This script reconstructs each channel individually. The script does
    % not implement a channel combine and is used for single channel 1.5T
    % PROSE exams in product.
    imageNumber = 0;
    for channel=1:(pfileHandle.channels)
        %% Read in multi-slice, multi-voxel data from pfile
        % Multi-Voxel spectroscopy data is stored in a pfile on a 
        % channel by channel basis using a flat view index. That is, 
        % the slice and echo dimensions remain constant. The view 
        % index is incremented with each acquired FID. Thus, to access 
        % the first acquired FID, use slice = 0, echo = 0, view = 0. 
        % To access the next acquired FID, use 
        % slice = 0, echo = 0, view = 1, and so on.
        % Consider a hypothetical 3x4x3 multi-voxel scan
        % (acquiredXRes = frequency res on UI = 3), (acquiredYRes =
        % phase res on UI = 4), (acquiredZRes = num slices on UI =
        % 3)
        %
        % For 1-based indices:
        % flatViewIndex = (xIndex-1) + (yIndex-1)*(acquiredXRes) + 
        %                 (slice-1)*(acquiredXRes*acquiredYRes) + 1
        %
        % * Voxel Indices  ----------  Pfile Indices
        % * x=1, y=1, z=1  ----------  slice=1, echo=1, view=1
        % * x=2, y=1, z=1  ----------  slice=1, echo=1, view=2
        % * x=3, y=1, z=1  ----------  slice=1, echo=1, view=3
        % * x=1, y=2, z=1  ----------  slice=1, echo=1, view=4
        % * x=2, y=2, z=1  ----------  slice=1, echo=1, view=5
        % * x=3, y=2, z=1  ----------  slice=1, echo=1, view=6
        % * x=1, y=3, z=1  ----------  slice=1, echo=1, view=7
        % * ...
        % * x=1, y=4, z=3  ----------  slice=1, echo=1, view=34
        % * x=2, y=4, z=3  ----------  slice=1, echo=1, view=35
        % * x=3, y=4, z=3  ----------  slice=1, echo=1, view=36
        %     
        slice = 1;
        echo = 1;
        singleChannelAllFids = GERecon('Pfile.KSpace', slice, echo, channel);    
        for slice = 1:params.acquiredZRes
            for yIndex = 1:params.acquiredYRes
                for xIndex = 1:params.acquiredXRes                
                    flatViewIndex = (xIndex-1) + (yIndex-1)*(params.acquiredXRes) + (slice-1)*(params.acquiredXRes*params.acquiredYRes) + 1;

                    % Add one to flatViewIndex on this line to account for
                    % 1-based indexing in Matlab
                    singleChannelSortedFids(:,xIndex,yIndex,slice) = singleChannelAllFids(:, flatViewIndex);                
                end
            end
        end
        
        %% Z-direction Transform
        % The Multi-Voxel Spectroscopy Z direction transform yields slices at
        % slice locations with the same spacing as the localizer scan. This
        % is accomplished by adding a conversion factor to the exponent of
        % the Fourier transform. The conversion factor is stored in
        % rhuser48 and equals:
        %
        % $$rhuser48=\frac{(LocalizerSliceThickness + SliceSpacing)}{(CSISliceThickness)}$$
        %
        % Transform with conversion factor:
        %
        % $$f(n)=\sum\limits_{k=0}^{n-1}F(k)e^{-i\frac{2\pi}{N}kn(rhuser48)}$$
        %
        % Using the z-transform above results in each transform location
        % 'n' being spaced by (LocalizerSliceThickness + SliceSpacing). The
        % number of output slices can be obtained from the pfile's
        % sliceCount.        
        slicesToReconstruct = GERecon('Spectro.PerChannelMultiVoxel.ZTransform', singleChannelSortedFids);
        
        for reconstructedSlice = 1:size(slicesToReconstruct,4)
        
            %% Spectral and X/Y Spatial Direction Transforms
            % The 'Spectro.PerChannelMultiVoxel.Transform' command will apply
            % the spectral and XY spatial direction transform (the Z-direction
            % transform can be applied with the ZTransform command, see above).
            % This command will also apply baseline correction to remove any DC
            % offset in the FIDs.
            transformedSliceData = GERecon('Spectro.PerChannelMultiVoxel.Transform', slicesToReconstruct(:,:,:,reconstructedSlice));
           
            %% Extract ROI
            % Extract the region of interest to include in postage stamps.
            % This command will interpolate the spectrum such that 256
            % points are extracted between ~4.3ppm and 0.4ppm. This command
            % will adjust the extracted region of the input spectrum based on
            % the temperature value in the currently active pfile's header.
            roiToPlot = GERecon('Spectro.PerChannelMultiVoxel.ExtractRoi', transformedSliceData);
            
            %% Orient
            % Orient the spectrums according to the rotate and transpose values
            % specified in the currently active pfile's header.
            orientedSpectrum = GERecon('Spectro.PerChannelMultiVoxel.Orient', roiToPlot);
            
            %% Postage Stamp Generation    
            % Postage stamps are a 2D representation of spectral data for multiple
            % voxels of a single slice. A postage stamp always consists of a 16 x
            % 16 grid of squares. Each square contains one point of a 256 point
            % spectrum for each spatial voxel in the slice. Thus, for an 8x8 set of
            % spatial locations, each of the 16 squares has a size of 8x8. The
            % resulting image size is 16*8 x 16*8 = 128x128. For a 16x16 set of
            % spatial locations the resulting image size is 16*16 x 16*16 =
            % 256x256.     
            matlabDicomPath = fullfile(pfilePath, 'matlabDicoms', filesep);    
            
            sliceString = ['Slice ' num2str(reconstructedSlice)];
            PlotSquareCuboid(imag(orientedSpectrum), [sliceString ' Imaginary Data']);drawnow;
            PlotSquareCuboid(real(orientedSpectrum), [sliceString ' Real Data']);drawnow;

            orientation = GERecon('Pfile.Orientation', reconstructedSlice);                    
            corners = GERecon('Pfile.Corners', reconstructedSlice);

            realImage = GERecon('Spectro.PerChannelMultiVoxel.GeneratePostageStamp', real(orientedSpectrum));
            if(imageNumber < 10)
                imageNumberString = ['00' num2str(imageNumber)];
            elseif(imageNumber < 100)
                imageNumberString = ['0' num2str(imageNumber)];
            else
                imageNumberString = num2str(imageNumber);
            end                
            GERecon('Dicom.Write', [matlabDicomPath 'Image_' imageNumberString '.dcm'], realImage, imageNumber, orientation, corners, (header.SeriesData.se_no*100));                        
            imageNumber = imageNumber + 1;
            postageStamps(:,:,imageNumber) = realImage;

            imaginaryImage = GERecon('Spectro.PerChannelMultiVoxel.GeneratePostageStamp', imag(orientedSpectrum));
            if(imageNumber < 10)
                imageNumberString = ['00' num2str(imageNumber)];
            elseif(imageNumber < 100)
                imageNumberString = ['0' num2str(imageNumber)];
            else
                imageNumberString = num2str(imageNumber);
            end            
            GERecon('Dicom.Write', [matlabDicomPath 'Image_' imageNumberString '.dcm'], imaginaryImage, imageNumber, orientation, corners, (header.SeriesData.se_no*100));        
            imageNumber = imageNumber + 1;             
            postageStamps(:,:,imageNumber) = imaginaryImage;
        end                               
    end     
end

function PlotSquareCuboid( squareCuboidInput, figureWindowTitle, localScale)
%% PlotSquareCuboid
% Plot the spectrums in the given cuboid (Spectrum x X x Y) to a 2D figure
% window. This function creates a set of MATLAB axes handles up front and
% then plots a spectrum in each axes. The third parameter is optional and
% specified whether the scale of each individual plot is scaled based on all
% voxels or only the data in a single voxel/individual plot. The default is
% to apply a global scale (localScale == 0) such that all plots are scaled
% uniformly
    
    if(size(squareCuboidInput, 2) ~= size(squareCuboidInput, 3))
       disp('Size mismatch, cuboid must be square!');
       return;
    end
    
    figure('name', figureWindowTitle);

    % Input data is expected to be square. Thus, obtain the num voxels in
    % each direction of the cuboid from the size of the second dimension.
    numVoxelsInOneDirection = size(squareCuboidInput, 2);    
                
    % Initialize a set of axes objects (one for each spectrum to be
    % plotted).
    axesIndex = 1;
    for i=1:numVoxelsInOneDirection
      xAxesLocation = (i-1)/numVoxelsInOneDirection;
      for j=1:numVoxelsInOneDirection
        yAxesLocation = (j-1)/numVoxelsInOneDirection;
        
        % Set units to normalized which maps the lower left corner of the
        % figure window to 0,0 and the upper right to 1,1
        % Then, set the position in the current window to the given x,y
        % location and set the width/height to 1/numVoxelsInOneDirection 
        squareHeightWidth = 1.0/numVoxelsInOneDirection;
        axesCollection(axesIndex) = axes('units','norm','pos',[xAxesLocation yAxesLocation squareHeightWidth squareHeightWidth]);
        axesIndex = axesIndex + 1;
      end
    end
    
    % Remove tick labels and set next plot to be added to the current
    % figure handle;
    set(axesCollection,'XtickLabel','','YTickLabel','','nextplot','add');    
    
    % Determine global min/max values for this slice. Scale all spectrums
    % to this global min/max.
    minVal = min(min(min(squareCuboidInput)));
    maxVal = max(max(max(squareCuboidInput)));
    
    if(nargin < 3)
       localScale = 0; 
    end
    
    % Postage stamps are oriented for display in FuncTool. Match FuncTool
    % display to the display generated with this code.
    i = 1;
    for y=1:numVoxelsInOneDirection
        for x=numVoxelsInOneDirection:-1:1
            kk = i;       
            
            plot(squareCuboidInput(:,y,x),'parent',axesCollection(kk));            
            
            xlim(axesCollection(kk), [1 size(squareCuboidInput,1)]);        
            
            if(localScale > 0)
                % Reset max/min values to the max/min from just this
                % spatial location
                minVal = min(min(min(squareCuboidInput(:,y,x))));
                maxVal = max(max(max(squareCuboidInput(:,y,x))));                        
            end
            
            if(minVal < maxVal)              
                ylim(axesCollection(kk), [minVal maxVal]);
            else
                ylim(axesCollection(kk), [-2000 2000]);                                   
            end
            
            i = i + 1;            
        end
    end
end
