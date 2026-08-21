function [ postageStamps ] = SpectroMultiChannelMultiVoxel( pfileFullPath, calData )
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
% Reconstruct multi-channel multi-voxel spectroscopy data. This script uses
% the GE product multi-coil combine algorithm. GE postage stamp pixel data
% will be generated and spectrums from each acquired voxel will be plotted.
%


    %% Load pfile
    % Load pfile and hold on to the pfile directory. The pfile directory
    % will be used to save postage stamp DICOM outputs.
    [pfilePath] = fileparts(pfileFullPath);
    pfileHandle = GERecon('Pfile.Load', pfileFullPath);
    header = GERecon('Pfile.Header', pfileHandle);
    
    % Extract spectroscopy parameters from the pfile. These parameters
    % include the excited voxel center and widths in RAS patient
    % coordinates.
    params = GERecon('Spectro.Mcsi.Parameters');        
            
    %% Load calibration
    % If a calibration pfile is supplied, the calibration processing will
    % be run on the pfile and the McsiCalibration file will be saved in the
    % calibration pfile's directory. The McsiCalibration.h5 file will then
    % be loaded from the calibration pfile directory.
    if(strcmp('.7', calData(end-1:end)))
        % Calibraiton pfile was supplied
        % Run calibration and read calibration from the cal pfile directory
        GERecon('Calibration.Process', calData, 'Standard');
        
        [calPfilePath] = fileparts(calData);
        calDataPath = fullfile(calPfilePath, filesep, 'McsiCalibration.h5');                    
        
        GERecon('Spectro.Mcsi.LoadCalibration', calDataPath); 
    else
        GERecon('Spectro.Mcsi.LoadCalibration', calData); 
    end    
    
    %% Channel by channel reconstruction
    % Allocate a matrix with space for the raw data for all slices and all
    % spatial voxels for a single channel. The data will be processed
    % channel by channel and accumulated into the GE multi-voxel
    % multi-channel channel combiner.
    singleChannelSortedFids = zeros(params.acquiredFidLength, params.acquiredXRes, params.acquiredYRes, params.acquiredZRes);    
    
    for channel=1:pfileHandle.channels
        %% Read in multi-slice, multi-voxel data from pfile
        % Multi-Voxel spectroscopy data is stored in a pfile on a 
        % channel by channel basis using a flat view index. That is, 
        % the slice and echo dimensions remain constant. The view 
        % index is incremented with each acquired FID. Thus, to access 
        % the first acquired FID, use slice = 1, echo = 1, view = 1. 
        % To access the next acquired FID, use 
        % slice = 1, echo = 1, view = 2, and so on.
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

        %% Chop along x-direction
        % chopX indicates to chop along x when it is true.
        % This is a negated version of the PSD RHTYPCHP bit which
        % indicates that rf chopping was turned on
        if(params.chopX)
            singleChannelSortedFids(:,2:2:end,:,:) = singleChannelSortedFids(:,2:2:end,:,:) .* -1.0;
        end
        
        %% Baseline correction
        % Inspect the end of each FID to determine an average baseline
        % offset. Subtract the average baseline offset from the raw
        % acquired data.
        singleChannelSortedFids = GERecon('Spectro.Mcsi.BaselineCorrection', singleChannelSortedFids);
        
        %% Orient
        % Rotate/Transpose the data according to the slice orientation
        % flags in the raw header for the pfile that is currently loaded.
        orientedData = GERecon('Spectro.Mcsi.Orient', singleChannelSortedFids);
       
        %% Zerofill
        % Zerofill data such that the spatial x/y dimensions have the 
        % same size and are equal to a power of 2. Example: An 8x10 
        % acquisition will be zerofilled to 16x16 spatial voxels. The FID 
        % spectral dimension is also zerofilled to twice the acquired 
        % FID length.
        zerofilledData = GERecon('Spectro.Mcsi.Zerofill', orientedData);
        
        %% Spectral transform 
        % Apply the spectral transform along the FID dimension of the
        % zerofilled data. A spectral filter is applied along
        % the FID dimension prior to transforming. The spectral filter can
        % be obtained with the Spectro.Mcsi.SpectralFilter command.
        if(channel == 0)
           % Display filter the first time through this loop
           spectralFilter = GERecon('Spectro.Mcsi.SpectralFilter');
           figure;plot(spectralFilter);title('Spectral Filter');drawnow;
        end
        spectralTransformedData = GERecon('Spectro.Mcsi.SpectralTransform', zerofilledData);
        
        %% Spatial Transform
        % Apply the spatial transform along the x,y,z spatial dimensions.
        % Note that this function does not apply any spatial apodization
        % filters.
        transformedData = GERecon('Spectro.Mcsi.SpatialTransform', spectralTransformedData);
               
        %% Individual Channel Plots
        % Plot intermediate data from 4.3ppm to -0.4ppm. Knowing that the
        % water peak is at 4.7ppm, the scan's center frequency
        % and the bandwidth of the acquisition, compute the indices to plot
        % -0.4ppm to 4.3ppm.
        
        % The number of voxels in the x/y directions are equal at this
        % point (because of the zerofilling step above). Thus, check only
        % dimension index 2 to obtain the centerVoxelIndex
        centerVoxelIndex = size(spectralTransformedData,2) / 2;
        
        if(size(spectralTransformedData,4) == 1)
           centerSliceIndex = 1; 
        else
           centerSliceIndex = size(spectralTransformedData,4) / 2; 
        end
        
        maxMagnitudeCenterVoxelCenterSlice = max(max(max(abs(transformedData(:,centerVoxelIndex,centerVoxelIndex,centerSliceIndex)))));
        tolerance = 0.0001;
        waterPeakIndex = find( abs(abs(transformedData(:,centerVoxelIndex,centerVoxelIndex,centerSliceIndex)) - maxMagnitudeCenterVoxelCenterSlice) < tolerance );        
                
        % Assume water is the largest peak and is at 4.7ppm. A shift of
        % 0.4ppm towards zero gets to 4.3ppm (the edge of the region that
        % will be plotted).
        zeroPointFourPpmShiftInHzFromWater = params.centerFrequencyMHz * 0.4;
        
        % 5.1ppm shift from water gets from 4.7ppm to -0.4ppm
        fivePointOnePpmShiftInHzFromWater = params.centerFrequencyMHz * 5.1;
        
        pointsPerHz = size(transformedData,1) / params.spectralBandwidthHz;
        startingPoint = ceil(waterPeakIndex + zeroPointFourPpmShiftInHzFromWater * pointsPerHz);
        endingPoint = ceil(waterPeakIndex + fivePointOnePpmShiftInHzFromWater * pointsPerHz);
        
        PlotSquareCuboid(real(transformedData(startingPoint:endingPoint,:,:,centerSliceIndex)), 'Channel Data', 1);                
        
        %% Accumulate data
        % Accumulate data using the GE multi-channel combine algorithm.
        % Note that this algorithm is specifically designed to be used 
        % with the product pulse sequence and processing steps above. This
        % algorithm will extract a subset of the spectrum to work with. The
        % algorithm is not compatible with 7.0T scans.
        % Part of the processing done in this function is to transform back
        % to kz and adjust the z-transform such that it yields slices at
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
        GERecon('Spectro.Mcsi.Accumulate', transformedData, channel);
    end
    
    %% Retrieve Accumulated Data
    % Retrieve the processed and accumulated data. The data returned here
    % will be a 256 point spectrum from ~ -0.4ppm to 4.3ppm. The ppm range
    % returned for plotting by GE multi-channel accumulation function is
    % not customizable at this time.
    accumulatedData = GERecon('Spectro.Mcsi.AccumulatedData');   
        
    %% Postage Stamp Generation    
    % Determine the number of reconstructed slices to make a postage stamp
    % for. The number of reconstructed slices is not equal to the number of
    % acquired slices because of the z-transform adjustment described
    % above. The slice dimension of the accumulated data is the fourth
    % dimension (SpectralDimension x X x Y x Slice).
    % Postage stamps are a 2D representation of spectral data for multiple
    % voxels of a single slice. A postage stamp always consists of a 16 x
    % 16 grid of squares. Each square contains one point of a 256 point
    % spectrum for each spatial voxel in the slice. Thus, for an 8x8 set of
    % spatial locations, each of the 16 squares has a size of 8x8. The
    % resulting image size is 16*8 x 16*8 = 128x128. For a 16x16 set of
    % spatial locations the resulting image size is 16*16 x 16*16 =
    % 256x256.     
    numReconStructedSlices = size(accumulatedData, 4);
       
    matlabDicomPath = fullfile(pfilePath, 'matlabDicoms', filesep);    
    imageNumber = 0;
    for slice = 1:numReconStructedSlices
        sliceString = ['Slice ' num2str(slice)];
        
        PlotSquareCuboid(imag(accumulatedData(:,:,:,slice)), [sliceString ' Imaginary Data']);drawnow;
        PlotSquareCuboid(real(accumulatedData(:,:,:,slice)), [sliceString ' Real Data']);drawnow;
                
        orientation = GERecon('Pfile.Orientation', slice);                    
        corners = GERecon('Pfile.Corners', slice);
                
        realImage = GERecon('Spectro.Mcsi.GeneratePostageStamp', real(accumulatedData(:,:,:,slice)));
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
        
        imaginaryImage = GERecon('Spectro.Mcsi.GeneratePostageStamp', imag(accumulatedData(:,:,:,slice)));
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

