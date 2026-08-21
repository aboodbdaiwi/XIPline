function [ finalImages ] = EpiMultiPhaseRecon(pfileFullPath)
%% EpiMultiPhaseRecon
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
%   Reconstruct multi-phase Epi data. This script takes additional
%   reference views acquired at the beginning of an Epi echo train into
%   account. These views are used to dynamically update phase correction
%   coefficients throughout a scan. The additional reference views can be
%   acquired using the GE fMRI (Brainwave) or Non-Brainwave sequences.
%   This script also reconstructs Hyperband data into multiple acquired
%   slices.
%
%   EpiMultiPhaseRecon(pfile)
%   will reconstruct the multi-phase Epi data in the pfile. The script will
%   look in the following locations for all secondary inputs.
%   Reference data - First look for ref.dat in the pfile directory. If
%   ref.dat cannot be found, look for a reference pfile in a 'ref'
%   subdirectory (fileFullPath/ref/Prrrrr.7) where Prrrrr.7 is the
%   reference pfile and the reference pfile run number is equal to the scan
%   pfile run number.
%   Sinc Interpolation Kernels - First look for vrgf.dat in the pfile
%   directory. If vrgf.dat cannot be found, look for vrgf_kernels.dat in
%   the pfile directory.
%   Asset Calibration - First look for AssetCalibration.h5 in the pfile
%   directory. If AssetCalibration.h5 cannot be found then use Sum of
%   Squares channel combination.
%   Row Flip Parameters - Look for rowflip.param in the pfile directory.
%   If found, apply rowflip to reference views (image views are already
%   flipped in the pfile). If not found, do not apply rowflip to reference
%   views. The reference view rowflip is not essential but is included for
%   completeness.

%% Determine secondary input locations
% If the number of arguments supplied is greater than 1 then the user
% supplied absolute locations for all secondary inputs. If the number
% of arguments is equal to 1 then look in the pfile directory for
% secondary inputs.
[pfilePath pfileName pfileExt] = fileparts(pfileFullPath);

% Look for reference data
refDotDatPath = fullfile(pfilePath, 'ref.dat');
refDotH5Path =  fullfile(pfilePath, 'ref.h5');
if(exist(refDotDatPath, 'file') == 2)
    % Found reference data, use ref.dat in pfile directory
    referenceData = refDotDatPath;
elseif(exist(refDotH5Path, 'file') == 2)
    referenceData = refDotH5Path;
else
    % ref.dat is not in pfile directory, look for reference pfile
    % in ref sub-directory
    refPfilePath = fullfile(pfilePath, 'ref', [pfileName pfileExt]);
    if(exist(refPfilePath, 'file') == 2)
        % Found reference data, use reference pfile
        referenceData = refPfilePath;
    else
        % Could not find reference data!
        disp('Could not find reference data.');
        return;
    end
end

% Look for vrgf data
vrgfDotDatPath = fullfile(pfilePath, 'vrgf.dat');
if(exist(vrgfDotDatPath, 'file') == 2)
    % Found vrgf.dat, use for vrgf interpolation kernels
    vrgfInterpKernels = vrgfDotDatPath;
else
    % Look for vrgf_kernels.dat
    vrgfKernelsDotDatPath = fullfile(pfilePath, 'vrgf_kernels.dat');
    if(exist(vrgfKernelsDotDatPath, 'file') == 2)
        vrgfInterpKernels = vrgfKernelsDotDatPath;
    else
        % Could not find vrgf interpolation kernels!
        disp('Could not find vrgf interpolation kernels.');
        return;
    end
end

% Look for AssetCalibration. Note that if you have a calibration pfile
% you can generate the AssetCalibration.h5 file using this command:
%
%   GERecon('Calibration.Process', 'pathToCalPfile/Pxxxxx.7').
%
% The calibration HDF5 file will be saved in the calibration pfile's
% directory.
assetCalibrationPath = fullfile(pfilePath, 'AssetCalibration.h5');
if(exist(assetCalibrationPath, 'file') == 2)
    useAsset = 1;
    assetCalibration = assetCalibrationPath;
else
    useAsset = 0;
end

% Look for RawCalibration. Note that if you have a calibration pfile
% you can generate the RawCalibration.h5 file using this command:
%
%   GERecon('Calibration.Process', 'pathToCalPfile/Pxxxxx.7').
%
% The calibration HDF5 file will be saved in the calibration pfile's
% directory.
rawCalibrationPath = fullfile(pfilePath, 'RawCalibration.h5');
if(exist(rawCalibrationPath, 'file') == 2)
    useHyperband = 1;
    rawCalibration = rawCalibrationPath;
else
    useHyperband = 0;
end

% This flag may be used to control the number of figures plotted during
% recon
verboseFigures = 1;

%% Load reference data
% Inspect the reference parameter to determine if a ref.dat file was
% provided or if a reference pfile was provided.
% If a ref.dat file is provided, load the ref.dat file into a phase
% correction reference handle.
% If a reference pfile is provided, use the pfile to compute phase
% correction coefficients. The coefficients are stored by the
% <matlab:doc('EpiComputeAndApplyPhaseCorrection') EpiComputeAndApplyPhaseCorrection> script in a phase correction
% reference handle which is returned to this caller.
% The phase correction reference handle will be used later in this
% script to apply phase correction.
if(strcmp('.7', referenceData(end-1:end)))
    % Reference pfile was provided, load that and compute coefficients
    phaseCorrectionHandle = EpiComputeAndApplyPhaseCorrection(referenceData);
    
    % Load pfile after loading and computing coefficients from ref pfile
    pfileHandle = GERecon('Pfile.Load', pfileFullPath);
    header = GERecon('Pfile.Header', pfileHandle);
else
    % Load pfile before loading ref.dat
    pfileHandle = GERecon('Pfile.Load', pfileFullPath);
    header = GERecon('Pfile.Header', pfileHandle);
    
    % ref.dat was provided, load it
    phaseCorrectionHandle = GERecon('Epi.LoadReference', referenceData);
end

%% Load VRGF interpolation kernels
% VRGF (variable readout gradient filtering) refers to sampling data on
% the gradient ramps during data acquisition. During recon, a sinc
% interpolation is performed to interpolate non-linearly sampled
% frequency data (i.e. data sampled on the gradient ramps) to linearly
% sampled frequency data.
% The sinc interpolation kernels that are used for each point in the
% interpolated output are contained in a vrgf.dat or vrgf_kernels.dat
% file. The vrgf.dat file contains a full sinc function (i.e. sinc
% function length = readout direction acquisition size). The
% vrgf_kernels.dat file contains a truncated sinc interpolation kernel
% (i.e. sinc function length < readout direction acquisition size).
% Either the full or truncated sinc functions may be used for recon. To
% match product functionality, use vrgf_kernels.dat for fMRI scans and
% vrgf.dat for all other scans.
vrgfHandle = GERecon('RampSampling.Load', vrgfInterpKernels);
kernelMatrix = GERecon('RampSampling.KernelMatrix', vrgfHandle);
figure;
subplot(2,1,1);
imagesc(kernelMatrix);colormap(gray);colorbar;axis off;title('VRGF Interpolation Kernel Matrix');
subplot(2,1,2);
numRows = size(kernelMatrix, 1);
plot(kernelMatrix(numRows/4, :));title(['Sinc Interpolation Kernel for Interpolated Point: ' num2str(numRows/4) ' of ' num2str(numRows)]);

%% Load Asset Calibration (if applicable)
% If the user specified an Asset calibration file and asset is enabled
% for this scan, read the asset calibration.
assetEnabled = (header.RawHeader.asset == 2 || header.RawHeader.asset == 7);
if((useAsset == 1) || assetEnabled)
    GERecon('Asset.LoadCalibration', assetCalibration);
end

%% Load Raw Calibration (if applicable)
% If the user specified an Raw calibration file and Hyperband is enabled
% for this scan, read the Raw calibration.
hyperbandEnabled = (header.RawHeader.mb_factor > 1);
multiBandFactor = header.RawHeader.mb_factor;
if((useHyperband == 1) && hyperbandEnabled)
    multibandCalHandle = GERecon('Epi.Hyperband.LoadCalibration', rawCalibration);
end

%% Load rowflip.param (if applicable)
% If a rowflip filepath was specified, load the rowflip file.
% The reference views acquired for use by dynamic
% phase correction are not flipped in the pfile. Image views are
% flipped in the pfile. Dynamic phase correction still works if the
% reference views are left un-flipped. The reference view flip is
% included here for completeness.
pfilePath = fileparts(pfileFullPath);
rowFlipFilePath = fullfile(pfilePath, 'rowflip.param');
if(exist(rowFlipFilePath, 'file') == 2 && ~pfileHandle.areReferenceViewsFlipped)
    rowFlipHandle = GERecon('Epi.LoadRowFlip', rowFlipFilePath);
    applyReferenceRowFlip = 1;
else
    applyReferenceRowFlip = 0;
end

%% Initialize workspace and header variables
% Pull values from the raw header and use these values to allocate
% space for intermediate data.
imageSize = header.RawHeader.im_size;
refViewsTop = header.RawHeader.extra_frames_top;
refViewsBottom = header.RawHeader.extra_frames_bot;
totalNumRefViews = refViewsTop + refViewsBottom;
totalImageViews = header.RawHeader.da_yres - totalNumRefViews - 1;

if (hyperbandEnabled)
    numAcquiredSlices = ceil(pfileHandle.slices/multiBandFactor);
else
    numAcquiredSlices = pfileHandle.slices;
end

% Allocate intermediate data space
channelImages = single(zeros(imageSize, imageSize, pfileHandle.channels));
finalImages = int16(zeros(imageSize, imageSize, numAcquiredSlices));
rawImageViews = single(zeros(header.RawHeader.da_xres, totalImageViews, pfileHandle.channels));
rawReferenceViews = single(zeros(header.RawHeader.da_xres, totalNumRefViews, pfileHandle.channels));

% Create figures and allocate dynamic coefficient space for visualizing
% dynamic phase correction coefficient changes over the course of a
% multi-phase scan
imagesFigureHandle = figure;
dynamicCoefficientsFigureHandle = figure;
updatedCoefficients = zeros(totalImageViews, pfileHandle.channels, 2);
constantDynamicCoefficients = zeros(pfileHandle.phases, 2);
linearDynamicCoefficients = zeros(pfileHandle.phases, 2);

% Initialize image number to zero and begin looping over all phases
% and slices in the scan
imageNumber = 0;
for phase = 1:(pfileHandle.phases)
    for logicalSlice = 1:numAcquiredSlices
        %% Extract Raw Data
        % Dynamic phase correction works across all channels, thus
        % extract all channel data for the current slice (both
        % reference and image views) into rawImageViews and
        % rawReferenceViews matrices.
        for channel = 1:(pfileHandle.channels)
            if(hyperbandEnabled)
                geometricSliceIndex = header.DataAcqTable.slice_in_pass(logicalSlice);
            else
                geometricSliceIndex = logicalSlice;
            end
            kSpace = zeros(pfileHandle.xRes, pfileHandle.yRes, logicalSlice, pfileHandle.channels);
            % Echo dimension is unused for baseic Epi scans
            echo = 1;
            if(hyperbandEnabled)
                sliceInfo.pass = phase;
                sliceInfo.sliceInPass = logicalSlice;
                kSpace = GERecon('Pfile.KSpace', sliceInfo, echo, channel);
            else
                kSpace = GERecon('Pfile.KSpace', logicalSlice, echo, channel, phase);
            end
            
            kSpaceTotalNumViews = size(kSpace,2);
            rawImageViews(:,:,channel) = kSpace(:,(refViewsTop+1):(kSpaceTotalNumViews-refViewsBottom));
            
            if(refViewsTop > 0)
                refViews = kSpace(:,1:refViewsTop);
            else
                refViews = kSpace(:,(kSpaceTotalNumViews-refViewsBottom+1):end);
            end
            
            % Reference views are not flipped in a pfile while image
            % views are flipped. Thus, flip reference views here
            if(applyReferenceRowFlip == 1 && totalNumRefViews > 0)
                rowFlippedRefViews = GERecon('Epi.ApplyReferenceRowFlip', rowFlipHandle, refViews);
                rawReferenceViews(:,:,channel) = rowFlippedRefViews;
            else
                rawReferenceViews(:,:,channel) = refViews;
            end
        end
        
        %% Apply Phase Correction
        % If reference views were acquired, apply dynamic phase
        % correction. If no reference views were acquired default to
        % static phase correction using coefficients from the reference
        % scan.
        if(totalNumRefViews > 0)
            % To apply dynamic phase correction pass the current phase
            % correction handle and slice index along with the
            % rawImageViews and  rawReferenceViews matrices (for the
            % current slice) to the dynamic phase correction command. The
            % static coefficients for the current slice (computed during
            % reference scan recon) will be extracted from the phase
            % correction handle and updated based on the rawImageViews and
            % rawReferenceViews.
            % The first time that the dynamic phase correction command is
            % called for each slice, the rawImageViews and
            % rawReferenceViews that are provided will be cached as the
            % baseline phase. All future calls to dynamic phase correction
            % for the given slice will compare the rawImageViews and
            % rawReferenceViews to the cached baseline phase to compute
            % updated phase correction coefficients.
            % The updated coefficients that are applied to the
            % rawImageViews are optionally returned (if the user specifies
            % a second return argument). Note that for product Epi scans
            % with Asset, the same coefficients are applied to all channels
            % of a given slice. For this case, the updated coefficients
            % returned to the user have a channel dimension size of 1.
            [phaseCorrectedKSpace, updatedCoefficients] = GERecon('Epi.ApplyDynamicPhaseCorrection', phaseCorrectionHandle, rawImageViews, rawReferenceViews, geometricSliceIndex);
        else
            % Apply static phase correction for each channel based on
            % the coefficients held in the phase correction handle
            for channel = 1:pfileHandle.channels
                phaseCorrectedKSpace(:,:,channel) = GERecon('Epi.ApplyPhaseCorrection', phaseCorrectionHandle, rawImageViews(:,:,channel), geometricSliceIndex, channel);
                currentChannelCoefs = GERecon('Epi.RetrieveCoefficients', phaseCorrectionHandle, geometricSliceIndex, channel);
                updatedCoefficients(:,channel,:) = currentChannelCoefs;
            end
        end
        
        %% Visualize dynamic coefficients
        % Plot the phase correction coefficients applied for each
        % phase. If dynamic phase correction is used, these
        % coefficients will change with each phase of the scan. If
        % dynamic phase correction is not used then these coefficients
        % will remain constant throughout the scan.
        % For demonstration purposes, plot the linear and constant
        % coefficients for the two center views of slice index 0. Note
        % that each view is shifted toward the center of kSpace. Thus,
        % one view shifts left while an adjacent view shifts right. The
        % shift is accomplished by applying a phase ramp along each
        % readout in x,ky space and then transforming back to kx,ky
        % space. Because each view shifts towards the center of kSpace,
        % the signs on coefficients for adjacent views are negated.
        % Dynamic phase correction is only supported for single shot
        % Epi scans.
        if(verboseFigures == 1 && (logicalSlice == floor(pfileHandle.slices / 2)))
            constantIndex = 1; % Constant coefficients are at third index 1 in the updatedCoefficients matrix
            linearIndex = 2; % Linear coefficients are at third index 2 in the updatedCoefficients matrix
            firstChannel = 1; % Only plot the first channel for visualization purposes
            centerView = size(updatedCoefficients, 1) / 2; % Only plot the center view for visualization purposes
            centerViewPlusOne = centerView + 1;
            constantDynamicCoefficients(phase, 1) = updatedCoefficients(centerView, firstChannel, constantIndex);
            linearDynamicCoefficients(phase, 1) = updatedCoefficients(centerView, firstChannel, linearIndex);
            constantDynamicCoefficients(phase, 2) = updatedCoefficients(centerViewPlusOne, firstChannel, constantIndex);
            linearDynamicCoefficients(phase, 2) = updatedCoefficients(centerViewPlusOne, firstChannel, linearIndex);
            
            figure(dynamicCoefficientsFigureHandle);
            subplot(2,2,1);
            plot(constantDynamicCoefficients(1:phase,1), '-o');title(['Constant Coefficients, View: ' num2str(centerView) ' (radians)']);xlabel('Phase Index');
            subplot(2,2,2);
            plot(constantDynamicCoefficients(1:phase,2), '-o');title(['Constant Coefficients, View: ' num2str(centerViewPlusOne) ' (radians)']);xlabel('Phase Index');
            subplot(2,2,3);
            plot(linearDynamicCoefficients(1:phase,1), '-o');title(['Linear Coefficients, View: ' num2str(centerView) ' (radians/pixel)']);xlabel('Phase Index');
            subplot(2,2,4);
            plot(linearDynamicCoefficients(1:phase,2), '-o');title(['Linear Coefficients, View: ' num2str(centerViewPlusOne) ' (radians/pixel)']);xlabel('Phase Index');
            drawnow;
        end
        
        %% For Hyperband scans
        if(hyperbandEnabled)
            GERecon('Epi.Hyperband.Calibrate',multibandCalHandle,logicalSlice);
            inPlaneAcceleration = 1/header.RawHeader.asset_R;
            inplaneUnpackedYRes = totalImageViews * inPlaneAcceleration;
            resampledKSpaceToUnalias = zeros(size(kernelMatrix,1),inplaneUnpackedYRes,pfileHandle.channels);
            for channel = 1:pfileHandle.channels
                % Sinc Interpolation for Ramp Sampled Data
                % Interpolate ramp sampled data to linearly spaced samples
                % in frequency space.
                kspaceToInterpolate = phaseCorrectedKSpace(:,:,channel);
                resampledKSpaceToUnalias(:,(1:inPlaneAcceleration:end),channel) = GERecon('RampSampling.Interpolate', vrgfHandle, kspaceToInterpolate);
            end
            
            [multibandUnaliasedSliceData,unaliasedGeometricSliceNumbers] = GERecon('Epi.Hyperband.Unalias', resampledKSpaceToUnalias,logicalSlice);
            geometricSliceLoopcount = size(unaliasedGeometricSliceNumbers,2);
        else
            geometricSliceLoopcount = 1;
        end
        
        for unaliasedSliceIndex = 1:geometricSliceLoopcount
            if (hyperbandEnabled)
                geometricSliceNumber = unaliasedGeometricSliceNumbers(unaliasedSliceIndex);
            else
                geometricSliceNumber = logicalSlice;
            end
            
            if (geometricSliceNumber <= pfileHandle.slices) % To discard extended slices for Multiband scans
                corners = GERecon('Pfile.Corners', geometricSliceNumber);
                orientation = GERecon('Pfile.Orientation', geometricSliceNumber);
                for channel = 1:(pfileHandle.channels)
                    if(hyperbandEnabled)
                        interpolatedData = multibandUnaliasedSliceData(:,:,unaliasedSliceIndex,channel);
                    else
                        phaseCorrectedKSpaceCurrentChannel = phaseCorrectedKSpace(:,:,channel);
                        interpolatedData = GERecon('RampSampling.Interpolate', vrgfHandle, phaseCorrectedKSpaceCurrentChannel);
                    end
                    
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
                    
                    channelImages(:,:,channel,:) = image;
                end
                
                %% Channel combination
                % Combine channel images using either Asset or a Sum of Squares
                % channel combination.
                if(useAsset)
                    channelCombinedImage = GERecon('Asset.Unalias', channelImages, corners);
                else
                    channelCombinedImage = GERecon('SumOfSquares', channelImages);
                end
                
                %% Image Finalization
                % Zero out kissoff views and apply gradwarp.
                kissoffViews = header.RawHeader.kissoff_views;
                channelCombinedImage(:,1:kissoffViews) = 0;
                channelCombinedImage(:,(end-kissoffViews+1):end) = 0;
                gradwarpImage = GERecon('Gradwarp', abs(channelCombinedImage), corners, 'XRMW');
                
                % Apply final scaling for partial ky scans (this scaling
                % matches product behavior)
                if(header.RawHeader.hnover > 0)
                    gradwarpImage = gradwarpImage * (256 / (header.RawHeader.rc_xres * header.RawHeader.rc_yres));
                end
                
                % Rotate/Transpose
                rotatedTransposedSlice = GERecon('Orient', gradwarpImage, orientation);
                rotatedTransposedCorners = GERecon('Orient', corners, orientation);
                
                % Clip to range of shorts, convert to short, and write Dicom
                % image
                rotatedTransposedSlice(rotatedTransposedSlice < 0) = 0;
                rotatedTransposedSlice(rotatedTransposedSlice > 32767) = 32767;
                finalImages(:,:,geometricSliceNumber) = int16(rotatedTransposedSlice);
                
                if( verboseFigures == 1 || ( (phase == pfileHandle.phases) && (geometricSliceNumber == pfileHandle.slices) ) )
                    figure(imagesFigureHandle);
                    imagesc(finalImages(:,:,geometricSliceNumber));colormap(gray);colorbar;axis off;title(['Geometric Slice: ' num2str(geometricSliceNumber) ' Phase: ' num2str(phase)]);
                    drawnow;
                end
                imageNumber = (phase- 1) * pfileHandle.slices  + (geometricSliceNumber -1) ;
                
                if(imageNumber < 10)
                    imageNumberString = ['000' num2str(imageNumber)];
                elseif(imageNumber < 100)
                    imageNumberString = ['00' num2str(imageNumber)];
                else
                    imageNumberString = num2str(imageNumber);
                end
                
                matlabDicomPath = fullfile(pfilePath, 'matlabDicoms', filesep);
                GERecon('Dicom.Write', [matlabDicomPath 'Image_' imageNumberString '.dcm'], finalImages(:,:,geometricSliceNumber), imageNumber, orientation, rotatedTransposedCorners, (header.SeriesData.se_no*100));
                
                
            end
        end % End Geometric Slice Loop
    end% End Logical Slice Loop
end % End Phase Loop
end

