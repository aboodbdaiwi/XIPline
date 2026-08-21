function [ finalImages ] = EpiRecon(pfileFullPath)
%% EpiRecon - Reconstruct EPI Images
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
%   EpiRecon(pfileFullPath)
%   will reconstruct the Epi data in the given pfile. The script will
%   look for all required secondary inputs in the pfile directory.
%   Reference data - First look for ref.dat in the pfile directory. If
%   ref.dat cannot be found, look for a reference pfile in a 'ref'
%   subdirectory (fileFullPath/ref/Prrrrr.7) where Prrrrr.7 is the
%   reference pfile and the reference pfile run number is equal to the scan
%   pfile run number.
%   Sinc Interpolation Kernels - First look for vrgf.dat in the pfile
%   directory. If vrgf.dat cannot be found, look for vrgf_kernels.dat in
%   the pfile directory.
%   Asset Calibration - Look for AssetCalibration.h5 in the pfile
%   directory. If AssetCalibration.h5 cannot be found then use Sum of
%   Squares channel combination.
%

    %% Determine secondary input locations
    % If the number of arguments supplied is greater than 1 then the user
    % supplied absolute locations for all secondary inputs. If the number
    % of arguments is equal to 1 then look in the pfile directory for
    % secondary inputs.
    [pfilePath pfileName pfileExt] = fileparts(pfileFullPath);
       
    % Look for reference data
    refDotDatPath = fullfile(pfilePath, 'ref.dat');
    if(exist(refDotDatPath, 'file') == 2)
        % Found reference data, use ref.dat in pfile directory
        referenceData = refDotDatPath;
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
    
    %% Load reference data
    % Inspect the reference parameter to determine if a ref.dat file was
    % provided or if a reference pfile was provided.
    % If a ref.dat file is provided, load the ref.dat file into a phase 
    % correction reference handle.
    % If a reference pfile is provided, use the pfile to compute phase
    % correction coefficients. The coefficients are stored by the
    % EpiComputeAndApplyPhaseCorrection script in a phase correction
    % reference handle which is returned to this caller.
    % The phase correction reference handle will be used later in this
    % script to apply phase correction.
    if(strcmp('.7', referenceData(end-1:end)))
        % Reference pfile was provided, load that and compute coefficients
        phaseCorrectionHandle = EpiComputeAndApplyPhaseCorrection(referenceData);    
        
        % Load pfile after loading and computing coefficients from ref pfile
        pfileHandle = GERecon('Pfile.Load', pfileFullPath);
    else
        % Load pfile before loading ref.dat
        pfileHandle = GERecon('Pfile.Load', pfileFullPath);
        
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
    % If the user specified an Asset calibration file, load the file and
    % use Asset for channel combination.
    if(useAsset == 1)        
        GERecon('Asset.LoadCalibration', assetCalibration);    
    end
    
    %% Initialize data storage
    % Extract information from pfile raw header and allocate space for 
    % reconstructed data.
    header = GERecon('Pfile.Header');
    imageSize = header.RawHeader.im_size;
    refViewsTop = header.RawHeader.extra_frames_top;
    refViewsBottom = header.RawHeader.extra_frames_bot;
    channelImages = single(zeros(imageSize, imageSize, pfileHandle.channels));    
    finalImages = int16(zeros(imageSize, imageSize, pfileHandle.slices));
        
    %% Loop and Recon
    % Loop over all slices/channels in the scan (echo dimension is
    % unused in basic Epi reconstructions).
    figure;    
    imageNumber = 0;    
    for slice = 1:pfileHandle.slices
        for channel = 1:pfileHandle.channels
            % Echo dimension is unused for basic Epi scans and this script
            % does not support multi-phase scans (see EpiMultiPhaseRecon)
            echo = 1;            
            phase = 1;
            kSpace = GERecon('Pfile.KSpace', slice, echo, channel, phase);
            corners = GERecon('Pfile.Corners', slice);
            orientation = GERecon('Pfile.Orientation', slice);

            % Account for additional reference views that may be
            % acquired with fMRI data
            kSpaceTotalNumViews = size(kSpace,2);                
            kSpaceWithoutRefViews = kSpace(:,(refViewsTop+1):(kSpaceTotalNumViews-refViewsBottom));

            % Apply phase correction
            phaseCorrectedKSpace = GERecon('Epi.ApplyPhaseCorrection', phaseCorrectionHandle, kSpaceWithoutRefViews, slice, channel);

            % Interpolate ramp sampled data
            interpolatedData = GERecon('RampSampling.Interpolate', vrgfHandle, phaseCorrectedKSpace);

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
            
            channelImages(:,:,channel, :) = image;
        end

        % ASSET Unalias
        if(useAsset)
            channelCombinedImage = GERecon('Asset.Unalias', channelImages, corners);
        else                
            channelCombinedImage = GERecon('SumOfSquares', channelImages);
        end            

        % Zero out kissoff views and apply gradwarp
        kissoffViews = header.RawHeader.kissoff_views;
        channelCombinedImage(:,1:kissoffViews) = 0;
        channelCombinedImage(:,(end-kissoffViews+1):end) = 0;
        gradwarpImage = GERecon('Gradwarp', abs(channelCombinedImage), corners, 'XRMW');

        if(header.RawHeader.hnover > 0)                
            gradwarpImage = gradwarpImage * (256 / (header.RawHeader.rc_xres * header.RawHeader.rc_yres));
        end

        % Rotate/Transpose
        rotatedTransposedSlice = GERecon('Orient', gradwarpImage, orientation);

        % Clip to range of shorts (match product functionality)
        rotatedTransposedSlice(rotatedTransposedSlice < 0) = 0;
        rotatedTransposedSlice(rotatedTransposedSlice > 32767) = 32767;
        finalImages(:,:,slice) = int16(rotatedTransposedSlice);

        imagesc(finalImages(:,:,slice));colormap(gray);colorbar;axis off;title(['Slice: ' num2str(slice) 'Phase: ' num2str(phase)]);
        drawnow;

        if(imageNumber < 10)
            imageNumberString = ['00' num2str(imageNumber)];
        elseif(imageNumber < 100)
            imageNumberString = ['0' num2str(imageNumber)];
        else
            imageNumberString = num2str(imageNumber);
        end

        matlabDicomPath = fullfile(pfilePath, 'matlabDicoms', filesep);

        GERecon('Dicom.Write', [matlabDicomPath 'Image_' imageNumberString '.dcm'], finalImages(:,:,slice), imageNumber, orientation, corners, (header.SeriesData.se_no * 100));

        imageNumber = imageNumber + 1;
    end
end

