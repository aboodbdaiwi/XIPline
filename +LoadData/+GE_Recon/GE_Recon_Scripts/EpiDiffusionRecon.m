function [ finalImages ] = EpiDiffusionRecon(pfileFullPath)
%% EpiDiffusionRecon
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
%   Reconstruct Epi Diffusion data from the given pfile. For scans with
%   multiple pfiles, supply the first pfile and the subsequent pfiles will
%   automatically be opened.
%   This script will reconstruct each B-Value/Diffusion Direction acquired
%   in the Epi Diffusion scan. This script does not generate combined
%   diffusion images.
%   If this diffusion scan was acquired with an integrated reference scan,
%   then the reference data acquired during the first phase of the scan
%   will be used to compute phase correction coefficients. If this scan was
%   not acquired with an integrated reference scan then this script will
%   look (in the following order) for the following sources of reference
%   data: ref.h5 in pfile directory, ref.dat in pfile directory, reference
%   pfile in pfileDirectory/ref/Prrrrr.7.
%
%   finalImages = EpiDiffusionRecon(pfile)
%

    %% Determine secondary input locations
    % Find the vrgf.dat or vrgf_kernels.dat file containing the ramp
    % sampling interpolation kernels to use with this scan. This code will
    % look in the pfile directory for either of these two files.
    [pfilePath pfileName pfileExt] = fileparts(pfileFullPath);

    % Look for vrgf data
    vrgfDotDatPath = fullfile(pfilePath, 'vrgf.dat');
    vrgfKernelsDotDatPath = fullfile(pfilePath, 'vrgf_kernels.dat');
    if(exist(vrgfDotDatPath, 'file') == 2)
        % Found vrgf.dat, use for vrgf interpolation kernels
        vrgfInterpKernels = vrgfDotDatPath;
        vrgfEnabled = 1;
    elseif(exist(vrgfKernelsDotDatPath, 'file') == 2)
        % Look for vrgf_kernels.dat        
        vrgfInterpKernels = vrgfKernelsDotDatPath;    
        vrgfEnabled = 1;
    else
        % Could not find vrgf interpolation kernels!
        disp('Could not find vrgf interpolation kernels.');
        disp('vrgf interpolation is disabled.');
        vrgfEnabled = 0;
    end
    
    % Load pfile and determine if asset was enabled for this scan
    pfileHandle = GERecon('Pfile.Load', pfileFullPath);
    header = GERecon('Pfile.Header');

    % Look for AssetCalibration. Note that if you have a calibration pfile
    % you can generate the AssetCalibration.h5 file using this command:
    %
    %   GERecon('Calibration.Process', 'pathToCalPfile/Pxxxxx.7').
    %
    % The calibration HDF5 file will be saved in the calibration pfile's
    % directory.
    assetCalibrationPath = fullfile(pfilePath, 'AssetCalibration.h5');
    assetEnabled = (header.RawHeader.asset == 2 || header.RawHeader.asset == 7);
    if(exist(assetCalibrationPath, 'file') == 2 && assetEnabled)
        useAsset = 1;
        assetCalibration = assetCalibrationPath;
    else
        useAsset = 0;
    end

    %% Load or process reference data    
    % Integrated reference scan is indicated by rhref == 5. If this
    % condition is true then set the referenceData variable to the first
    % pfile in the scan. This pfile will be passed to the
    % EpiComputeAndApplyPhaseCorrection function. 
    % If this is not an integrated reference scan then determine where the
    % reference data can be obtained. The reference data may be supplied 
    % as a ref.dat file, a ref.h5 file, or a separate reference scan pfile.    
    
    % Determine where the reference data can be obtained
    if(header.RawHeader.ref == 5)
       referenceData = pfileFullPath;
       firstImageAcquisitionPhase = 2;
    else
       % Reference data must be from either a ref.dat/ref.h5 file or a 
       % separate reference scan.
       refDotDatPath = fullfile(pfilePath, 'ref.dat');
       refDotH5Path = fullfile(pfilePath, 'ref.h5');
       refPfilePath = fullfile(pfilePath, 'ref', [pfileName pfileExt]);
       if(exist(refDotDatPath, 'file') == 2)
           referenceData = refDotDatPath;
       elseif(exist(refDotH5Path, 'file') == 2)
           referenceData = refDotH5Path;
       elseif(exist(refPfilePath, 'file') == 2)
           referenceData = refPfilePath;
       end       
       firstImageAcquisitionPhase = 1;
    end
    
    % Load or process the reference data. For an integrated reference scan
    % or a scan with a separate reference scan pfile supplied, compute
    % phase correction coefficients. If ref.dat or ref.h5 is supplied, load
    % the coefficients from the given file.
    if(strcmp('.7', referenceData(end-1:end)))
        % Reference pfile was provided, load that and compute coefficients
        phaseCorrectionHandle = EpiComputeAndApplyPhaseCorrection(referenceData);    
    else       
        % ref.dat or ref.h5 was provided, load the reference data file
        phaseCorrectionHandle = GERecon('Epi.LoadReference', referenceData);        
    end
            
    %% Load VRGF interpolation kernels
    % VRGF (variable readout gradient filtering) refers to sampling data on
    % the gradient ramps during data acquisition. During recon, a sinc
    % interpolation is performed to interpolate non-linearly sampled
    % frequency data (i.e. data sampled on the gradient ramps) to linearly
    % sampled frequency data.
    % The sinc interpolation kernels that are used for each point in the 
    % interpolated output are contained in a vrgf.dat for Epi Diffusion
    % scans.
    if(vrgfEnabled == 1)
        vrgfHandle = GERecon('RampSampling.Load', vrgfInterpKernels);    
        kernelMatrix = GERecon('RampSampling.KernelMatrix', vrgfHandle);
        figure;
        subplot(2,1,1);
        imagesc(kernelMatrix);colormap(gray);colorbar;axis off;title('VRGF Interpolation Kernel Matrix');
        subplot(2,1,2);
        numRows = size(kernelMatrix, 1);
        plot(kernelMatrix(numRows/4, :));title(['Sinc Interpolation Kernel for Interpolated Point: ' num2str(numRows/4) ' of ' num2str(numRows)]);   
    end
            
    %% Initialize workspace and header variables
    % Pull values from the raw header and use these values to allocate
    % space for intermediate data.
    GERecon('Pfile.SetActive', pfileHandle);
    header = GERecon('Pfile.Header');
    imageSize = header.RawHeader.im_size;
    channelImages = single(zeros(imageSize, imageSize, pfileHandle.channels));    
    finalImages = int16(zeros(imageSize, imageSize, (pfileHandle.slices*pfileHandle.phases)));        
    
    %% Load Asset Calibration (if applicable)
    % If the user specified an Asset calibration file, load the file and
    % use Asset for channel combination.
    if(useAsset == 1)
        GERecon('Asset.LoadCalibration', assetCalibration);    
    end     
    
    % Initialize image number to zero and begin looping over all phases 
    % and slices. Each B-Value/Diffusion Direction combination is a
    % separate phase in the scan
    finalImageFigure = figure;
    imageNumber = 0;     
    for phase = firstImageAcquisitionPhase:pfileHandle.phases
       
        %% Determine T2 or BValue/Direction Indices
        % The following code is for product EPI diffusion sequences.
        % Product sequences acquire the diffusion data in the
        % following order (for X number of B-Values and Y number of
        % directions per B-Value):
        %
        % * Reference Phase (if integrated ref scan)
        % * T2 Phase(s)
        % * B0, Dir0
        % * B0, Dir1
        % * ...
        % * BX, Dir0
        % * BX, DirY
        %
        % The BValue and direction indices may be used to generate combined
        % images for each B-Value (combined image generation is not 
        % included in this example).
        numBValues = header.RawHeader.numbvals;
        numDiffusionDirections = header.RawHeader.numdifdirs;
        numT2Images = pfileHandle.phases - (numBValues * numDiffusionDirections) - (firstImageAcquisitionPhase-1);
        isT2Phase = phase < (numT2Images + firstImageAcquisitionPhase);
        if(isT2Phase)            
            currentNexCount = header.RawHeader.difnext2;
        else            
            currentBValueIndex = floor((phase - firstImageAcquisitionPhase - numT2Images) / numDiffusionDirections) + 1;         
            currentNexCount = header.RawHeader.difnextab(currentBValueIndex);
            % currentDirectionIndex = mod((phase - firstImageAcquisitionPhase - numT2Images), numDiffusionDirections) + 1;
        end
        
        % Prevent any subsequent errors related to indexing or dividing by zero
        if (currentNexCount == 0)
            currentNexCount = 1;
        end
        
        for slice = 1:pfileHandle.slices
            
            corners = GERecon('Pfile.Corners', slice);
            orientation = GERecon('Pfile.Orientation', slice);            
            
            % Initialize space for NEX'd magnitude image to zero.
            % Epi diffusion scans use the echo index as the NEX index.
            magnitudeImage = zeros(imageSize, imageSize);
            
            % If complex image nex'ing is enabled, accumulate nex'd data on
            % a channel by channel basis. The output will be an 
            % accumulated kSpace matrix for each channel. This accumulated
            % channel data is used by the remaining recon steps as if 
            % this were a single nex scan.
            complexNexEnabled = bitget(header.RawHeader.data_collect_type1, 26); 
            if(complexNexEnabled > 0)
                complexNexCount = currentNexCount;
                magnitudeNexCount = 1;
            else
                complexNexCount = 1;
                magnitudeNexCount = currentNexCount;
            end
            
            if(complexNexEnabled > 0)
                magnitudeNexCount = 1;
                complexNexChannelData = zeros(header.RawHeader.vrgfxres, header.RawHeader.da_yres - 1, pfileHandle.channels);
                
                for channel = 1:pfileHandle.channels
                    for echo = 1:complexNexCount
                        % Retrieve input kSpace                        
                        kSpace = GERecon('Pfile.KSpace', slice, echo, channel, phase);    
                        
                        % Apply phase correction
                        phaseCorrectedKSpace = GERecon('Epi.ApplyPhaseCorrection', phaseCorrectionHandle, kSpace, slice, channel);

                        if(vrgfEnabled == 1)
                            % Interpolate ramp sampled data and accumulate
                            % into complex nex combiner
                            interpolatedData = GERecon('RampSampling.Interpolate', vrgfHandle, phaseCorrectedKSpace);                        
                            GERecon('Epi.ComplexImageNex', interpolatedData);
                        else
                            % No ramp sampling interpolation to do, 
                            % accumulate into nex comnbiner
                            GERecon('Epi.ComplexImageNex', phaseCorrectedKSpace);
                        end
                    end
                    
                    complexNexChannelData(:,:,channel) = GERecon('Epi.ComplexImageNex');
                end
            end
            
            for echo=1:magnitudeNexCount
                for channel=1:pfileHandle.channels
                    
                    % If complex image nex is enabled then phase correction
                    % and ramp sampling interpolation was done above.
                    % Retrieve the dataToTransform from the complex nex
                    % channel data that was processed above.
                    if(complexNexEnabled)
                        dataToTransform = complexNexChannelData(:,:,channel);
                    else                                                
                        kSpace = GERecon('Pfile.KSpace', slice, echo, channel, phase);

                        % Apply phase correction
                        phaseCorrectedKSpace = GERecon('Epi.ApplyPhaseCorrection', phaseCorrectionHandle, kSpace, slice, channel);

                        % Interpolate ramp sampled data
                        if(vrgfEnabled == 1)
                            dataToTransform = GERecon('RampSampling.Interpolate', vrgfHandle, phaseCorrectedKSpace);                        
                        else
                            dataToTransform = phaseCorrectedKSpace;
                        end
                    end
                                        
                    image = GERecon('Transform', dataToTransform);
                    
                    % If ASSET and Homodyne are both enabled for this scan
                    % then ASSET is run on both the high pass filtered and
                    % low pass filtered images generated by the Homodyne
                    % algorithm. To enable this use case, the Transform
                    % command returns the high pass filtered and low pass
                    % filtered images in indices one and two of the third
                    % dimension of the return image. For this case, the
                    % channel images array must have space to store the
                    % additional high pass filtered / low pass filtered
                    % images. Resize the channelImages matrix here to
                    % enable this use case. Note that this resize happens
                    % only once, the first time through this loop.
                    if( (size(image,3) > 1) && (size(channelImages,4) == 1) )
                        channelImages = single(zeros(size(channelImages,1), size(channelImages,2), size(channelImages,3), 2));
                    end
                    
                    channelImages(:,:,channel,:) = image;
                end
                
                % ASSET Unalias
                if(useAsset)
                    channelCombinedImage = GERecon('Asset.Unalias', channelImages, corners);
                else                
                    channelCombinedImage = GERecon('SumOfSquares', channelImages);
                end     
                
                % Accumulate magnitude NEX's
                magnitudeImage = magnitudeImage + abs(channelCombinedImage);                
            end
            
            magnitudeImage = magnitudeImage / magnitudeNexCount;
            
            % Zero out kissoff views prior to realtime field adjustment and
            % gradwarp (matches product functionality).
            kissoffViews = header.RawHeader.kissoff_views;
            magnitudeImage(:,1:kissoffViews) = 0;
            magnitudeImage(:,(end-kissoffViews+1):end) = 0;
            
            %% Apply Realtime Field Adjustment, if enabled
            % Realtime field adjustment interpolates image data along the
            % phase encoding direction to compensate for eddy currents that
            % exist when the diffusion gradients are turned on. This
            % utility will use the slice and phase indices to determine 
            % how much compensation to apply for the current slice based 
            % on the slice location and diffusion gradient strength for 
            % the current phase.
            % Note that scans acquired without realtime field adjustment
            % enabled lack the information required in the raw header to
            % apply this correction. Thus, only scans acquired with
            % realtime field adjustment enabled can take advantage of this
            % utility
            if(header.RawHeader.hoecc > 0) 
                % Apply Realtime Field Adjustment
                fieldAdjusted = GERecon('Epi.RealtimeFieldAdjustment', magnitudeImage, slice, phase);
            else
                % Do not apply Realtime Field Adjustment
                fieldAdjusted = magnitudeImage;
            end            
            
            %% Finalization
            % Apply gradwarp, orient image, scale and convert to shorts
            % prior to creating a DICOM image for this slice
            gradwarpImage = GERecon('Gradwarp', fieldAdjusted, corners, 'XRMW');           
            
            % Partial ky homodyne scans have additional scaling applied at
            % the end of recon (match product scaling here)
            partialKyHomodyne = bitget(header.RawHeader.data_collect_type, 5);
            if(partialKyHomodyne)                
                gradwarpImage = gradwarpImage * (256 / (header.RawHeader.rc_xres * header.RawHeader.rc_yres));
            end

            % Orient Image
            rotatedTransposedSlice = GERecon('Orient', gradwarpImage, orientation);

            % Clip to range of shorts (match product functionality)            
            rotatedTransposedSlice(rotatedTransposedSlice < 0) = 0;
            rotatedTransposedSlice(rotatedTransposedSlice > 32767) = 32767;
            finalImages(:,:,imageNumber+1) = int16(rotatedTransposedSlice);

            figure(finalImageFigure);
            imagesc(finalImages(:,:,imageNumber+1));colormap(gray);colorbar;axis off;title(['Slice: ' num2str(slice) 'Phase: ' num2str(phase)]);
            drawnow;

            if(imageNumber < 10)
                imageNumberString = ['00' num2str(imageNumber)];
            elseif(imageNumber < 100)
                imageNumberString = ['0' num2str(imageNumber)];
            else
                imageNumberString = num2str(imageNumber);
            end

            matlabDicomPath = fullfile(pfilePath, 'matlabDicoms', filesep);
            
            % Diffusion Image Type Annotation
            % Possible Values:
            %   DiffusionRightLeftDicomValue = 3
            %   DiffusionAnteriorPosteriorDicomValue = 4
            %   DiffusionSuperiorInferiorDicomValue = 5
            %   DiffusionT2DicomValue = 14
            %   DiffusionCombinedDicomValue = 15
            %   DiffusionDtiDicomValue = 16
            %   DiffusionDirection1DicomValue = 43
            %   DiffusionDirection2DicomValue = 44
            %   DiffusionDirection3DicomValue = 45
            %   DiffusionDirection4DicomValue = 46
            % The diffusion image type is not in the pool header. Thus, 
            % the diffusion image type must come from an external source.
            % By default, the diffusion image type is set to Dir 1.
            diffusionImageTypeTag.Group = hex2dec('0043');
            diffusionImageTypeTag.Element = hex2dec('1030');
            diffusionImageTypeTag.VRType = 'SS';
            diffusionImageTypeTag.Value = 43;
            
            % 1-based geometric index of slice
            geometricIndexTag.Group = hex2dec('0020');
            geometricIndexTag.Element = hex2dec('9057');
            geometricIndexTag.VRType = 'UL';
            geometricIndexTag.Value = slice;
            
            % BValue Bias Factor
            % Product DICOM images have a bias factor of 1e9 added to the
            % bValue dicom field (first integer in 0043,1039) for scans
            % with more than one bValue. The bias factor, if applied, is
            % stored in DICOM field (0043,107f). The bias factor
            % functionality is replicated here.
            % The bValue is not present in the pfile header. Thus, the
            % bValue must come from another external source. By default,
            % the bValue is set to 0.
            bValue = 0;
            bValueTag.Group = hex2dec('0043');
            bValueTag.Element = hex2dec('1039');
            bValueTag.VRType = 'IS';   
            bValueTag.Value = [num2str(bValue) '\ 0 \ 0 \ 0'];            
                        
            if(numBValues > 1)
                % Update bValue tag to include bias factor and include
                % b-Value bias factor in dicom image header.
                bValueBiasFactor = 1000000000;
                bValue = bValue + bValueBiasFactor;                
                bValueTag.Value = [num2str(bValue) '\ 0 \ 0 \ 0'];
                
                bValueBiasFactorTag.Group = hex2dec('0043');
                bValueBiasFactorTag.Element = hex2dec('107f');
                bValueBiasFactorTag.VRType = 'IS';   
                bValueBiasFactorTag.Value = num2str(bValueBiasFactor);

                GERecon('Dicom.Write', [matlabDicomPath 'Image_' imageNumberString '.dcm'], finalImages(:,:,imageNumber+1), imageNumber, orientation, ...
                        corners, (header.SeriesData.se_no*100), header.SeriesData.se_desc, diffusionImageTypeTag, geometricIndexTag, bValueTag, bValueBiasFactorTag);                
            else
                GERecon('Dicom.Write', [matlabDicomPath 'Image_' imageNumberString '.dcm'], finalImages(:,:,imageNumber+1), imageNumber, orientation, ...
                        corners, (header.SeriesData.se_no*100), header.SeriesData.se_desc, diffusionImageTypeTag, geometricIndexTag, bValueTag);
            end            

            imageNumber = imageNumber + 1;            
        end            
    end
end

