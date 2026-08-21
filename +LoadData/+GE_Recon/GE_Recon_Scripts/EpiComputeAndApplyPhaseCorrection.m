function [ phaseCorrectionHandle ] = EpiComputeAndApplyPhaseCorrection( referencePfile )
%% EpiComputeAndApplyPhaseCorrection
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
%    phaseCorrectionHandle = EpiComputeAndApplyPhaseCorrection('Prrrrr.7');
%    Given the reference pfile 'Prrrrr.7', run reference scan
%    reconstruction processing to compute phase correction coefficients.
%    The computed coefficients are applied to the raw reference scan data
%    to demonstrate that applying a phase ramp in x,ky space with the
%    computed linear and constant coefficients shifts and lines up the
%    kx,ky reference scan data.
%    A handle to the computed coefficients is returned. The handle can be
%    used with the GERecon function to apply the coefficients to raw data.
%

    %% Load
    % Load a reference scan pfile and initialize an empty phase correction
    % reference handle. The phase correction reference handle will be filled
    % with coefficients and returned to the caller. The filled handle can be
    % used to apply phase correction using the computed coefficients.
    pfileHandle = GERecon('Pfile.Load', referencePfile);               
    phaseCorrectionHandle = GERecon('Epi.LoadReference');
    
    %% Account for reference views
    % Some multi-phase EPI scans acquire additional non-phase-encoded views
    % at the top or bottom of the EPI echo train. These views are not
    % needed for reference scan processing. This code determines the number
    % of additional reference views acquired with the reference scan data
    % so they can be discarded when the data is processed.    
    header = GERecon('Pfile.Header');              
    refViewsTop = header.RawHeader.extra_frames_top;
    refViewsBottom = header.RawHeader.extra_frames_bot;    
    
    %% Compute coefficients
    % Compute linear and constant phase correction coefficients from
    % reference scan data. As stated previously, if additional
    % non-phase-encoded reference views are acquired then discard those
    % views here.
    coefficientsFigureHandle = figure;
    for slice = 1:pfileHandle.slices
        for channel = 1:pfileHandle.channels
            % Echo dimension is unused for baseic Epi scans
            echo = 1;            
            kSpace = GERecon('Pfile.KSpace', slice, echo, channel);
            
            % Account for additional reference views that may be
            % acquired with multi-phase scans
            kSpaceTotalNumViews = size(kSpace,2);                
            kSpaceWithoutRefViews = kSpace(:,(refViewsTop+1):(kSpaceTotalNumViews-refViewsBottom));            
            
            % Compute phase correction coefficients, store in reference
            % handle, and plot coefficients
            computedCoefficients = GERecon('Epi.ComputeCoefficients', phaseCorrectionHandle, kSpaceWithoutRefViews, slice, channel);
            figure(coefficientsFigureHandle);
            subplot(2,1,1,'replace');
            plot(computedCoefficients(:,1));title(['Constant Coefficients - Slice: ' num2str(slice) ', Channel: ' num2str(channel)]);
            subplot(2,1,2,'replace');
            plot(computedCoefficients(:,2));title(['Linear Coefficients - Slice: ' num2str(slice) ', Channel: ' num2str(channel)]);            
            drawnow;
        end
        % If Phase alignment bit is set in the pool header 
        % phase alignment is performed on the computed coefficients
        GERecon('Epi.PhaseAlignment', phaseCorrectionHandle, slice)
    end
    
    %% Apply Coefficients
    % Apply the computed coefficients to the reference scan data. This step
    % is included to demonstrate that applying reference scan coefficients
    % aligns the EPI echo train data.
    % The GERecon utility transforms to x,ky space, applies a linear phase
    % ramp using the computed linear and constant coefficients, and then
    % transforms back to kx,ky.
    % The shifted reference scan data is plotted before and after the
    % GERecon utility shifts the echoes.
    
    applyCoefficientsFigureHandle = figure;
    for slice = 1:pfileHandle.slices
        for channel = 1:pfileHandle.channels
            % Echo dimension is unused for baseic Epi scans
            echo = 1;            
            kSpace = GERecon('Pfile.KSpace', slice, echo, channel);
            
            % Account for additional reference views that may be
            % acquired with multi-phase scans
            kSpaceTotalNumViews = size(kSpace,2);                
            kSpaceWithoutRefViews = kSpace(:,(refViewsTop+1):(kSpaceTotalNumViews-refViewsBottom));                     
            
            % Apply phase correction
            phaseCorrectedKSpace = GERecon('Epi.ApplyPhaseCorrection', phaseCorrectionHandle, kSpaceWithoutRefViews, slice, channel);

            % Plot data before and after shifiting data
            figure(applyCoefficientsFigureHandle);
            subplot(2,1,1,'replace');
            imagesc(abs(kSpaceWithoutRefViews)');colormap(gray);colorbar;axis off;title(['Raw Slice ' num2str(slice) ', Channel ' num2str(channel)]);
            subplot(2,1,2,'replace');
            imagesc(abs(phaseCorrectedKSpace)');colormap(gray);colorbar;axis off;title(['Corrected Slice: ' num2str(slice) ', Channel: ' num2str(channel)]);
            drawnow;
        end
    end    
    
end

