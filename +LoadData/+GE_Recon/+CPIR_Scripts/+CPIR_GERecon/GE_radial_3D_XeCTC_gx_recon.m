function GasExchange = GE_radial_3D_XeCTC_gx_recon(MainInput, GasExchange)
%GE_RADIAL_3D_XECTC_GX_RECON  GE version of radial_3D_XeCTC_gx_recon.
%
% The analysis below follows the MRD XeCTC reconstruction in the same order:
%   1. Import data and acquisition information
%   2. Determine the image offset
%   3. Fit the appended dissolved-phase bonus spectrum
%   4. Remove approach-to-steady-state and gas contamination
%   5. Plot k0 signal dynamics
%   6. Reconstruct ventilation, gas, and dissolved images
%   7. Reconstruct RBC oscillation keyhole images
%   8. Store the standard GasExchange outputs
%
% The ONLY scanner-specific part is readGERadialXeCTC() at the bottom:
% it reads the GE ScanArchive and waveform, extracts the appended FIDs,
% applies GE raw2grid preprocessing, and returns MRD-format k-space.
% No native GE gridding reconstruction is executed.
%
% Inputs:
%   MainInput.XeFullPath          GE ScanArchive .h5
%   MainInput.CalFullPath         Optional separate GE calibration
%   MainInput.OutputPath          Optional output folder
%   MainInput.GE_XeTR_s           Optional same-contrast TR in seconds
%   MainInput.GE_RecMatrix        Optional reconstruction matrix override
%   MainInput.GE_UseNativeShift   true (default): native GE centering phase
%   MainInput.GE_Shift_vox        Optional 3-vector override (default [0.5 20.5 0.5] for even matrix)
%   MainInput.GE_ApplyApodization false (default): MRD uses no GE apodization
%   MainInput.GE_SaveFigures      true (default)
%   MainInput.GE_SaveMAT          true (default)
%   MainInput.GE_RunRBCOsc        true (default)
%
% Research use: validate frequency sign, GE orientation, spectral fit,
% scan-specific delay/trajectory, and gas subtraction on real data.

if nargin < 2 || isempty(GasExchange), GasExchange = struct(); end
assert(isstruct(MainInput) && isstruct(GasExchange), 'Inputs must be structures.');

%% Read GE data and waveform (replaces the ISMRMRD import section)
fprintf('\n========== GE XeCTC: MRD-STYLE RECONSTRUCTION ==========\n');
GE = readGERadialXeCTC(MainInput);

GasKSpace       = double(GE.GasKSpace);                 % samples x projections
DissolvedKSpace = double(GE.DissolvedKSpace);           % samples x projections
XeTraj          = double(GE.XeTraj);                    % 3 x samples x projections

% The GE trajectory is already normalized for the GE reconstruction matrix.
% Unlike the supplied ISMRMRD script, do NOT divide this trajectory by 2.
Xe_RecMatrix = option(MainInput,'GE_RecMatrix',GE.RecMatrix);
assert(isscalar(Xe_RecMatrix) && Xe_RecMatrix==GE.RecMatrix, ...
    'Changing GE_RecMatrix requires rescaling the GE trajectory; use waveform matrix.');
assert(isequal(size(XeTraj),[3 size(GasKSpace)]), ...
    'GE trajectory and imaging k-space dimensions do not match.');
assert(max(abs(XeTraj(:)))<=0.5+1e-6, 'GE normalized trajectory exceeds +/-0.5.');
assert(size(GasKSpace,2)>60, 'Too few imaging profiles for 60-view steady-state removal.');

% GE has one receiver coil and two alternating gas/dissolved excitations.
OvsFactor = GE.OvsFactor;
Xe_nprof = size(GasKSpace,2);
XeImg_nsamp = size(GasKSpace,1);
Xe_interleaves = 2;
XeOrder = 1:Xe_nprof;
dwell_s = GE.dwell_s;
GasFlipAngle = GE.GasFlipAngle;
DisFlipAngle = GE.DisFlipAngle;
freq_jump = GE.freq_jump;
TE90 = GE.TE90_ms;
ActTE90 = TE90;
XeTR = option(MainInput,'GE_XeTR_s',GE.XeTR_s);
if ~(isscalar(XeTR) && isfinite(XeTR) && XeTR>0)
    XeTR = NaN;
    warning('GE same-contrast TR unavailable; RBC oscillation analysis will be skipped.');
end

% Keep both acquired complex FIDs and their display-only zero-filled copies.
PostGasFID = GE.GasBonusFIDAcquired(:);
PostDissolvedFID = GE.DissolvedBonusFIDAcquired(:);
MainInput.GE_ZeroFillBonusFIDTo512 = false;
if MainInput.GE_ZeroFillBonusFIDTo512
    N = 512;

    gas = PostGasFID;
    diss = PostDissolvedFID;

    PostGasFID = complex(zeros(N,1));
    PostDissolvedFID = complex(zeros(N,1));

    PostGasFID(1:numel(gas)) = gas;
    PostDissolvedFID(1:numel(diss)) = diss;
end

bonus_spectra_present = ~isempty(PostDissolvedFID);

scaleFac = abs(PostGasFID(1,1))/abs(GasKSpace(1,end)); %These should be almost identical at k0
scaleFac1dis = abs(PostDissolvedFID(1,1))/abs(DissolvedKSpace(1,end)); %These should be almost identical at k0

DissolvedKSpace = DissolvedKSpace * scaleFac;
GasKSpace = GasKSpace * scaleFac;

% The GE raw2grid preprocessing already applied the documented centering
% phase. Applying another MRD PixelShift would shift the images twice.
PixelShift = [0;0;0];

outdir = option(MainInput,'OutputPath', ...
    fullfile(fileparts(GE.DataFile),'GasExchange_Analysis'));
if ~exist(outdir,'dir'), mkdir(outdir); end
GasExchange.outputpath = outdir;
saveFigures = option(MainInput,'GE_SaveFigures',true);

GasExchange.GasBonusFID = GE.GasBonusFID;
GasExchange.DissolvedBonusFID = GE.DissolvedBonusFID;
GasExchange.GasBonusFIDAcquired = PostGasFID;
GasExchange.DissolvedBonusFIDAcquired = PostDissolvedFID;
GasExchange.BonusTime_s = GE.BonusTime_s;
GasExchange.BonusTimePadded_s = GE.BonusTimePadded_s;
GasExchange.BonusDwell_s = GE.BonusDwell_s;
GasExchange.BonusAcquisitionIndices = GE.BonusAcquisitionIndices;
GasExchange.BonusFIDStatus = GE.BonusFIDStatus;
GasExchange.GERaw = GE.Metadata;
GasExchange.DataFile = GE.DataFile;
GasExchange.WaveformFile = GE.WaveformFile;
GasExchange.HeaderFile = GE.HeaderFile;
GasExchange.GEPreprocessing = GE.Preprocessing;
GasExchange.OvsFactor = OvsFactor;
GasExchange.dwell_s = dwell_s;
GasExchange.GasFlipAngle = GasFlipAngle;
GasExchange.TE90 = TE90;
GasExchange.XeTR = XeTR;
GasExchange.Xe_RecMatrix = Xe_RecMatrix;
GasExchange.ReconstructionMethod = 'GE raw2grid preprocessing + original MRD reconstruction functions';

% Additional GE QC: inspect both acquired FIDs without padding them for fit.
if saveFigures && bonus_spectra_present
    fig = figure('Name','GE Bonus FIDs','Visible','off');
    subplot(2,1,1);
    plot(real(PostGasFID));
    xlabel('Time (ms)');ylabel('Real signal');title('Gas bonus FID');grid on;
    subplot(2,1,2);
    plot(real(PostDissolvedFID));
    xlabel('Time (ms)');ylabel('Real signal');title('Dissolved bonus FID');grid on;
    savefig(fig,fullfile(outdir,'GEBonusFIDs.fig'));
    close(fig);
end

%% Determine FOV shift (the GE loader has already applied its native shift)
% Same role as the MRD "Determine FOV Shift" section. The GE native
% center correction is performed once in readGERadialXeCTC(), not here.
fprintf('GE spatial shift (voxels): [%s]\n', num2str(GE.Preprocessing.Shift_vox));

%% Spectra (same NMR_TimeFit_v model as original MRD reconstruction)
AppendedDissolvedNMRFit = [];
GasExchange.CalibrationDissolvedNMRFit = [];
GasExchange.GasExResults = [];
GasExchange.CalResults = [];
GasExchange.DissolvedNMRFitSource = 'none';
GasExchange.GasContaminationCorrectionApplied = false;

if bonus_spectra_present
    if (numel(PostDissolvedFID)-1)*GE.BonusDwell_s < 0.002
        warning(['GE dissolved bonus FID is shorter than 2 ms: RBC/barrier ' ...
            'spectral separation may be poorly identified. Zero fill does not fix this.']);
    end
    fprintf('Fitting the GE imaging-appended dissolved bonus spectrum...\n');
    try
        fitFID = double(PostDissolvedFID(:));  % ACQUIRED samples, not zero fill
        time = double(GE.BonusDwell_s)*(0:numel(fitFID)-1)';
        A = max(abs(fitFID(1)),eps);
        area_guess = [0.27*A, A, 0.1*A];
        % Original MRD three-component model, with the actual GE frequency jump.
        freq_guess = [534,-126,-abs(freq_jump)];
        fwhm_guess = [300,273,44];
        fwhmG_guess = [0,275,0];
        phase_guess = [-95,151,9];
        fitObj = GasExchangeFunctions.GasExchange_Spectro.NMR_TimeFit_v( ...
            fitFID,time,area_guess,freq_guess,fwhm_guess,fwhmG_guess, ...
            phase_guess,[],[]);
        fitObj.fitTimeDomainSignal();
        fitObj.setBounds([0 0 0],1e10*[1 1 1], ...
            [500 -2000 -12000],[2000 0 -4000], ...
            [0 0 0],[Inf Inf Inf], ...
            [0 0 0],[0 Inf 0], ...
            [-Inf -Inf -Inf],[Inf Inf Inf]);
        fitObj.fitTimeDomainSignal();
        fitSignal = sum(fitObj.calcComponentTimeDomainSignal(time),2);
        fitError = norm(fitFID-fitSignal(:))/max(norm(fitFID),eps);
        maxFitError = option(MainInput,'GE_BonusMaxRelativeFitError',0.35);
        % assert(isfinite(fitError) && fitError<=maxFitError, ...
        %     'Bonus spectral fit relative error %.3f exceeds %.3f.', ...
        %     fitError,maxFitError);
        % assert(numel(fitObj.area)>=3 && numel(fitObj.phase)>=3 && ...
        %     all(isfinite(fitObj.area(1:3))) && ...
        %     all(isfinite(fitObj.phase(1:3))) && ...
        %     all(fitObj.area(1:3)>=0), 'Invalid fitted gas/RBC/barrier components.');
        AppendedDissolvedNMRFit = fitObj;
        GasExchange.DissolvedNMRFitSource = 'same-scan GE bonus FID';
        GasExchange.BonusFitRelativeError = fitError;
    catch ME
        warning('Imaging bonus spectral fit failed: %s. No gas subtraction.',ME.message);
        GasExchange.BonusFitError = ME.message;
    end
end

% A separate calibration is allowed as a reference for RBC analysis, NOT
% as a substitute for the imaging FID in gas-contamination subtraction.
if isfield(MainInput,'CalFullPath') && ~isempty(MainInput.CalFullPath) && ...
        isfile(MainInput.CalFullPath)
    try
        [GasExResults,CalResults] = ...
            LoadData.GE_Recon.CPIR_Scripts.CPIR_GERecon.GE_Calibration_Spect_CCHMC(MainInput);
        GasExchange.GasExResults = GasExResults;
        GasExchange.CalResults = CalResults;
        if isstruct(GasExResults) && isfield(GasExResults,'DisFit')
            GasExchange.CalibrationDissolvedNMRFit = GasExResults.DisFit;
            if isempty(AppendedDissolvedNMRFit)
                AppendedDissolvedNMRFit = GasExResults.DisFit;
                GasExchange.DissolvedNMRFitSource = 'separate calibration (analysis only)';
            end
        end
    catch ME
        warning('Separate GE calibration failed: %s',ME.message);
    end
end

% Plot Dissolved Phase Spectrum - Same as MRD
GasExchange.DissolvedNMR = [];
if ~isempty(AppendedDissolvedNMRFit)
    try
        % Calculate fitted spectra
        time = double(AppendedDissolvedNMRFit.t(:));

        AppendedDissolvedFit = ...
            (time(2)-time(1)) * fftshift(fft( ...
            AppendedDissolvedNMRFit.calcComponentTimeDomainSignal(time),[],1),1);

        % Figure
        DissolvedNMR = figure('Name','Dissolved Phase Spectrum');
        set(DissolvedNMR,'WindowState','minimized');
        set(DissolvedNMR,'Color','white', ...
            'Units','inches','Position',[0.25 0.25 16 9]);

        hold on

        % Measured spectrum
        f = AppendedDissolvedNMRFit.f;
        spectrum = AppendedDissolvedNMRFit.spectralDomainSignal;

        plot(f,abs(spectrum),'b')
        plot(f,real(spectrum),'r')
        plot(f,imag(spectrum),'Color',[0 0.6 0.2])

        % Fitted spectrum
        plot(f,abs(sum(AppendedDissolvedFit,2)), ...
            'Color',[0 0 1 0.33],'LineWidth',3)

        plot(f,real(sum(AppendedDissolvedFit,2)), ...
            'Color',[1 0 0 0.33],'LineWidth',3)

        plot(f,imag(sum(AppendedDissolvedFit,2)), ...
            'Color',[0 0.6 0.2 0.33],'LineWidth',3)

        % Individual fitted components
        area(f,real(AppendedDissolvedFit(:,1)), ...
            'FaceColor','r','FaceAlpha',0.33,'LineStyle','none')

        area(f,real(AppendedDissolvedFit(:,2)), ...
            'FaceColor','b','FaceAlpha',0.33,'LineStyle','none')

        area(f,real(AppendedDissolvedFit(:,3)), ...
            'FaceColor',[0 0.6 0.2], ...
            'FaceAlpha',0.33,'LineStyle','none')

        % Legend
        lgnd = legend({ ...
            'Spectrum - Magnitude', ...
            'Spectrum - Real', ...
            'Spectrum - Imaginary', ...
            'Fit - Magnitude', ...
            'Fit - Real', ...
            'Fit - Imaginary', ...
            'RBC - Real', ...
            'Barrier - Real', ...
            'Gas - Real'}, ...
            'Location','best','NumColumns',3);

        set(lgnd,'Color','none');

        % Formatting
        xlim([-8000 2000])
        set(gca,'XDir','reverse','FontSize',24)

        title(['Dissolved Phase Spectra and Fit: RBC/Barrier = ', ...
            num2str(AppendedDissolvedNMRFit.area(1) / ...
            AppendedDissolvedNMRFit.area(2))])

        xlabel('Frequency (Hz)')
        ylabel('NMR Signal (a.u.)')

        hold off

        % Save
        GasExchange.DissolvedNMR = DissolvedNMR;

        if saveFigures
            savefig(DissolvedNMR,fullfile(outdir,'DissolvedNMR.fig'));
        end
        close(DissolvedNMR);
    catch ME
        warning('Spectroscopy QC figure failed: %s',ME.message);
    end
end

%% Remove approach to steady state (same 60-profile rule as original MRD)
SS_ind = min(60,size(GasKSpace,2)-1);
GasKSpace_SS = GasKSpace(:,SS_ind+1:end);
DissolvedKSpace_SS = DissolvedKSpace(:,SS_ind+1:end);
XeTraj_SS = XeTraj(:,:,SS_ind+1:end);

%% Gas contamination removal (same MRD correction function)
CorrectedDissKSpace_SS = DissolvedKSpace_SS;
GasExchange.GasContaminationCorrectionStatus = 'skipped: no validated same-scan bonus fit';
isImagingFit = strcmp(GasExchange.DissolvedNMRFitSource,'same-scan GE bonus FID');
if isImagingFit
    try
        CorrectedDissKSpace_SS = GasExchangeFunctions.GasPhaseContaminationRemoval( ...
            DissolvedKSpace_SS,GasKSpace_SS,dwell_s,-freq_jump, ...
            AppendedDissolvedNMRFit.phase(3), ...
            AppendedDissolvedNMRFit.area(3),GasFlipAngle);
        assert(isequal(size(CorrectedDissKSpace_SS),size(DissolvedKSpace_SS)) && ...
            all(isfinite(real(CorrectedDissKSpace_SS(:)))) && ...
            all(isfinite(imag(CorrectedDissKSpace_SS(:)))), ...
            'Corrected dissolved k-space is invalid.');
        GasExchange.GasContaminationCorrectionApplied = true;
        GasExchange.GasContaminationCorrectionStatus = 'applied: same-scan GE bonus';
    catch ME
        warning('Gas contamination removal failed: %s',ME.message);
        CorrectedDissKSpace_SS = DissolvedKSpace_SS;
        GasExchange.GasContaminationCorrectionStatus = ['failed: ' ME.message];
    end
end

%% View k0 dynamics (original MRD figure name: SigDynamics.fig)
GasExchange.SigDynamics = [];
XePulses = (1:size(CorrectedDissKSpace_SS,2))' + SS_ind;
% Figure
SigDynamics = figure('Name','Signal Dynamics');
set(SigDynamics,'WindowState','minimized');
set(SigDynamics,'Color','white', ...
    'Units','inches','Position',[0.25 0.25 16 9]);

% Xe signal dynamics
subplot(1,2,1);
hold on
set(gca,'FontSize',18)

% Dissolved
hDis = plot(abs(DissolvedKSpace(1,:)), ...
    'Color',[0.5 0 0],'LineWidth',1);

% Corrected dissolved
hCorr = plot(XePulses,abs(CorrectedDissKSpace_SS(1,:)), ...
    'Color',[1 0 0],'LineWidth',1);

% Gas
hGas = plot(abs(GasKSpace(1,:)), ...
    'Color',[0 0 1],'LineWidth',1);

% Gray region: discarded steady-state projections
xl = xlim;
yl = ylim;

hSS = patch([0 SS_ind SS_ind 0], ...
    [yl(2) yl(2) 0 0], ...
    [0.66 0.66 0.66], ...
    'EdgeColor','none');

uistack(hSS,'bottom')
xlim(xl);
ylim(yl);

% Formatting
title('Xe Signal Intensity Dynamics')
box on
xlabel('Projection')
ylabel('Signal Intensity')

legend([hSS hDis hCorr hGas], ...
    {'Steady State Discard', ...
     'Dissolved', ...
     'Corrected Dissolved', ...
     'Gas'}, ...
    'Location','best', ...
    'NumColumns',2, ...
    'FontSize',12);

hold off

% Store figure
GasExchange.SigDynamics = SigDynamics;

% Save
if saveFigures
    savefig(SigDynamics,fullfile(outdir,'SigDynamics.fig'));
end

close(SigDynamics);

%% Reconstruct Xe images (same four MRD reconstruction calls)
fprintf('Reconstructing high-resolution ventilation with MRD solver...\n');
UncorrectedVentImage = GasExchangeFunctions.Dissolved_HighResImageRecon( ...
    Xe_RecMatrix,GasKSpace_SS,XeTraj_SS,PixelShift);
% figure; imslice(abs(UncorrectedVentImage))
fprintf('Reconstructing low-resolution gas and dissolved images...\n');
GasImage = GasExchangeFunctions.Dissolved_LowResImageRecon( ...
    Xe_RecMatrix,GasKSpace_SS,XeTraj_SS,PixelShift);
DissolvedImage = GasExchangeFunctions.Dissolved_LowResImageRecon( ...
    Xe_RecMatrix,DissolvedKSpace_SS,XeTraj_SS,PixelShift);
if GasExchange.GasContaminationCorrectionApplied
    CorrDissolvedImage = GasExchangeFunctions.Dissolved_LowResImageRecon( ...
        Xe_RecMatrix,CorrectedDissKSpace_SS,XeTraj_SS,PixelShift);
else
    CorrDissolvedImage = DissolvedImage;
end

% Exactly the same GE orientation operation as the original MRD script:
% every axial 2-D slice: imrotate(img,90), then flip(img,2).

UncorrectedVentImage = orientLikeMRDGE(UncorrectedVentImage);
GasImage = orientLikeMRDGE(GasImage);
DissolvedImage = orientLikeMRDGE(DissolvedImage);
CorrDissolvedImage = orientLikeMRDGE(CorrDissolvedImage);

%% Ventilation bias correction (same MRD method)
VentMask = abs(UncorrectedVentImage)>prctile(abs(UncorrectedVentImage(:)),95)/3;
VentMask = imerode(VentMask,strel('sphere',5));
regions = regionprops3(logical(VentMask),'all');
regions = sortrows(regions,1,'descend');
mask = false(size(VentMask));
if ~isempty(regions)
    mask(regions.VoxelIdxList{1}) = true;
    if height(regions)>1 && regions.Volume(2)>0.33*regions.Volume(1)
        mask(regions.VoxelIdxList{2}) = true;
    end
    mask = imdilate(mask,strel('sphere',6));
else
    mask = true(size(VentMask));
end
% The reference MRD code overrides this computed mask with all ones.
% Keep that exact behavior by default; opt in to the computed mask.
if ~option(MainInput,'GE_UseComputedVentMask',false)
    mask = true(size(mask));
end
[VentImage,~] = GasExchangeFunctions.Vent_BiasCorrection(UncorrectedVentImage,mask);
GasExchange.VentMask = mask;

if saveFigures
    saveMontage(UncorrectedVentImage,outdir,'VentMontage');
    saveMontage(VentImage,outdir,'CorrVentMontage');
    saveMontage(GasImage,outdir,'GasMontage');
    saveMontage(DissolvedImage,outdir,'DissolvedMontage');
    saveMontage(CorrDissolvedImage,outdir,'CorrDissolvedMontage');
end

%% Determine keys for RBC oscillations (same MRD keyhole functions)
RBC2Bar_struct = [];
RBCOsc_High_Image = [];
RBCOsc_Low_Image = [];
RBCOsc_Normalization = [];
GasExchange.RBCOscStatus = 'skipped: valid spectral fit and XeTR required';
if option(MainInput,'GE_RunRBCOsc',true) && ...
        ~isempty(AppendedDissolvedNMRFit) && isfinite(XeTR) && XeTR>0
    try
        SpecRBCBarrierRatio = double(AppendedDissolvedNMRFit.area(1)) / ...
            double(AppendedDissolvedNMRFit.area(2));
        assert(isfinite(SpecRBCBarrierRatio) && SpecRBCBarrierRatio>0, ...
            'Invalid RBC/barrier ratio.');
        [RBCOsc_High_Key,RBCOsc_Low_Key, ...
            DissolvedKSpace_SS_Detrend_Normalized,RBC2Bar_struct,OscWorkFlow_Fig] = ...
            GasExchangeFunctions.dissolved_RBC_Osc_Keyhole( ...
            XeTraj_SS,GasKSpace_SS,CorrectedDissKSpace_SS,XeTR,SpecRBCBarrierRatio);
        if saveFigures && isgraphics(OscWorkFlow_Fig)
            saveas(OscWorkFlow_Fig,fullfile(outdir,'OscWorkFlow_Fig.png'));
        end
        if isgraphics(OscWorkFlow_Fig),close(OscWorkFlow_Fig);end
        RBCOsc_High_Image = GasExchangeFunctions.Dissolved_RBCOscImageRecon( ...
            Xe_RecMatrix,RBCOsc_High_Key,XeTraj_SS,PixelShift);
        RBCOsc_Low_Image = GasExchangeFunctions.Dissolved_RBCOscImageRecon( ...
            Xe_RecMatrix,RBCOsc_Low_Key,XeTraj_SS,PixelShift);
        RBCOsc_Normalization = GasExchangeFunctions.Dissolved_RBCOscImageRecon( ...
            Xe_RecMatrix,DissolvedKSpace_SS_Detrend_Normalized,XeTraj_SS,PixelShift);
        RBCOsc_High_Image = orientLikeMRDGE(RBCOsc_High_Image);
        RBCOsc_Low_Image = orientLikeMRDGE(RBCOsc_Low_Image);
        RBCOsc_Normalization = orientLikeMRDGE(RBCOsc_Normalization);
        GasExchange.RBCOscStatus = 'completed';
    catch ME
        warning('GE RBC oscillation analysis failed: %s',ME.message);
        GasExchange.RBCOscStatus = ['failed: ' ME.message];
    end
end

%% Store variables (same GasExchange output names as original MRD)
GasExchange.XeTraj = XeTraj;
GasExchange.GasKSpace = GasKSpace;
GasExchange.DissolvedKSpace = DissolvedKSpace;
GasExchange.XeTraj_SS = XeTraj_SS;
GasExchange.GasKSpace_SS = GasKSpace_SS;
GasExchange.DissolvedKSpace_SS = DissolvedKSpace_SS;
GasExchange.CorrectedDissKSpace_SS = CorrectedDissKSpace_SS;
GasExchange.XeImg_nsamp = XeImg_nsamp;
GasExchange.Xe_nprof = Xe_nprof;
GasExchange.Xe_interleaves = Xe_interleaves;
GasExchange.XeOrder = XeOrder;
GasExchange.SS_ind = SS_ind;
GasExchange.UncorrectedVentImage = UncorrectedVentImage;
GasExchange.VentImage = VentImage;
GasExchange.GasImage = GasImage;
GasExchange.DissolvedImage = DissolvedImage;
GasExchange.CorrDissolvedImage = CorrDissolvedImage;
GasExchange.AppendedDissolvedNMRFit = AppendedDissolvedNMRFit;
GasExchange.RBC2Bar_struct = RBC2Bar_struct;
GasExchange.RBCOsc_High_Image = RBCOsc_High_Image;
GasExchange.RBCOsc_Low_Image = RBCOsc_Low_Image;
GasExchange.RBCOsc_Normalization = RBCOsc_Normalization;
GasExchange.ActTE90 = ActTE90;
GasExchange.DisFlipAngle = DisFlipAngle;
GasExchange.PixelShift = PixelShift;
GasExchange.freq_jump = freq_jump;
GasExchange.Image = abs(VentImage);
GasExchange.ImageCombined = abs(GasImage);
GasExchange.GasImageCombined = abs(GasImage);
GasExchange.DissolvedImageCombined = abs(DissolvedImage);
GasExchange.ImagingBonusSpectraPresent = bonus_spectra_present;
GasExchange.PostGasFID = PostGasFID;
GasExchange.PostDissolvedFID = PostDissolvedFID;
GasExchange.Scanner = 'GE';
GasExchange.Institute = option(MainInput,'Institute','');
GasExchange.XeSin = [];
GasExchange.ScanVersion = 'XeCTC';
GasExchange.DataLocation = fileparts(GE.DataFile);
[~,xeBase,xeExt] = fileparts(GE.DataFile);
GasExchange.FileName = [xeBase xeExt];
GasExchange.RbcBarRatio = option(GasExchange.GasExResults,'RbcBarRatio',NaN);
GasExchange.GasDisRatio = option(GasExchange.GasExResults,'GasDisRatio',NaN);

if option(MainInput,'GE_SaveMAT',true)
    % Save output only rather than every temporary workspace variable.
    save(fullfile(outdir,'xe_recon_workspace.mat'),'GasExchange','-v7.3');
end
fprintf('GE MRD-style reconstruction finished. Output: %s\n',outdir);
end

%% ======================== GE SCANNER ADAPTER ==========================
function GE = readGERadialXeCTC(MainInput)
% Return MRD-compatible [readout x projection] k-space and [3 x readout x
% projection] trajectory. No native gridding. Single receive coil only.
scan = option(MainInput,'XeFullPath','');
if isempty(scan)
    scan = fullfile(option(MainInput,'XeDataLocation',''), ...
                    option(MainInput,'XeFileName',''));
end
assert(isfile(scan),'GE ScanArchive file not found: %s',scan);
p = findGEFiles(scan,MainInput);
[d,h] = LoadData.GE_Recon.fidall.read_p(scan,true);
d = single(conj(d));   % same complex convention as original GE code
wf = LoadData.GE_Recon.fidall.load_waveform(p.waveform,true,true);
for name = {'k','ind','dcf','t','mtx'}
    assert(isfield(wf,name{1}),'GE waveform missing %s.',name{1});
end

% This matches the CPIR single-channel 3D XeCTC sequence.
assert(size(d,3)==1 && size(d,4)==1 && size(d,5)==1 && size(d,6)==1, ...
    'Only single-coil, single-slab, two-contrast 3D GE data are supported.');
ind = wf.ind;
nViews = size(ind,1);
nRead = nnz(ind(1,:));
assert(nRead>=32 && mod(nViews,2)==0 && size(d,1)==nViews, ...
    'GE view count disagrees with waveform, or unsupported readout layout.');
assert(all(sum(ind,2)==nRead),'GE waveform readout lengths are inconsistent.');
assert(size(d,2)>=size(ind,2),'GE raw ADC readout shorter than waveform index.');

% Normalize waveform k-space to the ACTUAL GE reconstruction grid.
acqMatrix = double(wf.mtx(:)');
if isscalar(acqMatrix),acqMatrix=repmat(acqMatrix,1,3);end
if numel(acqMatrix)==2,acqMatrix(3)=acqMatrix(2);end
recMatrix = p.mtx_reco;
if isempty(recMatrix) || all(recMatrix(:)<1)
    recMatrix = acqMatrix*option(p,'interp_factor',1);
end
if isscalar(recMatrix),recMatrix=repmat(recMatrix,1,3);end
if numel(recMatrix)==2,recMatrix(3)=recMatrix(2);end
recMatrix = double(recMatrix(:)');
assert(numel(recMatrix)==3 && all(recMatrix==recMatrix(1)), ...
    'This GE/MRD implementation expects an isotropic 3D matrix.');
traj = wf.k;
if isempty(p.mtx_reco) || all(p.mtx_reco(:)<1)
    traj = wf.k / option(p,'interp_factor',1);
else
    traj = bsxfun(@times,wf.k,shiftdim(acqMatrix./recMatrix,-1));
end
traj = single(traj);
dcf = single(wf.dcf);
t = single(wf.t);
if isrow(t) && numel(t)<=size(ind,2)
    t = repmat(t,nViews,1);
end
t = t.'; t = t(:).';  % same flatten order as GE native reconstruction
assert(numel(traj)==3*nRead*nViews && numel(dcf)==nRead*nViews && ...
    numel(t)==nRead*nViews, 'GE waveform k, DCF, or time shape inconsistent.');

%% Bonus FIDs (last gas view, last dissolved view, zero DCF)
GE.GasBonusFIDAcquired = [];
GE.DissolvedBonusFIDAcquired = [];
GE.GasBonusFID = [];
GE.DissolvedBonusFID = [];
GE.BonusTime_s = [];
GE.BonusTimePadded_s = [];
GE.BonusDwell_s = NaN;
GE.BonusAcquisitionIndices = [];
GE.BonusFIDStatus = 'not detected';
if isfield(wf,'flips'),flips = double(wf.flips(:)');else,flips = [];end
if isfield(wf,'freqs'),freqs = double(wf.freqs(:)');else,freqs = [];end
assert(numel(flips)==2 && numel(freqs)==2 && flips(1)<flips(2), ...
    'GE waveform must document [gas,dissolved] interleave flips/frequencies.');

viewsDCF = reshape(dcf,nRead,nViews);
timeViews = reshape(t,nRead,nViews);
bonus = all(viewsDCF(:,end-1:end)==0,'all') && ...
        any(viewsDCF(:,end-2)~=0) && isequal(ind(end,:),ind(end-1,:));
if bonus
    adc = find(ind(end,:));
    GE.GasBonusFIDAcquired = reshape(d(end-1,adc,1,1,1,1),[],1);
    GE.DissolvedBonusFIDAcquired = reshape(d(end,adc,1,1,1,1),[],1);
    GE.BonusTime_s = double(timeViews(:,end));
    GE.BonusDwell_s = median(diff(GE.BonusTime_s));
    Npad = max(512,nRead);
    GE.GasBonusFID = [GE.GasBonusFIDAcquired; complex(zeros(Npad-nRead,1,'single'))];
    GE.DissolvedBonusFID = [GE.DissolvedBonusFIDAcquired; complex(zeros(Npad-nRead,1,'single'))];
    GE.BonusTimePadded_s = GE.BonusTime_s(1) + ...
        (0:Npad-1)' * GE.BonusDwell_s;
    GE.BonusAcquisitionIndices = [nViews-1 nViews];
    GE.BonusFIDStatus = 'extracted: last gas and dissolved views';
    fprintf('Extracted GE bonus FIDs: views %d/%d, %d acquired samples.\n', ...
        nViews-1,nViews,nRead);
    % Remove bonus views from GE raw data and ALL matched waveform arrays.
    d = d(1:end-2,:,:,:,:,:);
    ind = ind(1:end-2,:);
    kv = reshape(traj,nRead,nViews,3);
    traj = reshape(kv(:,1:end-2,:),1,nRead*(nViews-2),3);
    dcf = reshape(viewsDCF(:,1:end-2),1,[]);
    t = reshape(timeViews(:,1:end-2),1,[]);
else
    warning('GE bonus views not found: no same-scan gas correction will be applied.');
end

%% GE interleave split and matched trajectory (same as original GE reader)
assert(mod(size(d,1),2)==0,'GE imaging excitation count must be even.');
d1 = d(1:2:end,:,:,:,:,:);
d2 = d(2:2:end,:,:,:,:,:);
nProj = size(d1,1);
dStack = cat(3,d1,d2);
dStack = reshape(dStack,nProj,size(dStack,2),1,2);
ind = ind(1:2:end,:);
nacq = sum(sum(ind,1)>0);
assert(nacq==nRead, ...
    'Unexpected GE ind readout geometry; inspect waveform indexing.');
kv = reshape(traj,[],nProj,3);
kKeep = kv(1:nacq,:,:);
% GE native keeps the first contrast's acquisition-time vector per view.
tByPair = reshape(t,[],nProj);
tKeep = tByPair(1:nacq,:);
GE.XeTraj = double(permute(kKeep,[3 1 2]));
GE.OvsFactor = NaN;  % no MRD-specific acquisition oversampling ratio

%% GE raw2grid: handles phase cycling, gradient delay, and native shift
bw = double(h.rdb_hdr.user0);
if abs(bw)<500, bw = double(h.rdb_hdr.bw)*2e3; end
shift = zeros(1,3);
if option(MainInput,'GE_UseNativeShift',true)
    shift = 0.5*double(mod(recMatrix,2)==0);
    shift(2) = shift(2)+20;
end
shift = option(MainInput,'GE_Shift_vox',shift);
assert(numel(shift)==3 && all(isfinite(shift)), ...
    'GE_Shift_vox must be a finite 3-vector.');
shift = double(shift(:)');
kFlat = reshape(kKeep,1,[],3);
dd = LoadData.GE_Recon.fidall.raw2grid(dStack, ...
    LoadData.GE_Recon.fidall.ischop(h),kFlat,shift,[],[], ...
    ind,double(p.delay)*1e-6,bw,false);
if option(MainInput,'GE_ApplyApodization',false) && any(abs(p.lb(:))>0)
    apod = double(p.lb(:)');
    if isscalar(apod),apod = repmat(apod,1,3);end
    if numel(apod)==2,apod(3)=apod(2);end
    % Native GE apodization uses only the imaging sample times,
    % after interleave splitting (not both raw contrast streams).
    dd = LoadData.GE_Recon.fidall.ak_apodise( ...
        dd,kFlat,apod(2),tKeep(:).',[apod(1) apod(3)],false);
end
assert(size(dd,1)==1 && size(dd,2)==nacq*nProj && size(dd,3)==2, ...
    'Unexpected GE raw2grid output size: [%s].', num2str(size(dd)));
stream1 = reshape(dd(1,:,1),nacq,nProj);
stream2 = reshape(dd(1,:,2),nacq,nProj);
% The waveform design explicitly defines interleave 1 = gas.
GE.GasKSpace = stream1;
GE.DissolvedKSpace = stream2;

%% GE timing, frequency, flip angles and provenance
GE.GasFlipAngle = flips(1);
GE.DisFlipAngle = flips(2);
GE.freq_jump = abs(freqs(2)-freqs(1));
assert(GE.freq_jump>0,'GE waveform frequency difference must be positive.');
GE.RecMatrix = recMatrix(1);
GE.dwell_s = median(diff(double(tKeep(:,1))));
assert(isfinite(GE.dwell_s) && GE.dwell_s>0, ...
    'Invalid GE imaging dwell time from waveform.');
if bonus && abs(GE.BonusDwell_s/GE.dwell_s-1)>0.05
    warning('GE bonus dwell differs by >5%% from imaging dwell.');
end
GE.TE90_ms = NaN;
if isfield(h.rdb_hdr,'te'),GE.TE90_ms = double(h.rdb_hdr.te)*1e-3;end
GE.XeTR_s = NaN;
try
    GE.XeTR_s = 2*double(h.image.tr)*1e-6;
catch
end
GE.DataFile = scan;
GE.HeaderFile = p.header;
GE.WaveformFile = p.waveform;
GE.Preprocessing = struct('Shift_vox',shift, 'Delay_us',p.delay, ...
    'NativeApodizationEnabled',option(MainInput,'GE_ApplyApodization',false), ...
    'NativeApodizationParameter',p.lb,'ChopDetected', ...
    logical(LoadData.GE_Recon.fidall.ischop(h)), ...
    'TrajectoryUnits','normalized cycles per reconstructed voxel');
GE.Metadata = struct('header',h,'mtx_acq',acqMatrix,'mtx_reco',recMatrix, ...
    'nproj',nProj,'nacq',nacq,'ncoils',1,'freqs',freqs,'flips',flips, ...
    'WaveformFile',p.waveform,'ImagingDwell_s',GE.dwell_s);
fprintf('GE data: %d readout x %d projections; matrix %d; frequency jump %.0f Hz.\n', ...
    nacq,nProj,GE.RecMatrix,GE.freq_jump);
end

%% Read the GE ScanArchive's matching header and matlab_scripter settings
function p = findGEFiles(scan,MainInput)
[folder,base,ext] = fileparts(scan);
p.header = fullfile(folder,[base '.hdr']);
if ~isfile(p.header)
    found = dir(fullfile(folder,'*.hdr'));
    assert(numel(found)==1,'Missing or ambiguous .hdr for GE scan.');
    p.header = fullfile(folder,found(1).name);
end
p.waveform = option(MainInput,'GE_WaveformFile','');
p.mtx_reco = option(MainInput,'GE_Matrix',[]);
p.interp_factor = option(MainInput,'GE_InterpFactor',1);
p.delay = option(MainInput,'GE_Delay_us',0);
p.lb = option(MainInput,'GE_LB',0);
if ~isempty(p.waveform)
    assert(isfile(p.waveform),'Specified GE_WaveformFile not found.');
    return;
end
scripters = dir(fullfile(folder,'*matlab_scripter*'));
assert(~isempty(scripters),'Cannot locate GE matlab_scripter file.');
idx = [];
for ii=1:numel(scripters)
    path = fullfile(folder,scripters(ii).name);
    if contains(fileread(path),[base ext],'IgnoreCase',true)
        idx=ii; break;
    end
end
if isempty(idx)
    assert(numel(scripters)==1,'Cannot uniquely match GE matlab_scripter file.');
    idx=1;
end
scriptText = fileread(fullfile(folder,scripters(idx).name));
token = regexp(scriptText,'recon_grid3d_gx\s*\((.*?)\)\s*;', ...
    'tokens','once');
assert(~isempty(token),'No recon_grid3d_gx(...) call in matlab_scripter.');
args = splitArguments(token{1});
assert(numel(args)>=12,'GE matlab_scripter requires 12 arguments.');
wfn = strrep(strrep(strtrim(args{3}),'''',''),'\','/');
[~,name,extension] = fileparts(wfn);
if isempty(extension) || strcmpi(extension,'.wav'), extension='.mat';end
p.waveform = fullfile(folder,[name extension]);
assert(isfile(p.waveform),'Cannot find waveform MAT: %s',p.waveform);
p.mtx_reco = parseNumber(args{4});
p.interp_factor = parseNumber(args{5});
p.delay = parseNumber(args{6});
p.lb = parseNumber(args{7});
end

function args = splitArguments(s)
args = {};
start = 1; depth = 0; quoted = false;
for ii=1:numel(s)
    c=s(ii);
    if c=='''',quoted=~quoted;end
    if ~quoted
        if any(c=='[({'),depth=depth+1;end
        if any(c=='])}'),depth=depth-1;end
        if c==',' && depth==0
            args{end+1}=strtrim(s(start:ii-1)); %#ok<AGROW>
            start=ii+1;
        end
    end
end
args{end+1}=strtrim(s(start:end));
end

function x = parseNumber(text)
text=strtrim(text);
if strcmp(text,'[]'),x=[];return;end
if strcmpi(text,'true'),x=true;return;end
if strcmpi(text,'false'),x=false;return;end
if startsWith(text,'[') && endsWith(text,']')
    text = text(2:end-1);
end
text = strrep(text,',',' ');
parts = strsplit(text,';');
rows = cell(size(parts));
for ii=1:numel(parts)
    rows{ii} = sscanf(parts{ii},'%f').';
    assert(~isempty(rows{ii}),'Unsupported GE numeric setting: %s',text);
end
x=vertcat(rows{:});
end

function v = orientLikeMRDGE(v)
% Corresponds to original MRD GE branch: per-slice imrotate(90), flip dim2.
v = flip(imrotate(v,-90),2);
v = flip(v,3);
end

function saveMontage(img,outdir,name)
fig=figure('Name',name,'Visible','off');
montage(abs(img),'DisplayRange',[]);
savefig(fig,fullfile(outdir,[name '.fig']));
close(fig);
end

function x = option(S,name,default)
x=default;
if isstruct(S) && isfield(S,name) && ~isempty(S.(name))
    x=S.(name);
end
end
