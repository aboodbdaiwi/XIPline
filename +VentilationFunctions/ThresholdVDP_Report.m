
function ThresholdVDP_Report(Ventilation, Proton, MainInput)
    % ------------------------------------------------------------------------
    %  Ventilation Analysis Report - 
    % -------------------------------------------------------------------------
    clc;
    
    % DEFINE INPUTS
    pptDir         = Ventilation.outputpath;
    todayStr       = datestr(now,'yyyymmdd');
    pptxFileName = ['ThresholdVDP_Report_' todayStr];
    pptxName       = fullfile(pptDir, [pptxFileName,'.pptx']);
    
    % Reference Values 
    refStr = { '-', '-', ...
        Ventilation.HealthyRef.CoV,...
        Ventilation.HealthyRef.VHI,...
        Ventilation.HealthyRef.TH.VDPULN,...
        Ventilation.HealthyRef.TH.LVV, ...
        Ventilation.HealthyRef.TH.HVV, ...
        Ventilation.HealthyRef.skewness, ...
        Ventilation.HealthyRef.kurtosis,... 
        Ventilation.HealthyRef.DDI2D, ...
        Ventilation.HealthyRef.DDI3D
        };
    
    % colored rows
    rowColors = [
        0.9 0.9 0.9;    % SNRlung
        0.9 0.9 0.9;  % SNRvv
        0.8 0.8 1.0;    % CoV
        0.8 0.8 1.0;    % VHI
        0.9 0.4 0.4;    % VDP
        0.9 0.8 0.3;    % LVV
        0.3 0.6 0.9;    % HVV
        .8 .7 .6;    % Skew
        .8 .7 .6;    % TLV
        .9 0.4 0.7     % DDI .9 0.4 0.7
        .9 0.4 0.7     % DDI .9 0.4 0.7
        .7 0.7 0.7     % header
    ];
    
    % Open/Create PPTX
    isOpen = Global.exportToPPTX();
    if ~isempty(isOpen), Global.exportToPPTX('close'); end
    if isfile(pptxName)
        Global.exportToPPTX('open', pptxName);
        Global.exportToPPTX('switchslide', 1);
    else
        Global.exportToPPTX('new','Dimensions',[16 9],'Title','Ventilation Analysis','Author','CPIR @ CCHMC');
    end
    
    % Add Slide
    Global.exportToPPTX('addslide');
    
    % Add Title Bar
    Global.exportToPPTX('addtext','Ventilation Analysis Report', 'Position',[1.5 0 5 0.5], 'FontSize',25,'FontWeight','bold','Color',[1 0 0],'BackgroundColor',[1 1 1],'HorizontalAlignment','center');
    
    SubjInfo = {MainInput.SubjectID,MainInput.Age,MainInput.Sex,MainInput.Disease, MainInput.ScanDate};
    % Header row 
    SubjInfosummary = cell(2, 5);
    SubjInfosummary(1,:) = {
        {'Subject ID','BackgroundColor',rowColors(end,:),'FontWeight','bold'},
        {'Age(y)','BackgroundColor',rowColors(end,:),'FontWeight','bold'},
        {'Sex','BackgroundColor',rowColors(end,:),'FontWeight','bold'},
        {'Disease','BackgroundColor',rowColors(end,:),'FontWeight','bold'},   
        {'Scan Date','BackgroundColor',rowColors(end,:),'FontWeight','bold'}  
    };
    % Fill rows
    for i = 1:size(SubjInfosummary,2)
        SubjInfosummary{2,i} = {SubjInfo{i}, 'BackgroundColor', rowColors(2,:), 'FontWeight', 'bold'};
    end
    % Add to slide
    Global.exportToPPTX('addtable', SubjInfosummary, 'Position', [0.1 0.5 8.2 0.2], ...
        'Vert', 'middle', 'Horiz', 'center', 'FontSize', 16);
    
    
    
    settinglabels = {'Scanner','Scan Software','Sequence','Recon','XIPline Commit','Denoise','N4Bias','Method','AgeCor.','ProcessDate', 'Analyst Initials'};
    if Ventilation.N4Analysis == 0
        N4bias = 'no';
    else
        N4bias = 'yes';
    end
    settings = {MainInput.Scanner,MainInput.ScannerSoftware,MainInput.SequenceType,MainInput.Recon,MainInput.AnalysisCode_hash,MainInput.denoiseXe,N4bias,'Threshold',Ventilation.HealthyRef.AgeCorrected,todayStr, MainInput.Analyst};
    % Header row
    settingsssummary = cell(11, 2);
    settingsssummary(1,:) = {
        {'Setting','BackgroundColor',rowColors(end,:),'FontWeight','bold'},
        {'Value','BackgroundColor',rowColors(end,:),'FontWeight','bold'}
    };
    % Fill rows
    for i = 1:numel(settinglabels)
        settingsssummary{i+1,1} = {settinglabels{i}, 'BackgroundColor', rowColors(2,:), 'FontWeight', 'bold'};
        settingsssummary{i+1,2} = settings{i};
    end
    % Add to slide
    Global.exportToPPTX('addtable', settingsssummary, 'Position', [0.1 1.4 3.1 5], ...
        'Vert', 'middle', 'Horiz', 'center', 'FontSize', 15);
    
    if isempty(Ventilation.overallMeanCV)
        meanCV = 'NaN';
    else
        meanCV = num2str(Ventilation.overallMeanCV, '%.2f');
    end
    
    if isempty(Ventilation.overallVHI)
        vhi = 'NaN';
    else
        vhi = num2str(Ventilation.overallVHI, '%.2f');
    end   
    % DDI 2D Calculation
    if strcmp(Ventilation.DDI2D, 'yes')
        if strcmp(Ventilation.DDIDefectMap, 'Threshold')
            DDI2D = sprintf('%.2f ± %.2f', Ventilation.DDI2D_mean, Ventilation.DDI2D_std);
        else
            DDI2D = 'NaN';
        end
    else
        DDI2D = 'NaN';
    end
    
    % DDI 3D Calculation
    if strcmp(Ventilation.DDI3D, 'yes')
        if strcmp(Ventilation.DDIDefectMap, 'Threshold')
            DDI3D = sprintf('%.2f ± %.2f', Ventilation.DDI3D_mean, Ventilation.DDI3D_std);
        else
            DDI3D = 'NaN';
        end
    else
        DDI3D = 'NaN';
    end

    % Values and references
    valStr = { Ventilation.SNR_lung, Ventilation.SNR_vv, meanCV, vhi, Ventilation.Threshold.VDP, abs(Ventilation.Threshold.THBins(2)-Ventilation.Threshold.THBins(1)), Ventilation.Threshold.THBins(4), Ventilation.skewness, Ventilation.kurtosis, DDI2D, DDI3D };
    
    labels = {'SNRlung','SNRvv','CoV','VHI (%)','VDP (%)','LVV (%)','HVV (%)','skewness','kurtosis','DDI(2D)','DDI(3D)'};
    
    
    % Header row
    resultssummary = cell(10, 3);
    resultssummary(1,:) = {
        {'Metric','BackgroundColor',rowColors(end,:),'FontWeight','bold'},
        {'Value','BackgroundColor',rowColors(end,:),'FontWeight','bold'},
        {'Reference','BackgroundColor',rowColors(end,:),'FontWeight','bold'}
    };
    
    % Fill rows
    for i = 1:numel(labels)
        resultssummary{i+1,1} = {labels{i}, 'BackgroundColor', rowColors(i,:), 'FontWeight', 'bold'};
        resultssummary{i+1,2} = valStr{i};
        resultssummary{i+1,3} = refStr{i};
    end
    
    % Add to slide
    Global.exportToPPTX('addtable', resultssummary, 'Position', [3.3 1.4 5 4], ...
        'Vert', 'middle', 'Horiz', 'center', 'FontSize', 15);
    
    % Insert Histogram Image
    histImg = fullfile(Ventilation.outputpath,'histogram.png');
    if isfile(histImg)
        Global.exportToPPTX('addpicture',histImg,'Position',[2 6.7 6 2]);
    end
    
    % 8. CREATE MONTAGES WITH 2 ROWS
    close all;
    x0 = 8.5; w = 7.3; yStep = 2.1; y0 = 0.1;
    figXe  = makeTwoRowMontage(Ventilation.Image, '129Xe Ventilation');
    % ----------- Proton montage preparation with auto-detect color or grayscale
    Pr = Proton.ProtonRegistered; %Proton.ProtonRegisteredColored;  % fallback
    if ndims(Pr) == 4  && size(Pr,4) == 3% grayscale: [X Y Z 1]
        Pr = permute(Pr, [1 2 4 3]);  % 
        Pr = squeeze(Pr);  % remove singleton dimension
    end
    figPr = makeTwoRowMontage(Pr, 'Proton');
    
    figVDP = makeTwoRowMontage(Ventilation.Threshold.VentDefectmap, 'VDP Map');
    
    MVent = Ventilation.Mask_Vent_Reg;  % fallback
    if ndims(MVent) == 4  && size(MVent,4) == 3% grayscale: [X Y Z 1]
        MVent = permute(MVent, [1 2 4 3]);  % 
        MVent = squeeze(MVent);  % remove singleton dimension
    end
    figMask_Vent = makeTwoRowMontage(MVent, 'MVent');
    
    XePos  = get(figXe, 'Position');
    PrPos  = get(figPr, 'Position');
    DefPos = get(figVDP, 'Position');
    Mask_VentPos = get(figMask_Vent, 'Position');
    
    foldername = "VDP_Analysis\";
    Venthist = imread([pptDir,'\', char(foldername) ,'TH_Histogram.png']);    
    Global.exportToPPTX('addpicture',Venthist,'Position',[3.3 5.7 5.5 3]);
    Global.exportToPPTX('addtext','Histogram', 'Position',[3.3 5.7 2 0.1], 'FontSize',14,'FontWeight','bold','Color',[1 0 0],'HorizontalAlignment','left');
    
    % 9. ADD IMAGES TO SLIDE
    Global.exportToPPTX('addpicture', figXe,  'Position', [x0 y0         w w*XePos(4)/XePos(3)]);
    Global.exportToPPTX('addpicture', figPr,  'Position', [x0 y0+yStep   w w*PrPos(4)/PrPos(3)]);
    Global.exportToPPTX('addpicture', figMask_Vent, 'Position', [x0 y0+2*yStep w w*Mask_VentPos(4)/Mask_VentPos(3)]);
    Global.exportToPPTX('addpicture', figVDP, 'Position', [x0 y0+3*yStep w w*DefPos(4)/DefPos(3)]);
    
    
    Global.exportToPPTX('addtext','129Xe', 'Position',[x0 y0 2 0.1], 'FontSize',14,'FontWeight','bold','Color',[1 0 0],'HorizontalAlignment','left');
    Global.exportToPPTX('addtext','Anatomical', 'Position',[x0 y0+yStep 2 0.1], 'FontSize',14,'FontWeight','bold','Color',[1 0 0],'HorizontalAlignment','left');
    Global.exportToPPTX('addtext','Mask', 'Position',[x0 y0+2*yStep 2 0.1], 'FontSize',14,'FontWeight','bold','Color',[1 0 0],'HorizontalAlignment','left');
    Global.exportToPPTX('addtext','Defects', 'Position',[x0 y0+3*yStep 2 0.1], 'FontSize',14,'FontWeight','bold','Color',[1 0 0],'HorizontalAlignment','left');
    
    Global.exportToPPTX('addtext',['Report location:',pptDir], 'Position',[0 8.6 16 0.3], 'FontSize',12,'FontWeight','bold','Color',[0 0 1],'HorizontalAlignment','left');
    
    % Save & close
    Global.exportToPPTX('save',pptxName);
    Global.exportToPPTX('close');
    fprintf('Report saved to %s\n',pptxName);
    close all;
    
    % save report as a PDF
    ppt = actxserver('PowerPoint.Application');
    
    presentation = ppt.Presentations.Open(pptxName);
    PDFoutputPath = fullfile(pptDir,[pptxFileName,'.pdf']);
    if exist(PDFoutputPath, 'file')
        delete(PDFoutputPath);  % remove existing PDF to avoid overwrite conflict
    end    
    presentation.SaveAs(PDFoutputPath, 32);
    presentation.Close();
    ppt.Quit();
    delete(ppt);
end

% Local function
function figHandle = makeTwoRowMontage(vol, figName)
    % Detect grayscale or RGB
    sz = size(vol);
    
    if numel(sz) == 3
        % Grayscale: X × Y × Z
        isRGB = false;
        numSlices = sz(3);
    elseif numel(sz) == 4 && sz(3) == 3
        % RGB: X × Y × Z × 3 → handled as indexed RGB
        isRGB = true;
        numSlices = sz(4);
    else
        error('Unsupported volume dimensions: must be grayscale or RGB 4D');
    end

    % Pad if odd number of slices
    if mod(numSlices, 2) == 1
        if isRGB
            vol2 = permute(vol, [1 2 4 3]); 
            padSlice = zeros(sz(1), sz(2), 1, 3);
            vol2 = cat(3, vol2, padSlice);  % add slice at end
            vol = permute(vol2, [1 2 4 3]); 
            numSlices = numSlices + 1;
        else
            padSlice = zeros(sz(1), sz(2), 1, 'like', vol);
            vol = cat(3, vol, padSlice);
            numSlices = numSlices + 1;
        end
    end

    % Create figure
    figHandle = figure('Name', figName, 'Visible', 'off');

    % Display montage
    if isRGB
        montage(vol, 'Size', [2, numSlices/2]);
    else
        montage(reshape(vol, [sz(1), sz(2), 1, numSlices]), ...
                'Size', [2, numSlices/2], 'DisplayRange', []);
    end

    % Format figure size
    set(gca, 'Units', 'pixels');
    axPos = get(gca, 'Position');
    set(figHandle, 'Units', 'pixels', 'Position', [100 100 axPos(3) axPos(4)]);
    set(gca, 'Units', 'normalized', 'Position', [0 0 1 1]);
end
