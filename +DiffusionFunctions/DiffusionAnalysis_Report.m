
function DiffusionAnalysis_Report(Diffusion, MainInput)
    % ------------------------------------------------------------------------
    %  Diffusion Analysis Report - 
    % -------------------------------------------------------------------------
    clc;
    
    % DEFINE INPUTS
    pptDir         = Diffusion.outputpath;
    todayStr       = datestr(now,'yyyymmdd');
    if strcmp(MainInput.CCHMC_DbDiffAnalysis, 'yes')
        pptxFileName = 'DiffusionAnalysis_Report';
        pptxName       = fullfile(pptDir, [pptxFileName,'.pptx']);
        pptoutputPattern = fullfile(pptDir, 'DiffusionAnalysis_Report.pptx');
    else
        pptxFileName = ['DiffusionAnalysis_Report_' todayStr];
        pptxName       = fullfile(pptDir, [pptxFileName,'.pptx']);
        pptoutputPattern = fullfile(pptDir, 'DiffusionAnalysis_Report_*.pptx');
    end
    pptFiles = dir(pptoutputPattern);   
    for k = 1:length(pptFiles)
        delete(fullfile(pptDir, pptFiles(k).name));  % delete each matching file
    end

    if strcmp(MainInput.CCHMC_DbDiffAnalysis, 'yes')
        pptoutputPattern = fullfile(pptDir, 'DiffusionAnalysis_Report.pdf');
    else
        pptoutputPattern = fullfile(pptDir, 'DiffusionAnalysis_Report_*.pdf');
    end

    pptFiles = dir(pptoutputPattern);   
    for k = 1:length(pptFiles)
        delete(fullfile(pptDir, pptFiles(k).name));  % delete each matching file
    end

    % Reference Values 
    refStr = { '-', '-', ...
        Diffusion.HealthyRef.ADC,...
        Diffusion.HealthyRef.LDR,...
        Diffusion.HealthyRef.NDR,...
        Diffusion.HealthyRef.HDR, ...
        Diffusion.HealthyRef.ADCskewness, ...
        Diffusion.HealthyRef.ADCkurtosis,... 
        Diffusion.HealthyRef.ADCFitR2,... 
        Diffusion.HealthyRef.DDC, ...
        Diffusion.HealthyRef.LmD, ...
        Diffusion.HealthyRef.alpha,...
        Diffusion.HealthyRef.R,...
        Diffusion.HealthyRef.r,...
        Diffusion.HealthyRef.h, ...
        Diffusion.HealthyRef.Lm, ...
        Diffusion.HealthyRef.SVR,... 
        Diffusion.HealthyRef.Na, ...
        Diffusion.HealthyRef.LmD        
        };
    
    % colored rows
    rowColors = [
        0.90 0.90 0.90;  % Nb-values (s/cm^2)
        0.90 0.90 0.90;  % SNR [b0,b5]
        0.70 0.70 0.88;  % ADC (cm^2/s)
        0.42 0.44 0.84;  % LDR (%)
        0.50 0.79 0.43;  % NDR (%)
        0.84 0.50 0.47;  % HDR (%)
        0.72 0.64 0.55;  % ADC skewness
        0.72 0.64 0.55;  % ADC kurtosis
        0.72 0.64 0.55;  % ADC Fit R2
        0.83 0.37 0.65;  % SEM-DDC (cm^2/s)
        0.83 0.37 0.65;  % SEM-LmD (µm)
        0.83 0.37 0.65;  % SEM-alpha
        0.54 0.75 0.81;  % CM-R (µm)
        0.54 0.75 0.81;  % CM-r (µm)
        0.54 0.75 0.81;  % CM-h (µm)
        0.54 0.75 0.81;  % CM-Lm (µm)
        0.54 0.75 0.81;  % CM-SVR (cm^-1)
        0.54 0.75 0.81;  % CM-Na (mm^-3)
        0.70 0.70 0.70     % header
    ];

    
    % Open/Create PPTX
    isOpen = Global.exportToPPTX();
    if ~isempty(isOpen), Global.exportToPPTX('close'); end
    if isfile(pptxName)
        Global.exportToPPTX('open', pptxName);
        Global.exportToPPTX('switchslide', 1);
    else
        Global.exportToPPTX('new','Dimensions',[16 9],'Title','Diffusion Analysis','Author','CPIR @ CCHMC');
    end
    
    % Add Slide
    Global.exportToPPTX('addslide');
    
    % Add Title Bar
    Global.exportToPPTX('addtext','Diffusion Analysis Report', 'Position',[1.5 0 5 0.5], 'FontSize',25,'FontWeight','bold','Color',[1 0 0],'BackgroundColor',[1 1 1],'HorizontalAlignment','center');
    
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
    % addtable
    Global.exportToPPTX('addtable', SubjInfosummary, 'Position', [0.1 0.5 8.2 0.2], ...
        'Vert', 'middle', 'Horiz', 'center', 'FontSize', 16);
    
    
    settinglabels = {'Scanner','Scan Software','Sequence','Recon','XIPline Commit','Denoise','N4Bias',...
        'ADC Fit', 'SEM-Analysis','CM-Analysis', 'RefAgeCor.','Image Quality','Note','ProcessDate', 'Analyst Initials'};
    if MainInput.N4Bias == 0
        N4bias = 'no';
    else
        N4bias = 'yes';
    end
    settings = {MainInput.Scanner,MainInput.ScannerSoftware,MainInput.SequenceType,MainInput.Recon,MainInput.AnalysisCode_hash,MainInput.denoiseXe,N4bias,...
        Diffusion.ADCFittingType,Diffusion.SEMMorphometry, Diffusion.CMMorphometry, Diffusion.HealthyRef.AgeCorrected,MainInput.ImageQuality, MainInput.Note, todayStr, MainInput.Analyst};
    % Header row
    settingsssummary = cell(13, 2);
    settingsssummary(1,:) = {
        {'Setting','BackgroundColor',rowColors(end,:),'FontWeight','bold'},
        {'Value','BackgroundColor',rowColors(end,:),'FontWeight','bold'}
    };
    % Fill rows
    for i = 1:numel(settinglabels)
        settingsssummary{i+1,1} = {settinglabels{i}, 'BackgroundColor', rowColors(2,:), 'FontWeight', 'bold'};
        settingsssummary{i+1,2} = settings{i};
    end
    % addtable
    Global.exportToPPTX('addtable', settingsssummary, 'Position', [0.1 1.4 3.1 5], ...
        'Vert', 'middle', 'Horiz', 'center', 'FontSize', 15);
    

    % Values and references
    valStr = { '','',...
        [num2str(round(Diffusion.meanADC,4)),'±',num2str(round(Diffusion.stdADC,4))],...
        Diffusion.LDR, Diffusion.NDR, Diffusion.HDR, ...
        Diffusion.ADC_skewness,Diffusion.ADC_kurtosis,...
        Diffusion.meanRmap, ... 
        [num2str(Diffusion.DDC_mean),'±',num2str(Diffusion.DDC_std)],...
        [num2str(Diffusion.LmD_mean),'±',num2str(Diffusion.LmD_std)],...
        [num2str(Diffusion.alpha_mean),'±',num2str(Diffusion.alpha_std)],...
        [num2str(Diffusion.R_mean),'±',num2str(Diffusion.R_std)],...
        [num2str(Diffusion.r_mean),'±',num2str(Diffusion.r_std)],...
        [num2str(Diffusion.h_mean),'±',num2str(Diffusion.h_std)],...
        [num2str(Diffusion.Lm_mean),'±',num2str(Diffusion.Lm_std)],...
        [num2str(Diffusion.SVR_mean),'±',num2str(Diffusion.SVR_std)],...
        [num2str(Diffusion.Na_mean),'±',num2str(Diffusion.Na_std)],...
        };
    
    labels = {
        'b-value (s/cm²)', ...
        'SNR', ...
        'ADC (cm²/s)', ...
        'LDR (%)', ...
        'NDR (%)', ...
        'HDR (%)', ...
        'ADC skewness', ...
        'ADC kurtosis', ...
        'ADC Fit R²', ...
        'SEM-DDC(cm²/s)', ...
        'SEM-LmD (µm)', ...
        'SEM-alpha', ...
        'CM-R (µm)', ...
        'CM-r (µm)', ...
        'CM-h (µm)', ...
        'CM-Lm (µm)', ...
        'CM-SVR(cm⁻¹)', ...
        'CM-Na (mm⁻³)'
    };

    
    % Header row
    resultssummary = cell(size(labels,2), 3);
    resultssummary(1,:) = {
        {'Metric','BackgroundColor',rowColors(end,:),'FontWeight','bold'},
        {'Value','BackgroundColor',rowColors(end,:),'FontWeight','bold'},
        {'Reference','BackgroundColor',rowColors(end,:),'FontWeight','bold'}
    };
    
    % Colors
    abnormalColor = [1 0.8 0.8]; % light red
    normalColor = [1 1 1];       % white
    
    % Fill table rows
    for i = 1:numel(labels)
        label = labels{i};
        if isnumeric(valStr{i}) || islogical(valStr{i})
            val = round(valStr{i}, 2);
        else
            val = valStr{i};  % Leave unchanged if not numeric
        end
        ref = refStr{i};
    
        % Default metric cell style
        resultssummary{i+1,1} = {label, 'BackgroundColor', rowColors(i,:), 'FontWeight', 'bold'};
    
        % Default color for value cell
        bgColor = normalColor;
    
        % Check if numeric and if reference is in 'mean±std' format
        if ischar(ref) && contains(ref, '±') && isnumeric(val) && ~isnan(val)
            try
                % Extract mean and std using regular expression
                tokens = regexp(ref, '([\d.]+)\s*±\s*([\d.]+)', 'tokens');
                if ~isempty(tokens)
                    refMean = str2double(tokens{1}{1});
                    refSD   = str2double(tokens{1}{2});
                    lowerBound = refMean - refSD;
                    upperBound = refMean + refSD;
                    % Flag abnormal value
                    if i == 3 || i == 4 || i == 8 || i == 9
                        if val < lowerBound || val > upperBound
                            bgColor = abnormalColor;
                        end
                    else
                        if val > upperBound
                            bgColor = abnormalColor;
                        end                        
                    end
                end
            catch
                warning('Could not parse reference value for %s: %s', label, ref);
            end
        elseif ischar(ref) && contains(ref, '≤') && isnumeric(val) && ~isnan(val)
            % Handle one-sided reference (e.g., '≤4')
            try
                maxVal = str2double(erase(ref, '≤'));
                if val > maxVal
                    bgColor = abnormalColor;
                end
            catch
                warning('Could not parse ≤ reference for %s: %s', label, ref);
            end
        end
    
        % Write value and reference to result table
        resultssummary{i+1,2} = {val, 'BackgroundColor', bgColor};
        resultssummary{i+1,3} = ref;
    end

    % Insert Histogram Image
    try
        histImg = fullfile(Diffusion.outputpath,'LB_ADC_Histogram.png'); %LB_ADC_Histogram
    catch
        histImg = fullfile(Diffusion.outputpath,'ADC_Histogram.png');
    end
    if isfile(histImg)
        Global.exportToPPTX('addpicture',histImg,'Position',[8 6 8 3]);
    end
    % Global.exportToPPTX('addtext','Histogram', 'Position',[8.3 5.8  2 0.1], 'FontSize',14,'FontWeight','bold','Color',[1 0 0],'HorizontalAlignment','left');


    % Add to slide
    Global.exportToPPTX('addtable', resultssummary, 'Position', [3.3 1.4 5 4], ...
        'Vert', 'middle', 'Horiz', 'center', 'FontSize', 15);


    % 8. CREATE MONTAGES WITH 2 ROWS
    close all;
    x0 = 8.5; w = 7.3; yStep = 1.5; y0 = 0.05; 
    
    figXe  = makeOneRowMontage(Diffusion.Image(:,:,:,1), '129Xe Diffusion (b0)');
    % numslice = size(Diffusion.Image,3);
    % if numslice < 10
    %     scaleFactore = 1;
    % elseif numslice < 12
    %     scaleFactore = 0.75;
    % elseif numslice < 16
    %     scaleFactore = 0.9;
    % else 
    %     scaleFactore = 1;
    % end
    scaleFactore = 1;

    Pr = Diffusion.Mask_Diff_boundaries;
    if ndims(Pr) == 4  && size(Pr,4) == 3% grayscale: [X Y Z 1]
        Pr = permute(Pr, [1 2 4 3]);  % 
        Pr = squeeze(Pr);  % remove singleton dimension
    end
    figMsk = makeOneRowMontage(Pr, 'Mask Overlay');
    
    figADC = makeOneRowMontage(Diffusion.ADCcoloredmap, 'ADC Map');
    
    LB_Diffmap = Diffusion.LB_colorDiffmap;  % fallback
    if ndims(LB_Diffmap) == 4  && size(LB_Diffmap,4) == 3% grayscale: [X Y Z 1]
        LB_Diffmap = permute(LB_Diffmap, [1 2 4 3]);  % 
        LB_Diffmap = squeeze(LB_Diffmap);  % remove singleton dimension
    end

    figLB_Diffmap = makeOneRowMontage(LB_Diffmap, 'LB_Diffmap');
    
    XePos  = get(figXe, 'Position');
    MskPos  = get(figMsk, 'Position');
    ADCPos = get(figADC, 'Position');
    LB_DiffmapPos = get(figLB_Diffmap, 'Position');
    
    % foldername = "ADC_Analysis\";
    % Venthist = imread([pptDir,'\', char(foldername) ,'TH_Histogram.png']);    
    % Global.exportToPPTX('addpicture',Venthist,'Position',[3.3 5.7 5.5 3]);
    % Global.exportToPPTX('addtext','Histogram', 'Position',[3.3 5.7 2 0.1], 'FontSize',14,'FontWeight','bold','Color',[1 0 0],'HorizontalAlignment','left');
    % 
    % 9. ADD IMAGES TO SLIDE
    Global.exportToPPTX('addpicture', figXe,  'Position', [x0 y0         scaleFactore*w scaleFactore*w*XePos(4)/XePos(3)]);
    Global.exportToPPTX('addpicture', figMsk,  'Position', [x0 y0+yStep   scaleFactore*w scaleFactore*w*MskPos(4)/MskPos(3)]);
    Global.exportToPPTX('addpicture', figADC, 'Position', [x0 y0+2*yStep scaleFactore*w scaleFactore*w*ADCPos(4)/ADCPos(3)]);
    Global.exportToPPTX('addpicture', figLB_Diffmap, 'Position', [x0 y0+3*yStep scaleFactore*w scaleFactore*w* LB_DiffmapPos(4)/LB_DiffmapPos(3)]);
    
    
    Global.exportToPPTX('addtext','129Xe Image (b0)', 'Position',[x0 y0 2 0.1], 'FontSize',14,'FontWeight','bold','Color',[1 0 0],'HorizontalAlignment','left');
    Global.exportToPPTX('addtext','Mask Overlay', 'Position',[x0 y0+yStep 2 0.1], 'FontSize',14,'FontWeight','bold','Color',[1 0 0],'HorizontalAlignment','left');
    Global.exportToPPTX('addtext','ADC Map', 'Position',[x0 y0+2*yStep 2 0.1], 'FontSize',14,'FontWeight','bold','Color',[1 0 0],'HorizontalAlignment','left');
    Global.exportToPPTX('addtext','LB ADC Map', 'Position',[x0 y0+3*yStep 2 0.1], 'FontSize',14,'FontWeight','bold','Color',[1 0 0],'HorizontalAlignment','left');
    
    Global.exportToPPTX('addtext',['Report location:',pptDir], 'Position',[0 8.75 16 0.3], 'FontSize',12,'FontWeight','bold','Color',[0 0 1],'HorizontalAlignment','left');
    
    % Red box (as an empty text box with red fill/border)
    Global.exportToPPTX('addtext',' ', ...
        'Position',[0.1 7.5 0.3 0.2], ...     % [x y w h] – adjust as needed
        'BackgroundColor',abnormalColor, ...
        'LineColor',[1 1 1], 'LineWidth',1);
    
    % Label next to it
    Global.exportToPPTX('addtext','Abnormal', ...
        'Position',[0.35 7.4 1.5 0.3], ...
        'FontSize',16, 'FontWeight','bold', ...
        'Color',[0 0 0], 'HorizontalAlignment','left');

    % add Nb-values
    vals = Diffusion.bvalues;          % e.g. [0 6.25 12.5 18.75 25]
    valStr = strjoin(string(vals), ',');   % => "0,6.25,12.5,18.75,25"
    bvaluelabel = sprintf('%d [%s]', numel(vals), valStr);
    Global.exportToPPTX('addtext',bvaluelabel, ...
        'Position',[5 1.8 3.2 0.3], ...
        'FontSize',16, 'FontWeight','normal', ...
        'Color',[0 0 0],'BackgroundColor',[1 1 1], 'HorizontalAlignment','left');

    % add SNR
    SNRvals  = round(table2array(Diffusion.SNR_table(:,2)),1);
    SNRStr   = strjoin(string(SNRvals), ',');   
    SNRlabel = sprintf('%.1f [%s]', mean(SNRvals), SNRStr);
    Global.exportToPPTX('addtext',SNRlabel, ...
        'Position',[5 2.2 3.2 0.3], ...
        'FontSize',16, 'FontWeight','normal', ...
        'Color',[0 0 0],'BackgroundColor',[1 1 1], 'HorizontalAlignment','left');


    % Get full path of the current function file
    fullFuncPath = mfilename('fullpath');
    funcFolder = fileparts(fullFuncPath);
    % funcFolder = 'D:\Github\XIPline\+DiffusionFunctions';
    colorbarPath = fullfile(funcFolder, 'LBmapColorbar.png');
    colorbarImg = imread(colorbarPath);    
    Global.exportToPPTX('addpicture',colorbarImg,'Position',[15.5 y0+3*yStep 0.5 scaleFactore*w*LB_DiffmapPos(4)/LB_DiffmapPos(3)]);    

    colorbarPath = fullfile(funcFolder, 'ADCmapColorbar.png');
    colorbarImg = imread(colorbarPath);    
    Global.exportToPPTX('addpicture',colorbarImg,'Position',[15.5 y0+2*yStep 0.5 scaleFactore*w*ADCPos(4)/ADCPos(3)]);    


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
    pause(2);  % <-- allow time for file to be written
    presentation.Close();
    ppt.Quit();
    delete(ppt);
end


function figHandle = makeOneRowMontage(vol, figName)

    % Detect grayscale or RGB
    sz = size(vol);

    if numel(sz) == 3
        % Grayscale: X × Y × Z
        isRGB = false;
        numSlices = sz(3);
    elseif numel(sz) == 4 && sz(3)==3
        % RGB stored as X × Y × 3 × Z
        isRGB = true;
        numSlices = sz(4);
    elseif numel(sz) == 4 && sz(4)==3
        % RGB stored as X × Y × Z × 3
        isRGB = true;
        numSlices = sz(3);
        vol = permute(vol, [1 2 4 3]); % convert to X×Y×3×Z
    else
        error('Unsupported volume dimensions. Must be grayscale or RGB.');
    end

    % Create figure
    figHandle = figure('Name', figName, 'Visible', 'off');

    % ---- Display montage in 1 row ----
    if isRGB
        montage(vol, 'Size', [1 numSlices]);

    else
        % reshape grayscale to X × Y × 1 × Z
        vol4d = reshape(vol, [sz(1), sz(2), 1, numSlices]);
        montage(vol4d, 'Size', [1 numSlices], 'DisplayRange', []);
    end

    % ---- Adjust figure to fit montage ----
    set(gca, 'Units', 'pixels');
    axPos = get(gca, 'Position');

    set(figHandle, 'Units', 'pixels', ...
                   'Position', [100 100 axPos(3) axPos(4)]);

    set(gca, 'Units', 'normalized', 'Position', [0 0 1 1]);

end

