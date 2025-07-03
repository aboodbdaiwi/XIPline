function [Ventilation] = GLRLM_Analysis(Ventilation)
%GRAYCOPROPS Properties of gray-level run-length matrix.
%  -------------------------------------------
%  STATS = GRAYCOPROPS(GLRLM,PROPERTIES) Each element in  GLRLM, (r,c),
%   is the probability occurrence of pixel having gray level values r, run-length c in the image.
%   GRAYCOPROPS is to calculate PROPERTIES.
%  -------------------------------------------
%  Requirements:
%  -------------------------------------------
%   GLRLM mustbe an cell array of valid gray-level run-length
%   matrices.Recall that a valid glrlm must be logical or numerical.
%  -------------------------------------------
%  Current supported statistics include:
%  -------------------------------------------
%   Short Run Emphasis (SRE)
%   Long Run Emphasis (LRE)
%   Gray-Level Nonuniformity (GLN)
%   Run Length Nonuniformity (RLN)
%   Run Percentage (RP)
%   Low Gray-Level Run Emphasis (LGRE)
%   High Gray-Level Run Emphasis (HGRE)
%   Short Run Low Gray-Level Emphasis (SRLGE)
%   Short Run High Gray-Level Emphasis (SRHGE)
%   Long Run Low Gray-Level Emphasis (LRLGE)
%   Long Run High Gray-Level Emphasis (LRHGE)
%  --------------------------------------------
%                    Angle     OFFSET
%                    -----     ------
%                    0          1
%                    45         2
%                    90         3
%                    135        4
%  --------------------------------------------
%  Reference:
%  --------------------------------------------
%   Xiaoou Tang,Texture Information in Run-Length Matrices
%   IEEE TRANSACTIONS ON IMAGE PROCESSING, VOL.7, NO.11,NOVEMBER 1998

% this code is based on the GLRLM analysis 
% https://www.mathworks.com/matlabcentral/fileexchange/17482-gray-level-run-length-matrix-toolbox
% Abdullah Bdaiwi
parentPath = Ventilation.parentPath;
cd(parentPath)
maskarray = double(Ventilation.LungMask + Ventilation.VesselMask);
maskarray(maskarray > 1) = 0;
maskarray = double(maskarray);

% remove defects from image
VentLungMask = double(Ventilation.defectMap_forGLRLM == 1);

Image = Ventilation.Image.*maskarray.*VentLungMask;
%figure; imslice(Image);
foldername = "GLRLM_Analysis\";
mkdir(foldername)
outputPath = char(foldername);
cd(outputPath)
disp('GLRLM analysis starting...');
Image = (Image - min(Image(:)))/(max(Image(:)) - min(Image(:)));
GLRLM_output.GLRLM_slice = cell(1,4,size(Image,3));  
GLRLM_output.SI_slice = zeros(size(Image));
GLRLM_output.stats_slice = zeros(4,11,size(Image,3));
GLRLM = cell(1,4);
SI = zeros(size(Image,1),size(Image,2));
stats = zeros(4,11);
for i = 1:size(Image,3)
    I = Image(:,:,i);
    disp(['process slice number ',num2str(i)])
    if sum(I(:)) > 0  
        [GLRLM,SI]= VentilationFunctions.grayrlmatrix(I,'OFFSET',[1;2;3;4],'NumLevels',8,'G',[]);
        stats = VentilationFunctions.grayrlprops(GLRLM);  
    else 
        GLRLM = cell(1,4);  
        SI = zeros(size(Image(:,:,1)));
        stats = zeros(4,11);
    end
    GLRLM_output.GLRLM_slice(:,:,i) = GLRLM;
    GLRLM_output.SI_slice(:,:,i) = SI;
    GLRLM_output.stats_slice(:,:,i) = stats;       
end

% find the mean across directions 
Mean_dir = mean(GLRLM_output.stats_slice,1);
GLRLM_output.Mean_dir = reshape(Mean_dir,[11,size(Image,3)]);

% find the mean across slices
[ii,~,v] = find(GLRLM_output.Mean_dir);
GLRLM_output.stats_average = accumarray(ii,v,[],@mean);
stats_average = GLRLM_output.stats_average;
%% save results
% folderName = "Texture Analysis\";
state_slices = GLRLM_output.stats_slice;
Mean_dir = GLRLM_output.Mean_dir;
for j = 1:4
    for i = 1:11
        matriceslice = state_slices(j,i,:);
        matriceslice = matriceslice - min(matriceslice(:));
        state_slices(j,i,:) = matriceslice / max(matriceslice(:));
        
        matricemean = Mean_dir(i,:);
        matricemean = matricemean - min(matricemean(:));
        Mean_dir(i,:) = matricemean / max(matricemean(:));        
    end 
end
cmap = colormap (jet);
cmap(1,:)=0;
close all;
%------------------------Angle 0 --------------------------
Angel1 = squeeze(state_slices(1,:,:));
figure('position',[350 350 900 350]);
imagesc(Angel1)
title('Angle = 0','Color','black')
xlabel('Slices')
yticks([1 2 3 4 5 6 7 8 9 10 11])
yticklabels({'Short Run Emphasis (SRE)',...
    'Long Run Emphasis (LRE)',...
    'Gray-Level Nonuniformity (GLN)',...
    'Run Length Nonuniformity (RLN)',...
    'Run Percentage (RP)',...
    'Low Gray-Level Run Emphasis (LGRE)',...
    'High Gray-Level Run Emphasis (HGRE)',...
    'Short Run Low Gray-Level Emphasis (SRLGE)',...
    'Short Run High Gray-Level Emphasis (SRHGE)',...
    'Long Run Low Gray-Level Emphasis (LRLGE)',...
    'Long Run High Gray-Level Emphasis (LRHGE)'})
colormap(cmap); colorbar;
% print([folderName + 'GLRLM_Angle0'],'-dpng','-r300');
saveas(gca,'GLRLM_Angle0.png');
close all;
%------------------------Angle 45 --------------------------
Angel2 = squeeze(state_slices(2,:,:));
figure('position',[350 350 900 350]);
imagesc(Angel2)
title('Angle = 45','Color','black')
xlabel('Slices')
yticks([1 2 3 4 5 6 7 8 9 10 11])
yticklabels({'Short Run Emphasis (SRE)',...
    'Long Run Emphasis (LRE)',...
    'Gray-Level Nonuniformity (GLN)',...
    'Run Length Nonuniformity (RLN)',...
    'Run Percentage (RP)',...
    'Low Gray-Level Run Emphasis (LGRE)',...
    'High Gray-Level Run Emphasis (HGRE)',...
    'Short Run Low Gray-Level Emphasis (SRLGE)',...
    'Short Run High Gray-Level Emphasis (SRHGE)',...
    'Long Run Low Gray-Level Emphasis (LRLGE)',...
    'Long Run High Gray-Level Emphasis (LRHGE)'})
colormap(cmap); colorbar;
% print([folderName + 'GLRLM_Angle45'],'-dpng','-r300');
saveas(gca,'GLRLM_Angle45.png');
close all;
%------------------------Angle 90  --------------------------
Angel3 = squeeze(state_slices(3,:,:));
figure('position',[350 350 900 350]);
imagesc(Angel3)
title('Angle = 90','Color','black')
xlabel('Slices')
yticks([1 2 3 4 5 6 7 8 9 10 11])
yticklabels({'Short Run Emphasis (SRE)',...
    'Long Run Emphasis (LRE)',...
    'Gray-Level Nonuniformity (GLN)',...
    'Run Length Nonuniformity (RLN)',...
    'Run Percentage (RP)',...
    'Low Gray-Level Run Emphasis (LGRE)',...
    'High Gray-Level Run Emphasis (HGRE)',...
    'Short Run Low Gray-Level Emphasis (SRLGE)',...
    'Short Run High Gray-Level Emphasis (SRHGE)',...
    'Long Run Low Gray-Level Emphasis (LRLGE)',...
    'Long Run High Gray-Level Emphasis (LRHGE)'})
colormap(cmap); colorbar;
% print([folderName + 'GLRLM_Angle90'],'-dpng','-r300');
saveas(gca,'GLRLM_Angle90.png');

close all;
%------------------------Angle 135  --------------------------
Angel4 = squeeze(state_slices(4,:,:));
figure('position',[350 350 900 350]);
imagesc(Angel4)
title('Angle = 135','Color','black')
xlabel('Slices')
yticks([1 2 3 4 5 6 7 8 9 10 11])
yticklabels({'Short Run Emphasis (SRE)',...
    'Long Run Emphasis (LRE)',...
    'Gray-Level Nonuniformity (GLN)',...
    'Run Length Nonuniformity (RLN)',...
    'Run Percentage (RP)',...
    'Low Gray-Level Run Emphasis (LGRE)',...
    'High Gray-Level Run Emphasis (HGRE)',...
    'Short Run Low Gray-Level Emphasis (SRLGE)',...
    'Short Run High Gray-Level Emphasis (SRHGE)',...
    'Long Run Low Gray-Level Emphasis (LRLGE)',...
    'Long Run High Gray-Level Emphasis (LRHGE)'})
colormap(cmap); colorbar;
% print([folderName + 'GLRLM_Angle135'],'-dpng','-r300');
saveas(gca,'GLRLM_Angle135.png');

close all;
%--------------------- Mean angles --------------------
figure('position',[350 350 900 350]);
imagesc(Mean_dir)
title('Mean Angles','Color','black')
xlabel('Slices')
yticks([1 2 3 4 5 6 7 8 9 10 11])
yticklabels({'Short Run Emphasis (SRE)',...
    'Long Run Emphasis (LRE)',...
    'Gray-Level Nonuniformity (GLN)',...
    'Run Length Nonuniformity (RLN)',...
    'Run Percentage (RP)',...
    'Low Gray-Level Run Emphasis (LGRE)',...
    'High Gray-Level Run Emphasis (HGRE)',...
    'Short Run Low Gray-Level Emphasis (SRLGE)',...
    'Short Run High Gray-Level Emphasis (SRHGE)',...
    'Long Run Low Gray-Level Emphasis (LRLGE)',...
    'Long Run High Gray-Level Emphasis (LRHGE)'})
colormap(cmap); colorbar;
% print([folderName + 'GLRLM_Mean Angles'],'-dpng','-r300');
saveas(gca,'GLRLM_Mean Angles.png');

close all;
%--------------------- Average of all slices and angles --------------------
figure('position',[350 350 600 350]);
imagesc(stats_average)
title('Average All','Color','black')
yticks([1 2 3 4 5 6 7 8 9 10 11])
yticklabels({['Short Run Emphasis (SRE)=',num2str(stats_average(1))],...
    ['Long Run Emphasis (LRE)=',num2str(stats_average(2))],...
    ['Gray-Level Nonuniformity (GLN)=',num2str(stats_average(3))],...
    ['Run Length Nonuniformity (RLN)=',num2str(stats_average(4))],...
    ['Run Percentage (RP)=',num2str(stats_average(5))],...
    ['Low Gray-Level Run Emphasis (LGRE)=',num2str(stats_average(6))],...
    ['High Gray-Level Run Emphasis (HGRE)=',num2str(stats_average(7))],...
    ['Short Run Low Gray-Level Emphasis (SRLGE)=',num2str(stats_average(8))],...
    ['Short Run High Gray-Level Emphasis (SRHGE)=',num2str(stats_average(9))],...
    ['Long Run Low Gray-Level Emphasis (LRLGE)=',num2str(stats_average(10))],...
    ['Long Run High Gray-Level Emphasis (LRHGE)=',num2str(stats_average(11))]})
colormap(cmap); colorbar;
% print([folderName + 'GLRLM_AverageAll'],'-dpng','-r300');
saveas(gca,'GLRLM_AverageAll.png');

close all;
%% Save summary into xlsx file
GLRLMSummary = table({...
    'Short Run Emphasis (SRE)';...
    'Long Run Emphasis (LRE)';...
    'Gray-Level Nonuniformity (GLN)';...
    'Run Length Nonuniformity (RLN)';...
    'Run Percentage (RP)';...
    'Low Gray-Level Run Emphasis (LGRE)';...
    'High Gray-Level Run Emphasis (HGRE)';...
    'Short Run Low Gray-Level Emphasis (SRLGE)';...
    'Short Run High Gray-Level Emphasis (SRHGE)'; ...
    'Long Run Low Gray-Level Emphasis (LRLGE)'; ...
    'Long Run High Gray-Level Emphasis (LRHGE)'},...
    [ stats_average(1); stats_average(2); stats_average(3); stats_average(4);...
    stats_average(5); stats_average(6); stats_average(7); stats_average(8); stats_average(9);...
    stats_average(10); stats_average(11)]);
headers =({'summary','Mean'});
GLRLMSummary.Properties.VariableNames = headers;
writetable(GLRLMSummary,[parentPath, outputPath,'GLRLMSummary.xlsx'],'Sheet',1)

%% save ppt 
cd(parentPath)
%Start new presentation
isOpen  = Global.exportToPPTX(); 
if ~isempty(isOpen) %If PowerPoint already started, then close first and then open a new one
    Global.exportToPPTX('close');
end
% Generate filename with today's date
today = datetime('today');
date_str = datestr(today, 'yyyymmdd');
ppt_file_name = ['Ventilation_Analysis_', date_str, '.pptx'];
ReportTitle = ['Ventilation_Analysis_', date_str];
% Create or open the presentation
if isfile(ppt_file_name)
    disp('File existed')
    Global.exportToPPTX('open', ppt_file_name);
    Global.exportToPPTX('switchslide', 1);
else            
    Global.exportToPPTX('new', 'Dimensions', [16 9], ...
        'Title', ReportTitle, ...
        'Author', 'CPIR @ CCHMC');
end
%Add slides
Global.exportToPPTX('addslide'); % angles
Global.exportToPPTX('addtext',sprintf(foldername),'Position',[6 0 7 1],'Color','b','FontSize',25);

GLRLM_Angle0 = imread([parentPath, outputPath ,'GLRLM_Angle0.png']);    
Global.exportToPPTX('addpicture',GLRLM_Angle0,'Position',[0 0.5 8.5 8.5*(size(GLRLM_Angle0,1)/size(GLRLM_Angle0,2))]);

GLRLM_Angle45 = imread([parentPath, outputPath ,'GLRLM_Angle45.png']);    
Global.exportToPPTX('addpicture',GLRLM_Angle45,'Position',[8 0.5 8.5 8.5*(size(GLRLM_Angle45,1)/size(GLRLM_Angle45,2))]);

GLRLM_Angle90 = imread([parentPath, outputPath ,'GLRLM_Angle90.png']);    
Global.exportToPPTX('addpicture',GLRLM_Angle90,'Position',[0 5 8.5 8.5*(size(GLRLM_Angle90,1)/size(GLRLM_Angle90,2))]);

GLRLM_Angle135 = imread([parentPath, outputPath ,'GLRLM_Angle135.png']);    
Global.exportToPPTX('addpicture',GLRLM_Angle135,'Position',[8 5 8.5 8.5*(size(GLRLM_Angle135,1)/size(GLRLM_Angle135,2))]);

Global.exportToPPTX('addslide'); % average
Global.exportToPPTX('addtext',sprintf(foldername),'Position',[6 0 7 1],'Color','b','FontSize',25);

GLRLM_MeanAngles = imread([parentPath, outputPath ,'GLRLM_Mean Angles.png']);    
Global.exportToPPTX('addpicture',GLRLM_MeanAngles,'Position',[0.5 0.5 8.5 8.5*(size(GLRLM_MeanAngles,1)/size(GLRLM_MeanAngles,2))]);

GLRLM_AverageAll = imread([parentPath, outputPath ,'GLRLM_AverageAll.png']);    
Global.exportToPPTX('addpicture',GLRLM_AverageAll,'Position',[0.5 4 8.5 8.5*(size(GLRLM_AverageAll,1)/size(GLRLM_AverageAll,2))]);


Global.exportToPPTX('save',fullfile(parentPath, ReportTitle));
Global.exportToPPTX('close');
fprintf('PowerPoint file has been saved\n');
            Global.exportToPPTX('save',fullfile(parentPath, ReportTitle));
            Global.exportToPPTX('close');
            fprintf('PowerPoint file has been saved\n');               

close all;

%% %% %% save maps in mat file
Ventilation.GLRLM = GLRLM_output;
Ventilation.SRE = stats_average(1);
Ventilation.LRE = stats_average(2);
Ventilation.GLN = stats_average(3);
Ventilation.RLN = stats_average(4);
Ventilation.RP = stats_average(5);
Ventilation.LGRE = stats_average(6);
Ventilation.HGRE = stats_average(7);
Ventilation.SRLGR = stats_average(8);
Ventilation.SRHGE = stats_average(9);
Ventilation.LRLGE = stats_average(10);
Ventilation.LRHGE = stats_average(11);
    
save_data=[parentPath, outputPath,'GLRLM_Analysis','.mat'];
save(save_data,'GLRLM_output');  

disp('GLRLM analysis completed');
end


