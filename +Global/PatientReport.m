
function PatientReport(MainInput)
%   Inputs:
%          
%   Outputs:
%     
%
%   Example: 
%   Package: 
%
%   Author: Abdullah Bdaiwi 
%   Work email: abdullah.bdaiwi@cchmc.org
%   Personal email: abdaiwi89@gmail.com
%   Website: https://cpir.cchmc.org/

%%

cd(MainInput.XeDataLocation);
%Start new presentation
isOpen  = Global.exportToPPTX();
if ~isempty(isOpen) %If PowerPoint already started, then close first and then open a new one
    Global.exportToPPTX('close');
end

switch MainInput.AnalysisType
    case 'Ventilation'
        ReportTitle='Ventilation_Analysis';
        ppt_file_name = 'Ventilation_Analysis.pptx';
        analysis_folder = '\Ventilation Analysis';
    case 'Diffusion'
        ReportTitle='Diffusion_Analysis';
        ppt_file_name = 'Diffusion_Analysis.pptx';
        analysis_folder = '\Diffusion Analysis';
    case 'GasExchange'
        ReportTitle='GasExchange_Analysis';
        ppt_file_name = 'GasExchange_Analysis.pptx';
        analysis_folder = '\Gas Exchange Analysis';
end 
outputpath = [MainInput.XeDataLocation,analysis_folder];
cd(outputpath)

if exist(ppt_file_name, 'file') == 2            
    Global.exportToPPTX('open',ppt_file_name);
    Global.exportToPPTX('switchslide',1);
else            
    Global.exportToPPTX('new','Dimensions',[16 9], ...
        'Title',ReportTitle, ...
        'Author','CPIR @ CCHMC');
end 
Global.exportToPPTX('addslide','Position',1); %Title Slide
% PatientReportSummary = imread([MainInput.XeDataLocation,'\PatientReportSummary.png']);    
% Global.exportToPPTX('addpicture',PatientReportSummary,'Position',[0.5 0.5 13 13*(size(PatientReportSummary,1)/size(PatientReportSummary,2))]);
% add titles 
Global.exportToPPTX('addtext',sprintf(['Patient Report (',analysis_folder(2:end),')']),'Position', [4 0.1 10 0.7],'Color','b','FontSize',35);
Global.exportToPPTX('addtext',sprintf('Study ID: '),'Position',             [0.5 1 3 0.5],'Color','r','FontSize',25);
Global.exportToPPTX('addtext',sprintf('Study Type: '),'Position',           [0.5 1.5 3 0.5],'Color','r','FontSize',25);
Global.exportToPPTX('addtext',sprintf('Patient ID: '),'Position',           [0.5 2 3 0.5],'Color','r','FontSize',25);
Global.exportToPPTX('addtext',sprintf('Scan Date: '),'Position',            [0.5 2.5 3 0.5],'Color','r','FontSize',25);
Global.exportToPPTX('addtext',sprintf('Xe Dose Volume (mL): '),'Position',  [0.5 3 4 0.5],'Color','r','FontSize',25);
Global.exportToPPTX('addtext',sprintf('Patient Name: '),'Position',         [0.5 3.5 3 0.5],'Color','r','FontSize',25);
Global.exportToPPTX('addtext',sprintf('MRM number: '),'Position',           [0.5 4 3 0.5],'Color','r','FontSize',25);
Global.exportToPPTX('addtext',sprintf('Sex: '),'Position',                  [0.5 4.5 3 0.5],'Color','r','FontSize',25);
Global.exportToPPTX('addtext',sprintf('Age (year): '),'Position',           [0.5 5 3 0.5],'Color','r','FontSize',25);
Global.exportToPPTX('addtext',sprintf('Height (cm): '),'Position',          [0.5 5.5 3 0.5],'Color','r','FontSize',25);
Global.exportToPPTX('addtext',sprintf('Weight (lbs): '),'Position',         [0.5 6 3 0.5],'Color','r','FontSize',25);
Global.exportToPPTX('addtext',sprintf('Summary of Findings: '),'Position',  [0.5 6.5 4 0.5],'Color','r','FontSize',25);
Global.exportToPPTX('addtext',sprintf('Notes: '),'Position',                [0.5 7 3 0.5],'Color','r','FontSize',25);
Global.exportToPPTX('addtext',sprintf('Data Analyst: '),'Position',         [0.5 7.5 3 0.5],'Color','r','FontSize',25);
Global.exportToPPTX('addtext',sprintf('Processing Date: '),'Position',      [0.5 8 3 0.5],'Color','r','FontSize',25);
% add patient info (answers)
Global.exportToPPTX('addtext',sprintf(MainInput.studyID),'Position',        [4 1 10 0.5],'Color','k','FontSize',25);
Global.exportToPPTX('addtext',sprintf(MainInput.studyType),'Position',      [4 1.5 10 0.5],'Color','k','FontSize',25);
Global.exportToPPTX('addtext',sprintf(MainInput.patientID),'Position',      [4 2 10 0.5],'Color','k','FontSize',25);
Global.exportToPPTX('addtext',sprintf(MainInput.scanDate),'Position',       [4 2.5 10 0.5],'Color','k','FontSize',25);
Global.exportToPPTX('addtext',sprintf(MainInput.xeDoseVolume),'Position',   [4 3 10 0.5],'Color','k','FontSize',25);
Global.exportToPPTX('addtext',sprintf(MainInput.PatientName),'Position',    [4 3.5 10 0.5],'Color','k','FontSize',25);
Global.exportToPPTX('addtext',sprintf(MainInput.MRMnumber),'Position',      [4 4 10 0.5],'Color','k','FontSize',25);
Global.exportToPPTX('addtext',sprintf(MainInput.sex),'Position',            [4 4.5 10 0.5],'Color','k','FontSize',25);
Global.exportToPPTX('addtext',sprintf(MainInput.age),'Position',            [4 5 10 0.5],'Color','k','FontSize',25);
Global.exportToPPTX('addtext',sprintf(MainInput.height),'Position',         [4 5.5 10 0.5],'Color','k','FontSize',25);
Global.exportToPPTX('addtext',sprintf(MainInput.weight),'Position',         [4 6 10 0.5],'Color','k','FontSize',25);
Global.exportToPPTX('addtext',sprintf(MainInput.summaryofFindings),'Position', [4 6.5 10 0.5],'Color','k','FontSize',25);
Global.exportToPPTX('addtext',sprintf(MainInput.notes),'Position',          [4 7 10 0.5],'Color','k','FontSize',25);
Global.exportToPPTX('addtext',sprintf(MainInput.dataAnalyst),'Position',    [4 7.5 10 0.5],'Color','k','FontSize',25);
Global.exportToPPTX('addtext',sprintf(MainInput.processingDate),'Position', [4 8 10 0.5],'Color','k','FontSize',25);


% save all
Global.exportToPPTX('save',fullfile(outputpath, ReportTitle));
Global.exportToPPTX('close');      

%% save report as a PDF
ppt = actxserver('PowerPoint.Application');
presentation = ppt.Presentations.Open([outputpath,'\',ppt_file_name]);
outputPath = [outputpath,'\',ReportTitle,'_Patient_Report.pdf'];
presentation.SaveAs(outputPath, 32);
presentation.Close();
ppt.Quit();
delete(ppt);

end






