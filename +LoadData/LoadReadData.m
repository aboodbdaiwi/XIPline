
function [Ventilation, Diffusion, GasExchange, Proton] = LoadReadData(MainInput)
%   Inputs:
%      
%   Outputs:
%                   
%   Package: https://github.com/aboodbdaiwi/HP129Xe_Analysis_App
%
%   Author: Abdullah S. Bdaiwi
%   Work email: abdullah.bdaiwi@cchmc.org
%   Personal email: abdaiwi89@gmail.com
%   Website: https://www.cincinnatichildrens.org/research/divisions/c/cpir
%
%   Please add updates at the end. Ex: 3/10/24 - ASB: update .... 

% Initialize output variables
Ventilation.Image = [];
Ventilation.filename = [];
Ventilation.folder = [];

Diffusion.Image = [];
Diffusion.filename = [];
Diffusion.folder = [];

GasExchange.UncorrectedVentImage = [];
GasExchange.VentImage = [];
GasExchange.GasImage = [];
GasExchange.DissolvedImage = [];
GasExchange.CorrDissolvedImage = [];
GasExchange.AppendedDissolvedNMRFit = [];
GasExchange.RBC2Bar_struct = [];
GasExchange.RBCOsc_High_Image = [];
GasExchange.RBCOsc_Low_Image = [];
GasExchange.RBCOsc_Normalization = [];
GasExchange.ActTE90 = [];
GasExchange.DisFlipAngle = [];
GasExchange.PixelShift = [];
GasExchange.freq_jump = [];
GasExchange.DissolvedNMR = [];
GasExchange.PixelShift = [];

Proton.Image = [];
Proton.filename = [];
Proton.folder = [];
Proton.H_RecMatrix = [];
Proton.ProtonMax = [];

%% Load/Read Xenon data 
cd(MainInput.XeDataLocation)
if strcmp(MainInput.AnalysisType,'Ventilation')                 
    mkdir([MainInput.XeDataLocation '\Ventilation Analysis']);
    outputpath = [MainInput.XeDataLocation '\Ventilation Analysis'];
    Ventilation.outputpath = outputpath;    
elseif strcmp(MainInput.AnalysisType,'Diffusion')
    mkdir([MainInput.XeDataLocation '\Diffusion Analysis']);
    outputpath = [MainInput.XeDataLocation '\Diffusion Analysis'];
    Diffusion.outputpath = outputpath;
elseif strcmp(MainInput.AnalysisType,'GasExchange')
    mkdir([MainInput.XeDataLocation '\Gas Exchange Analysis']);

    mkdir([MainInput.XeDataLocation '\Gas Exchange Analysis']);
    outputpath = [MainInput.XeDataLocation '\Gas Exchange Analysis'];
    GasExchange.outputpath = outputpath;
end

if strcmp(MainInput.XeDataext,'.dcm') == 1         
    [Image, file_folder, FileNames] = LoadData.DICOM_Load(MainInput.XeDataLocation);    
    if strcmp(MainInput.AnalysisType,'Ventilation') == 1                 
        Ventilation.Image = Image;
        Ventilation.filename = FileNames;
        Ventilation.folder = file_folder;
    elseif strcmp(MainInput.AnalysisType,'Diffusion') == 1 
        if length(size(Image)) == 3
            try 
                Nb = MainInput.Nbvalues;
                switch MainInput.DiffAcqOrder 
                    case 'b-value interleave'
                        Image = reshape(Image, [size(Image,1),size(Image,2),Nb,size(Image,3)/Nb]);
                        Image = permute(Image, [1 2 4 3]);
                    case 'slice interleave'
                        Image = reshape(Image, [size(Image,1),size(Image,2),size(Image,3)/Nb,Nb]);
                end   
%                 for i = 1:length(FileNames)
%                     img = squeeze(Image(:,:,i));
%                     SignalMean(i) = mean(img(:));
%                 end
%                 [~,locs] =findpeaks(SignalMean);
%                 Nb = length(locs)+1; 
%                 Image = reshape(Image, [size(Image,1),size(Image,2),size(Image,3)/Nb,Nb]);
%                 Image = permute(Image, [1 2 4 3]);
            catch
               disp('b value dimension cannot be determined');
            end
        end
        Diffusion.Image = Image;
        Diffusion.filename = FileNames;
        Diffusion.folder = file_folder;

    elseif strcmp(MainInput.AnalysisType,'GasExchange') == 1 
        % not supported yet
    end
elseif strcmp(MainInput.XeDataext,'.mat') == 1 
%         DataFiles = dir([MainInput.XeDataLocation,'\*.mat']);
        file_name = MainInput.XeFileName;
        file_folder = MainInput.XeDataLocation;
%         file_with_path = strcat(file_folder,'\',file_name);  % join path and filename to open
        load(MainInput.XeFullPath);
        file_name2 = file_name; 
        file_name2(end-3:end)=[];    
    if strcmp(MainInput.AnalysisType,'Ventilation') == 1                 
        Ventilation.Image = eval(file_name2);
        Ventilation.filename = file_name;
        Ventilation.folder = file_folder;       
    elseif strcmp(MainInput.AnalysisType,'Diffusion') == 1 
        Diffusion.Image = eval(file_name2);
        Diffusion.filename = file_name;
        Diffusion.folder = file_folder;
    elseif strcmp(MainInput.AnalysisType,'GasExchange') == 1 && strcmp(MainInput.Institute,'CCHMC') == 1  
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
        GasExchange.DissolvedNMR = DissolvedNMR;
        GasExchange.PixelShift = SigDynamics;
    end    
elseif strcmp(MainInput.XeDataext,'.nii') == 1 || strcmp(MainInput.XeDataext,'.gz') == 1 
%         DataFiles = MainInput.XeDataLocation; 
%         if length(DataFiles) == 1
%             DataFiles = dir([MainInput.XeDataLocation,'\*.nii.gz']); 
%         else 
%            DataFiles = dir([MainInput.XeDataLocation,'\*.nii']); 
%         end     
        file_name = MainInput.XeFileName;
        file_folder = MainInput.XeDataLocation;
        try
            A1 = LoadData.load_nii(MainInput.XeFullPath); % Original
        catch
            A1 = LoadData.load_untouch_nii(MainInput.XeFullPath); % Original
        end
        
        A=A1.img;
        A = double(squeeze(A));
        A=imrotate(A,90);
        A=flip(A,2); % Original                   
        Image=A;        
    if strcmp(MainInput.AnalysisType,'Ventilation') == 1                 
        Ventilation.Image = Image;
        Ventilation.filename = file_name;
        Ventilation.folder = file_folder;       
    elseif strcmp(MainInput.AnalysisType,'Diffusion') == 1 
        Diffusion.Image = Image;
        Diffusion.filename = file_name;
        Diffusion.folder = file_folder;
    elseif strcmp(MainInput.AnalysisType,'GasExchange') == 1 
        % not supported yet
    end
elseif (strcmp(MainInput.XeDataext,'.h5') || strcmp(MainInput.XeDataext,'.mrd')) 
    if strcmp(MainInput.AnalysisType,'Ventilation') == 1   
        MainInput.ReconImageMode = 'xenon';
        [Image] = LoadData.ismrmrd.cartesian_2D_recon(MainInput);
        Ventilation.Image = Image;  
    elseif strcmp(MainInput.AnalysisType,'Diffusion') 
        MainInput.ReconImageMode = 'xenon';
        [Image] = LoadData.ismrmrd.cartesian_2D_recon(MainInput);
        Diffusion.Image = Image;
    elseif strcmp(MainInput.AnalysisType,'GasExchange') && strcmp(MainInput.Institute,'XeCTC') && strcmp(MainInput.SequenceType, '3D Radial')
        [GasExchange] = LoadData.ismrmrd.radial_3D_XeCTC_gx_recon(MainInput,GasExchange);
    end
elseif strcmp(MainInput.XeDataext,'.data')  && strcmp(MainInput.Scanner, 'Philips')  
    if strcmp(MainInput.AnalysisType,'Ventilation') == 1 && strcmp(MainInput.SequenceType, '2D GRE') == 1 ...
            && (strcmp(MainInput.ScannerSoftware, '5.3.1') == 1 || strcmp(MainInput.ScannerSoftware, '5.6.1') == 1)
        [Image, file_folder, file_name] = LoadData.LoadData_Gas_VentDiff_Philips_GRE(MainInput.XeDataLocation);
        Ventilation.Image = Image;
        Ventilation.filename = file_name;
        Ventilation.folder = file_folder;     

    elseif strcmp(MainInput.AnalysisType,'Ventilation') == 1 && strcmp(MainInput.SequenceType, '2D GRE') == 1 && strcmp(MainInput.ScannerSoftware, '5.9.0') == 1                
        [Image, file_folder, file_name] = LoadData.LoadData_Gas_VentDiff_Philips_GRE_R590(MainInput);
        Ventilation.Image = Image;
        Ventilation.filename = file_name;
        Ventilation.folder = file_folder;          
        
    elseif strcmp(MainInput.AnalysisType,'Diffusion') == 1  && strcmp(MainInput.SequenceType, '2D GRE') == 1 ...
            && (strcmp(MainInput.ScannerSoftware, '5.3.1') == 1 || strcmp(MainInput.ScannerSoftware, '5.6.1') == 1)
        [Image, file_folder, file_name] = LoadData.LoadData_Gas_VentDiff_Philips_GRE(MainInput);
        Diffusion.Image = Image;
        Diffusion.filename = file_name;
        Diffusion.folder = file_folder;

    elseif strcmp(MainInput.AnalysisType,'Diffusion') == 1  && strcmp(MainInput.SequenceType, '2D GRE') == 1 && strcmp(MainInput.ScannerSoftware, '5.9.0') == 1 
        [Image, file_folder, file_name] = LoadData.LoadData_Gas_VentDiff_Philips_GRE_R590(MainInput);
        Diffusion.Image = Image;
        Diffusion.filename = file_name;
        Diffusion.folder = file_folder;

    elseif strcmp(MainInput.XeDataext,'.data') && strcmp(MainInput.AnalysisType,'GasExchange')  && strcmp(MainInput.SequenceType, '3D Radial')  ...
            && (strcmp(MainInput.ScannerSoftware, '5.6.1')  || strcmp(MainInput.ScannerSoftware, '5.9.0')) &&...
            (strcmp(MainInput.Institute,'CCHMC')  || strcmp(MainInput.Institute,'XeCTC')   )
           [GasExchange] = LoadData.LoadData_Gas_GasExchange_Philips_Sin(MainInput.XeDataLocation,MainInput.Institute,GasExchange);                      
    end

elseif strcmp(MainInput.XeDataext,'.dat') && strcmp(MainInput.Scanner,'Siemens') 
            XeFullPath = MainInput.XeFullPath;
            XeDataLocation = MainInput.XeDataLocation;
            filename = MainInput.XeFileName;
            Xe_name = MainInput.Xe_name;
            ext = MainInput.XeDataext;
            cd(XeDataLocation)
        if strcmp(MainInput.AnalysisType,'Ventilation') == 1    
            LoadData.ismrmrd.Siemens.gre_to_ismrmrd(filename,Xe_name,fullfile(XeDataLocation,[Xe_name '.h5']));
            MainInput.XeFileName = [Xe_name '.h5'];
            MainInput.ReconImageMode = 'xenon';
            [Image] = LoadData.ismrmrd.cartesian_2D_recon(MainInput);
            Ventilation.Image = Image;  
        elseif strcmp(MainInput.AnalysisType,'Diffusion') == 1 
            LoadData.ismrmrd.Siemens.diff_to_ismrmrd(filename,Xe_name,fullfile(XeDataLocation,[Xe_name '.h5']));
            MainInput.XeFileName = [Xe_name '.h5'];
            MainInput.ReconImageMode = 'xenon';
            [Image] = LoadData.ismrmrd.cartesian_2D_recon(MainInput);
            Diffusion.Image = Image;
        elseif strcmp(MainInput.AnalysisType,'GasExchange') == 1 
            try
                LoadData.ismrmrd.Siemens.xpdixon_2303_2_mrd(filename,Xe_name,fullfile(XeDataLocation,[Xe_name '.h5']));
            catch
                LoadData.ismrmrd.Siemens.dissolved_to_ismrmrd(filename,Xe_name,fullfile(XeDataLocation,[Xe_name '.h5']));
            end
            LoadData.ismrmrd.Siemens.calibration_to_ismrmrd(MainInput.CalFileName,MainInput.Cal_name,fullfile(XeDataLocation,[MainInput.Cal_name '_Calibration.h5']));
            
            MainInput.XeFileName = [Xe_name '.h5'];
            MainInput.CalFileName = [MainInput.Cal_name '_Calibration.h5'];            
            [GasExchange] = LoadData.ismrmrd.radial_3D_XeCTC_gx_recon(MainInput,GasExchange);
        end
elseif strcmp(MainInput.XeDataext,'.7') && strcmp(MainInput.Scanner,'GE') 
            XeFullPath = MainInput.XeFullPath;
            XeDataLocation = MainInput.XeDataLocation;
            filename = MainInput.XeFileName;
            Xe_name = MainInput.Xe_name;
            ext = MainInput.XeDataext;
            cd(XeDataLocation)
        if strcmp(MainInput.AnalysisType,'Ventilation')  
            MainInput.ReconImageMode = 'xenon';
            LoadData.ismrmrd.GE.gre_to_ismrmrd(filename,fullfile(XeDataLocation,[Xe_name '.h5']));
            MainInput.XeFileName = [Xe_name '.h5'];
            MainInput.ReconImageMode = 'xenon';
            [Image] = LoadData.ismrmrd.cartesian_2D_recon(MainInput);
            Ventilation.Image = Image;  
        elseif strcmp(MainInput.AnalysisType,'Diffusion') 
            LoadData.ismrmrd.GE.diff_to_ismrmrd(filename,Xe_name,fullfile(XeDataLocation,[Xe_name '.h5']));
            MainInput.XeFileName = [Xe_name '.h5'];
            MainInput.ReconImageMode = 'xenon';
            [Image] = LoadData.ismrmrd.cartesian_2D_recon(MainInput);
            Diffusion.Image = Image;
        elseif strcmp(MainInput.AnalysisType,'GasExchange')  
            mrdfile = [Xe_name '.h5']; %mrdfile,dissolvedarchive,dissolvedmat,calmat
            LoadData.ismrmrd.GE.dissolved_to_ismrmrd(mrdfile,dissolvedarchive,dissolvedmat,calmat);
            LoadData.ismrmrd.GE.calibration_to_ismrmrd(MainInput.CalFileName,MainInput.Cal_name,fullfile(XeDataLocation,[MainInput.Cal_name '_Calibration.h5']));
            MainInput.XeFileName = [Xe_name '.h5'];
            MainInput.CalFileName = [MainInput.Cal_name '_Calibration.h5'];            
            [GasExchange] = LoadData.ismrmrd.radial_3D_XeCTC_gx_recon(MainInput,GasExchange);
        end

%--------------------- add new read load function here --------------------
% elseif strcmp(MainInput.XeDataType,'add DataType') == 1
%     if strcmp(MainInput.AnalysisType,'Ventilation') == 1                 
%       %  add load/read function here 
%     elseif strcmp(MainInput.AnalysisType,'Diffusion') == 1 
%       %  add load/read function here
%     elseif strcmp(MainInput.AnalysisType,'GasExchange') == 1 
%       %  add load/read function here  
%     end    

    % apply denoising 
    if strcmp(MainInput.denoiseXe,'yes')
        if strcmp(MainInput.AnalysisType,'Ventilation')
            Ventilation.Image = (Ventilation.Image - min(Ventilation.Image(:)))./(max(Ventilation.Image(:)) - min(Ventilation.Image(:)));
            for i = 1:size(Ventilation.Image,3)
                    Ventilation.Image(:,:,i) = Global.bm3d.BM3D(squeeze(Ventilation.Image(:,:,i)), MainInput.denoiseSD);
            end
        elseif strcmp(MainInput.AnalysisType,'Diffusion')
            Image = Diffusion.Image;
            Image = (Image - min(Image(:)))./(max(Image(:)) - min(Image(:)));
            for i = 1:size(Image,3)
                for j = 1:size(Image,4)
                    Image(:,:,i,j) = Global.bm3d.BM3D(squeeze(Image(:,:,i,j)), 0.01);
                end
            end
            Diffusion.Image = Image;
        elseif strcmp(MainInput.AnalysisType,'GasExchange')
        % we can't apply denoising on complex data         
        %             GasExchange.VentImage
        %             GasExchange.GasImage
        %             GasExchange.Dissolved Image
        
        end
    end

end 


%% Load/Read Proton data 
%% 
if MainInput.NoProtonImage == 0    
    try 
        cd(MainInput.XeDataLocation)
        if strcmp(MainInput.AnalysisType,'Ventilation')                 
            mkdir([MainInput.HDataLocation '\Ventilation Analysis']);  
        elseif strcmp(MainInput.AnalysisType,'Diffusion')
        elseif strcmp(MainInput.AnalysisType,'GasExchange')          
            mkdir([MainInput.HDataLocation '\Gas Exchange Analysis']);
        end
        if strcmp(MainInput.HDataext,'.dcm')               
            [HImage, file_folder, file_name] = LoadData.DICOM_Load(MainInput.HDataLocation);
            Proton.Image = double(HImage);
            Proton.filename = file_name;
            Proton.folder = file_folder;
            
        elseif strcmp(MainInput.HDataext,'.mat') == 1 
            DataFiles = dir([MainInput.HDataLocation,'\*.mat']);
            file_name = DataFiles.name;
            file_folder = DataFiles.folder;
            file_with_path = strcat(file_folder,'\',file_name);  % join path and filename to open
            load(file_with_path);
            file_name2 = file_name; 
            file_name2(end-3:end)=[];    
            Proton.Image = eval(file_name2);
            Proton.filename = file_name;
            Proton.folder = file_folder;
            
        elseif strcmp(MainInput.HDataext,'.nii') == 1 || strcmp(MainInput.HDataext,'.gz') == 1 
            DataFiles = dir([MainInput.HDataLocation,'\*.nii.gz']); 
            if length(DataFiles) == 1
                DataFiles = dir([MainInput.HDataLocation,'\*.nii.gz']); 
            else 
               DataFiles = dir([MainInput.HDataLocation,'\*.nii']); 
            end        
            file_name = DataFiles.name;
            file_folder = DataFiles.folder;
            A1 = LoadData.load_nii([file_folder ,'\',file_name]); % Original
            A=A1.img;
            A = double(squeeze(A));
            I90=imrotate(A,90);
            Ifv=flip(I90,2); % Original                   
            HImage=Ifv;        
            Proton.Image = HImage;
            Proton.filename = file_name;
            Proton.folder = file_folder;
        elseif strcmp(MainInput.XeDataext,'.data')  && strcmp(MainInput.Scanner, 'Philips')  
            MainInput.ReconImageMode = 'proton';
            if strcmp(MainInput.AnalysisType,'Ventilation') && strcmp(MainInput.SequenceType, '2D GRE') ...
                    && (strcmp(MainInput.ScannerSoftware, '5.3.1') || strcmp(MainInput.ScannerSoftware, '5.6.1'))
                [Image, file_folder, file_name] = LoadData.LoadData_Gas_VentDiff_Philips_GRE(MainInput.XeDataLocation);
                Proton.Image = Image;
                Proton.filename = file_name;
                Proton.folder = file_folder;     
        
            elseif strcmp(MainInput.AnalysisType,'Ventilation') && strcmp(MainInput.SequenceType, '2D GRE') && strcmp(MainInput.ScannerSoftware, '5.9.0')               
                [Image, file_folder, file_name] = LoadData.LoadData_Gas_VentDiff_Philips_GRE_R590(MainInput);
                Proton.Image = Image;
                Proton.filename = file_name;
                Proton.folder = file_folder; 
            end

        elseif strcmp(MainInput.HDataext,'.data') == 1 && strcmp(MainInput.Scanner, 'Philips') == 1 &&...
               strcmp(MainInput.SequenceType, '3D Radial')        
                PixelShift = GasExchange.PixelShift; 
                if sum(PixelShift) > 1
                    PixelShift = GasExchange.PixelShift;
                else
                    PixelShift = [0; 0; 0];
                end
                [Proton] = LoadData.LoadData_Proton_GasExchange_Philips_Sin(MainInput.HDataLocation,PixelShift,MainInput.Institute);
        elseif (strcmp(MainInput.XeDataext,'.h5') || strcmp(MainInput.XeDataext,'.mrd'))              
            if strcmp(MainInput.AnalysisType,'Ventilation')  
                MainInput.ReconImageMode = 'proton';
                [Image] = LoadData.ismrmrd.cartesian_2D_recon(MainInput);
                Proton.Image = Image;  
            elseif strcmp(MainInput.AnalysisType,'Diffusion') 
                %
            elseif strcmp(MainInput.AnalysisType,'GasExchange') && strcmp(MainInput.Institute,'XeCTC') && strcmp(MainInput.SequenceType, '3D Radial')
                [Proton] = LoadData.ismrmrd.radial_3D_XeCTC_H_recon(MainInput, GasExchange, Proton);
            end
        elseif strcmp(MainInput.HDataext,'.dat') && strcmp(MainInput.Scanner,'Siemens') 
                    HFullPath = MainInput.HFullPath;
                    HDataLocation = MainInput.XeDataLocation;
                    HFileName = MainInput.HFileName;
                    H_name = MainInput.H_name;
                    H_ext = MainInput.HDataext;
                    cd(HDataLocation)                  
                if strcmp(MainInput.AnalysisType,'Ventilation')    
                    LoadData.ismrmrd.Siemens.gre_to_ismrmrd(HFileName,H_name,fullfile(HDataLocation,[H_name '.h5']));
                    MainInput.HFileName = [H_name '.h5'];
                    MainInput.ReconImageMode = 'proton';
                    [Image] = LoadData.ismrmrd.cartesian_2D_recon(MainInput);
                    Proton.Image = Image;  
                elseif strcmp(MainInput.AnalysisType,'Diffusion') 
                    % not supported yet
                elseif strcmp(MainInput.AnalysisType,'GasExchange')
                    LoadData.ismrmrd.Siemens.ute_to_ismrmrd(HFileName,H_name,fullfile(HDataLocation,[H_name '.h5']));
                    MainInput.HFileName = [H_name '.h5'];         
                    [GasExchange] = LoadData.ismrmrd.radial_3D_XeCTC_H_recon(MainInput,GasExchange);
                end
        elseif strcmp(MainInput.HDataext,'.7') && strcmp(MainInput.Scanner,'GE') 
                    HFullPath = MainInput.HFullPath;
                    HDataLocation = MainInput.XeDataLocation;
                    HFileName = MainInput.HFileName;
                    H_name = MainInput.H_name;
                    H_ext = MainInput.HDataext;
                    cd(HDataLocation) 
                if strcmp(MainInput.AnalysisType,'Ventilation')    
                    MainInput.ReconImageMode = 'proton';
                    LoadData.ismrmrd.GE.gre_to_ismrmrd(HFileName,H_name,fullfile(HDataLocation,[H_name '.h5']));
                    MainInput.HFileName = [H_name '.h5'];
                    [Image] = LoadData.ismrmrd.cartesian_2D_recon(MainInput);
                    Proton.Image  = Image;  
                elseif strcmp(MainInput.AnalysisType,'Diffusion')  
                % not supported yet
                elseif strcmp(MainInput.AnalysisType,'GasExchange')  
                    LoadData.ismrmrd.GE.ute_to_ismrmrd(HFileName,H_name,fullfile(HDataLocation,[H_name '.h5']));                    
                    MainInput.HFileName = [H_name '.h5'];         
                    [GasExchange] = LoadData.ismrmrd.radial_3D_XeCTC_H_recon(MainInput,GasExchange);
                end

    %--------------------- add new read load function here --------------------
                
    %     elseif strcmp(MainInput.HDataext,'add DataType') == 1
            %  add load/read function here 
        end 
    catch
        disp('no proton images is selected')
        MainInput.NoProtonImage = 1;
    end
end % end of Load/Read Proton data 

close all;
end

% Please add updates here
% 
% 
% 
% 
