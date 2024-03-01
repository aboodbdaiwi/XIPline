function [Dis_Fid,Gas_Fid,Dis_Traj,Gas_Traj,Params,Post_Cal] = pull_dis_data(Xe_file)

% Function to generate dissolved data nicely based on the methods used to
% date at KUMC. If these don't work for your purposes, you can create your
% own function to generate data. Outputs are:

%Dis_Fid - dissolved FIDs (in a matrix that is nPoints x NProjections)
%Gas_Fid - gas FIDs (in a matrix that is nPoints x NProjections)
%Dis_Traj - dissolved trajectories (in a matrix that is 3 x nPoints x NProjections)
%Gas_Traj - gas trajectories (in a matrix that is 3 x nPoints x NProjections)
%Post_Cal - Raw data for a post-imaging calibration (i.e. spectroscopy appended to the acquisition) - Use nan if not collected.
%Params - Structure with fields:
%       scandatestr (Scan date)
%       TR (Repetition Time)
%       GasFA (Flip angle for gas)
%       DisFA (Flip angle for dissolved)
%       TE (Echo Time)
%       GE_FOV (scalar Field of View in mm)
%       imsize (scalar Matrix size)
%       Dwell (dwell time (in s))
%       Ramp_Time (gradient ramp time (usually 100))
%       freq (center frequency in Hz)
%       freq_offset (offset of dissolved from gas)

parent_path = which('DataImport.pull_dis_data');
idcs = strfind(parent_path,filesep);%determine location of file separators
parent_path = parent_path(1:idcs(end-1)-1);%remove file


Xe_Dat_twix = DataImport.mapVBVD(Xe_file,'ignoreSeg');
Xe_Dat_twix.flagIgnoreSeg = 1;
Xe_Dat_twix.image.flagIgnoreSeg = 1;
Xe_Dat_twix.image.flagAverageReps = 1;
Xe_Dat_twix.flagAverageReps = 1;

%Identify sequence:
Seq_Name = Xe_Dat_twix.hdr.Config.SequenceFileName;

% Handle Mugler and Niedbalski sequences:
if contains(Seq_Name,'xe_radial_Dixon')
    if str2double(Seq_Name((end-3):end)) >= 2102 && Xe_Dat_twix.image.NPar > 1
        Ver = 'Mugler_spec';
    else
        Ver = 'Mugler_basic';
    end
elseif contains(Seq_Name,'SPIRAL_Rad_DIXON')
    if contains(Seq_Name,'20220627')
        Ver = 'Niedbalski_spec';
    else
        Ver = 'Niedbalski_basic';
    end
else
    error('Unable to identify Dissolved Phase Sequence Used');
end

%% Read in data based on version
Post_Cal = nan;
%For my parameter settings, this is a good way to do things... but probably
%not perfectly generic.
switch Ver
    case 'Mugler_spec'
            %Issues with memory using the twix object. Use this one instead
        Xe_Raw = DataImport.ReadSiemensMeasVD13_idea(Xe_file);
        fid = Xe_Raw.rawdata;
        loop = Xe_Raw.loopcounters;
        %Relevant Columns are 1, 5, and 7 - Reshape data to be in the shape
        %pts, projections,contrast(gas/dis),repetitions
        first_nrep = find(loop(:,7)==1,2); %index 1 will be first index of pre-spectra
        first_nrep = first_nrep(2); %this index will be first index of images
        last_rep = find(loop(:,7)==1,1,'last'); %this index will be first index of post-spectra
        Spec_Pre = zeros(max(loop(1:first_nrep,7)),loop(1,16));
        Im_Raw = zeros(max(loop(:,1)),loop(first_nrep,16),2);
        Spec_Post = zeros(size(loop,1)-last_rep+1,loop(last_rep,16));
        for i = 1:size(loop,1)
            %get pre-spectroscopy data
            if loop(i,4) == 1 && loop(i,29) < 0.5*size(loop,1)
                Spec_Pre(loop(i,7),:) = fid(i,1:loop(1,16));
            elseif loop(i,4) == 33
                Im_Raw(loop(i,1),:,loop(i,5)) = fid(i,1:loop(first_nrep,16));
            else
                Spec_Post(loop(i,7),:) = fid(i,1:loop(last_rep,16));
            end
        end
        Xe_Raw = Im_Raw;
        Dis_Fid = transpose(Xe_Raw(:,:,2));
        Gas_Fid = transpose(Xe_Raw(:,:,1));
        Post_Cal = Spec_Post;
        Gas_Traj = DataImport.gas_exchange_traj_gen(size(Gas_Fid,1),size(Gas_Fid,2),Xe_Dat_twix);
        Dis_Traj = DataImport.gas_exchange_traj_gen(size(Dis_Fid,1),size(Dis_Fid,2),Xe_Dat_twix);
        Params.imsize = Xe_Dat_twix.hdr.MeasYaps.sKSpace.lBaseResolution;
        Params.TR = (Xe_Dat_twix.hdr.MeasYaps.alTR{1}/1000)*2;
        Params.TE = (Xe_Dat_twix.hdr.MeasYaps.alTE{1}/1000);
        Params.GasFA = Xe_Dat_twix.hdr.MeasYaps.adFlipAngleDegree{2};
        Params.DisFA = Xe_Dat_twix.hdr.MeasYaps.adFlipAngleDegree{1};
        Params.freq_offset = Xe_Dat_twix.hdr.MeasYaps.sWipMemBlock.alFree{17};
        Params.freq = Xe_Dat_twix.hdr.Dicom.lFrequency;
        Params.Dwell = Xe_Dat_twix.hdr.MeasYaps.sRXSPEC.alDwellTime{1,1}*1e-9;
        scanDate = Xe_Dat_twix.hdr.Phoenix.tReferenceImage0; 
        scanDate = strsplit(scanDate,'.');
        scanDate = scanDate{end};
        scanDateStr = [scanDate(1:4),'-',scanDate(5:6),'-',scanDate(7:8)];
        Params.scandatestr = scanDateStr;
        Params.FOV = Xe_Dat_twix.hdr.Config.ReadFoV;
        Params.GE_FOV = 400;
        Params.GE_Voxel = 6.25;  
        Params.Cal_Dwell = Xe_Dat_twix.hdr.MeasYaps.sWipMemBlock.adFree{15}*1e-6/2;
        Params.Ramp_Time = Xe_Dat_twix.hdr.MeasYaps.sWipMemBlock.alFree{8};
    case 'Mugler_basic'
        
        Xe_Raw = double(squeeze(Xe_Dat_twix.image()));
        Dis_Fid = Xe_Raw(:,:,2);
        Gas_Fid = Xe_Raw(:,:,1);
        Gas_Traj = DataImport.gas_exchange_traj_gen(size(Xe_Raw,1),size(Xe_Raw,2),Xe_Dat_twix);
        Dis_Traj = DataImport.gas_exchange_traj_gen(size(Xe_Raw,1),size(Xe_Raw,2),Xe_Dat_twix);
        
        Params.imsize = Xe_Dat_twix.hdr.MeasYaps.sKSpace.lBaseResolution;
        Params.TR = (Xe_Dat_twix.hdr.MeasYaps.alTR{1}/1000)*2;
        Params.TE = (Xe_Dat_twix.hdr.MeasYaps.alTE{1}/1000);
        Params.GasFA = Xe_Dat_twix.hdr.MeasYaps.adFlipAngleDegree{2};
        Params.DisFA = Xe_Dat_twix.hdr.MeasYaps.adFlipAngleDegree{1};
        Params.freq_offset = Xe_Dat_twix.hdr.MeasYaps.sWipMemBlock.alFree{17};
        Params.freq = Xe_Dat_twix.hdr.Dicom.lFrequency;
        Params.Dwell = Xe_Dat_twix.hdr.MeasYaps.sRXSPEC.alDwellTime{1,1}*1e-9;
        scanDate = Xe_Dat_twix.hdr.Phoenix.tReferenceImage0; 
        scanDate = strsplit(scanDate,'.');
        scanDate = scanDate{end};
        scanDateStr = [scanDate(1:4),'-',scanDate(5:6),'-',scanDate(7:8)];
        Params.scandatestr = scanDateStr;
        Params.FOV = Xe_Dat_twix.hdr.Config.ReadFoV;
        Params.GE_FOV = 400;
        Params.GE_Voxel = 6.25;    
        Params.Ramp_Time = Xe_Dat_twix.hdr.MeasYaps.sWipMemBlock.alFree{9};
    case 'Niedbalski_basic'
        Xe_Raw = DataImport.ReadSiemensMeasVD13_idea(Xe_file);
        fid = Xe_Raw.rawdata;
        Dis_Fid = fid(1:2:(end-1),:);
        Dis_Fid(:,65:end) = [];
        Dis_Fid = Dis_Fid';
        Gas_Fid = fid(2:2:end,:);
        Gas_Fid = Gas_Fid';
        Dis_Traj = DataImport.get_allinone_distraj(Dis_Fid);
        traj_file = fullfile(parent_path,'\Traj\Vent_GasExchange_20210819_Traj.dat');
        traj_twix = DataImport.mapVBVD(traj_file);

        Gas_Traj = DataImport.spiral_coords_from_dat(traj_twix,Xe_Dat_twix);
        
        Params.imsize = Xe_Dat_twix.hdr.MeasYaps.sKSpace.lBaseResolution;
        Params.TR = ((Xe_Dat_twix.hdr.MeasYaps.alTR{1}+Xe_Dat_twix.hdr.MeasYaps.alTR{2})/1000);
        Params.TE = (Xe_Dat_twix.hdr.MeasYaps.alTE{1}/1000);
        Params.GasFA = Xe_Dat_twix.hdr.MeasYaps.adFlipAngleDegree{2};
        Params.DisFA = Xe_Dat_twix.hdr.MeasYaps.adFlipAngleDegree{1};
        Params.freq_offset = Xe_Dat_twix.hdr.MeasYaps.sWipMemBlock.alFree{5};
        Params.freq = Xe_Dat_twix.hdr.Dicom.lFrequency;
        Params.Dwell = Xe_Dat_twix.hdr.MeasYaps.sRXSPEC.alDwellTime{1,1}*1e-9;
        scanDate = Xe_Dat_twix.hdr.Phoenix.tReferenceImage0; 
        scanDate = strsplit(scanDate,'.');
        scanDate = scanDate{end};
        scanDateStr = [scanDate(1:4),'-',scanDate(5:6),'-',scanDate(7:8)];
        Params.scandatestr = scanDateStr;
        Params.GE_FOV = Xe_Dat_twix.hdr.Config.ReadFoV;
        Params.scale_factor = Xe_Dat_twix.hdr.MeasYaps.sWipMemBlock.adFree{6};
        Params.GE_Voxel = Params.GE_FOV/Params.imsize;
        Params.Vent_Voxel = Params.GE_FOV/(Params.imsize*Params.scale_factor);
        [Gas_Fid,Gas_Traj] = DataImport.truncate_spirals(Gas_Fid,Gas_Traj,Params.scale_factor);
        Params.Ramp_Time = 100;
    case 'Niedbalski_spec'
        Xe_Raw = DataImport.ReadSiemensMeasVD13_idea(Xe_file);
        fid = Xe_Raw.rawdata;
        Dis_Fid = fid(1:2:(end-1),:);
        %Again could do something more elegant, but start with easy
        Dis_Fid(:,65:end) = [];
        Dis_Fid = Dis_Fid';
        Gas_Fid = fid(2:2:end,:);
        Gas_Fid = Gas_Fid';
        Post_Cal = fid(end,1:512);
        Dis_Traj = DataImport.get_allinone_distraj(Dis_Fid);
        %Need to fix this path so it's more generic!
        traj_file = fullfile(parent_path,'\Traj\Vent_GasExchange_20210819_Traj.dat');
        traj_twix = DataImport.mapVBVD(traj_file);

        Gas_Traj = DataImport.spiral_coords_from_dat(traj_twix,Xe_Dat_twix);
        Params.imsize = Xe_Dat_twix.hdr.MeasYaps.sKSpace.lBaseResolution;
        Params.TR = ((Xe_Dat_twix.hdr.MeasYaps.alTR{1}+Xe_Dat_twix.hdr.MeasYaps.alTR{2})/1000);
        Params.TE = (Xe_Dat_twix.hdr.MeasYaps.alTE{1}/1000);
        Params.GasFA = Xe_Dat_twix.hdr.MeasYaps.adFlipAngleDegree{2};
        Params.DisFA = Xe_Dat_twix.hdr.MeasYaps.adFlipAngleDegree{1};
        Params.freq_offset = Xe_Dat_twix.hdr.MeasYaps.sWipMemBlock.alFree{5};
        Params.freq = Xe_Dat_twix.hdr.Dicom.lFrequency;
        Params.Dwell = Xe_Dat_twix.hdr.MeasYaps.sRXSPEC.alDwellTime{1,1}*1e-9;
        Params.Cal_Dwell = Params.Dwell;
        scanDate = Xe_Dat_twix.hdr.Phoenix.tReferenceImage0; 
        scanDate = strsplit(scanDate,'.');
        scanDate = scanDate{end};
        scanDateStr = [scanDate(1:4),'-',scanDate(5:6),'-',scanDate(7:8)];
        Params.scandatestr = scanDateStr;
        Params.GE_FOV = Xe_Dat_twix.hdr.Config.ReadFoV;
        Params.scale_factor = Xe_Dat_twix.hdr.MeasYaps.sWipMemBlock.adFree{6};
        Params.GE_Voxel = Params.GE_FOV/Params.imsize;
        Params.Vent_Voxel = Params.GE_FOV/(Params.imsize*Params.scale_factor);
        
        [Gas_Fid,Gas_Traj] = DataImport.truncate_spirals(Gas_Fid,Gas_Traj,Params.scale_factor);
        Params.Ramp_Time = 100;
end