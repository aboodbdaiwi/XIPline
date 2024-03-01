function [MultiRAID,pathname] = ReadSiemensMeasVD13_idea(measfile, dispopt)
%function MultiRAID = ReadSiemensMeasVD13_idea(measfile, dispopt)
%
% Reads the VD13 measurement data (.dat) got from the Siemens scanners RAID file system using
% TWIX
%
% dispopt, - 'on' or 'off'. The default is 'off'. The FFTed images for slices of the first 3D
% is displayed only for your reference. if your data were segmented, PF, ReadOutOffcentre, or
% arranged in the special way, you might not see the correct images. You have to turn this
% option off, and work on the rawdata output from this function for the right raw line order.
%
% Output MultiRAID is an array of structures each cotaining the following
% fields corresponding to one RAID file:
% rawdata, - raw data in matrix; oversmapled (x2) in complex (A+iB), indexed by loopcounters and
% ulChannelId;
% loopcounters,- all the sMDH counters and ChannelId corresponding to all the ADC data lines;
% All dimension are in the same order as saved in sMDH.sLC;
% sMDH, - the MDH of the 1st ADC acquision except that extended sMDH.sLC is set to maximum
% loop values.
% MrProt, - Some, not all information in the header file. 
% Yaps, - Not effectively used yet.
% lc_names, - tells you what are in the Loopcounters, e.g., 
% sCHA, - 1st channel structure
% RaidInfo, - information of one RAID file
%
% EXAMPLE:
% 
%     MultiRAD = ReadSiemensMeasVD13_idea('Z:\n4\pkg\MeasurementData\FLASH\meas.dat');
%
%     for iRAID = 1:size(MultiRAID, 2) 
%         
%         rawdata         = MultiRAID(iRAID).rawdata;
%         loopcounters    = MultiRAID(iRAID).loopcounters;
%         sMDH            = MultiRAID(iRAID).sMDH;
%         MrProt          = MultiRAID(iRAID).MrProt;
%         Yaps            = MultiRAID(iRAID).Yaps;
%         lc_names        = MultiRAID(iRAID).lc_names;
%         sCHA            = MultiRAID(iRAID).sCHA;
%         RaidInfo        = MultiRAID(iRAID).RaidInfo;
%         
%         protname = deblank(RaidInfo.protname');
%         
%         % ......
%         
%     end
%
% Maolin Qiu YALE 06-25-2015 (maolin.qiu(at)yale(dot)edu)

if nargin < 1 | ~exist(measfile)
%     help ReadSiemensMeasVD13_idea;
%     [filename pathname] = uigetfile( ...
%        {'*.dat';'*.out';'*.*'}, ...
%         'Pick a Siemens MEASUREMENT file');
%     if ~filename & ~pathname
%         disp(['You selected no file.']);
%         return;
    error('Measurement File Not Found')
%     else
%         measfile = fullfile(pathname, filename);
%     end
else
    pathname = fileparts(measfile);
end

if nargin < 2
    dispopt = 'off';
end

% measfile

[pa na ex] = fileparts(measfile);
measfile = fullfile(pa, [na ex]);
if exist(measfile, 'file')
    disp(['Measurement data file: ' measfile]);
else
    disp(['Measurement data file does not exist: ' measfile]);
    return;
end

%% init output

rawdata = []; data = []; dimensions = []; loopcounters = []; sMDH = []; MrProt = []; Yaps = []; lc_names = []; sCHA = []; %MultiRAID = [];

MAX_NUMBER_OF_CHANNELS = 1024; % if a channel id is set greater than this number in the meas data it will not be considered as a channel id

%% RAID structure

% for i = 1:NRAID % read multi-RAID info 
%     RaidInfo(i).measID      = fread(fid,1,'uint32'); % 3:5450
%     RaidInfo(i).fileID      = fread(fid,1,'uint32'); % 4:55097
%     RaidInfo(i).measOffset  = fread(fid,1,'uint64'); % 5+6:10240
%     RaidInfo(i).measLength  = fread(fid,1,'uint64'); % 7+8:8628288
%     RaidInfo(i).patname     = fread(fid,64,'uint8=>char');
%     RaidInfo(i).protname    = fread(fid,64,'uint8=>char');
% end % NRAID


aRaidInfo = struct(...
    'measID', [], ...
    'fileID', [], ...
    'measOffset', [], ...
    'measLength', [], ...
    'patname', [], ...
    'protname', [] ...
    );

%     MultiRAID(iRAID).rawdata = rawdata;
%     MultiRAID(iRAID).loopcounters = loopcounters;
%     MultiRAID(iRAID).sMDH = sMDH;
%     MultiRAID(iRAID).MrProt = MrProt;
%     MultiRAID(iRAID).Yaps = Yaps;
%     MultiRAID(iRAID).lc_names = lc_names;
%     MultiRAID(iRAID).sCHA = sCHA;
%     MultiRAID(iRAID).RaidInfo = RaidInfo(iRAID);


aRAID = struct(...
    'rawdata', [], ...
    'loopcounters', [], ...
    'sMDH', [], ...
    'MrProt', [], ...
    'Yaps', [], ...
    'lc_names', [], ...
    'sCHA', [], ...
    'RaidInfo', [] ...
    );

%% Constants used in sMDH

MDH_NUMBEROFEVALINFOMASK   = 2;
MDH_NUMBEROFICEPROGRAMPARA = 24;
MDH_RESERVEDHDRPARA        = 4;
% MDH_FREEHDRPARA          = 4;

IDX_DIM                    = 64;

%%--------------------------------------------------------------------------%%
%% Definition of loop counter structure                                     %%
%% Note: any changes of this structure affect the corresponding swapping    %%
%%       method of the measurement data header proxy class (MdhProxy)       %%
%%--------------------------------------------------------------------------%%

sLoopCounter = struct( ...
  'ushLine',0,...                  %% unsigned short  line index                   %%
  'ushAcquisition',0,...           %% unsigned short  acquisition index            %%
  'ushSlice',0,...                 %% unsigned short  slice index                  %%
  'ushPartition',0,...             %% unsigned short  partition index              %%
  'ushEcho',0,...                  %% unsigned short  echo index                   %%
  'ushPhase',0,...                 %% unsigned short  phase index                  %%
  'ushRepetition',0,...            %% unsigned short  measurement repeat index     %%
  'ushSet',0,...                   %% unsigned short  set index                    %%
  'ushSeg',0,...                   %% unsigned short  segment index  (for TSE)     %%
  'ushIda',0,...                   %% unsigned short  IceDimension a index         %%
  'ushIdb',0,...                   %% unsigned short  IceDimension b index         %%
  'ushIdc',0,...                   %% unsigned short  IceDimension c index         %%
  'ushIdd',0,...                   %% unsigned short  IceDimension d index         %%
  'ushIde',0 ...                   %% unsigned short  IceDimension e index         %%
);                                 %% sizeof : 28 byte             %%

%%--------------------------------------------------------------------------%%
%%  Definition of slice vectors                                             %%
%%--------------------------------------------------------------------------%%

sVector = struct( ...
  'flSag',0.0,...       %% float
  'flCor',0.0,...       %% float
  'flTra',0.0 ...       %% float
);

sSliceData = struct( ...
  'sSlicePosVec',sVector,...                   %% slice position vector               %%
  'aflQuaternion',zeros(1,4) ...               %% float rotation matrix as quaternion %%
);                                              %% sizeof : 28 byte                    %%

%%--------------------------------------------------------------------------%%
%%  Definition of cut-off data                                              %%
%%--------------------------------------------------------------------------%%

sCutOffData = struct( ...
  'ushPre',0,...               %% unsigned short  write ushPre zeros at line start %%
  'ushPost',0 ...              %% unsigned short  write ushPost zeros at line end  %%
); %% 4 bytes

%%--------------------------------------------------------------------------%%
%% Definition of Channel data for VD13                                      %%
%%--------------------------------------------------------------------------%%


sChannelHeader = struct( ...
  ... %'ulTypeAndChannelLength',0,...   %uint32_t,         ///< 0: (4) bit  0.. 7: type (0x02 => ChannelHeader)
  ...                                   %                  ///<        bit  8..31: channel length (header+data) in byte
  ...                                   %                  ///<        type   := ulTypeAndChannelLength & 0x000000FF
  ...                                   %                  ///<        length := ulTypeAndChannelLength >> 8
  'ulType',0,...                        %= fread(fid, 1, 'ubit8'); 
  'ulChannelLength',0,...               %= fread(fid, 1, 'ubit24'); 
  'lMeasUID',0,...                      %int32_t,          ///< 4: (4) measurement user ID
  'ulScanCounter',0,...                 %uint32_t,         ///< 8: (4) scan counter [1...]
  'ulReserved1',0,...                   %uint32_t,         ///< 12:(4) reserved
  'ulSequenceTime',0,...                %uint32_t,         ///< 16:(4) Sequence readout starting time bit 31..9 time in [10us]
  ...                                   %                  ///<                                       bit  8..0 time in [25ns]
  'ulUnused2',0,...                     %uint32_t,         ///< 20:(4) unused
  'ulChannelId',0,...                   %uint16_t,         ///< 24:(2) unused
  'ulUnused3',0,...                     %uint16_t,         ///< 26:(2) unused
  'ulCRC',0 ...                         %uint32_t,         ///< 28:(4) CRC32 checksum of channel header
);   %% // total length:  32 byte


lChannelHeaderLength = 32; % bytes; This is new for VD13


%%--------------------------------------------------------------------------%%
%%  Definition of measurement data header                                   %%
%%--------------------------------------------------------------------------%%

sMDH = struct( ...
  ... % 'ulFlagsAndDMALength',0,...       %uint32_t,                        ///<  0: ( 4) bit  0..24: DMA length [bytes]
  ...                               %                                       ///<          bit     25: pack bit
  ...                               %                                       ///<          bit 26..31: pci_rx enable flags
  'ulDMALength',0,...               %= fread(fid, 1, 'ubit25');          
  'ulPackBit',0,...                 %= fread(fid, 1, 'ubit1');
  'ulPCI_rx',0,...                  %= fread(fid, 1, 'ubit6');           
  'lMeasUID',0,...                  %int32_t,                               ///<  4: ( 4) measurement user ID
  'ulScanCounter',0,...             %uint32_t,                              ///<  8: ( 4) scan counter [1...]
  'ulTimeStamp',0,...               %uint32_t,                              ///< 12: ( 4) time stamp [2.5 ms ticks since 00:00]
  'ulPMUTimeStamp',0,...            %uint32_t,                              ///< 16: ( 4) PMU time stamp [2.5 ms ticks since last trigger]
  'ushSystemType',0,...             %uint16_t,                              ///< 20: ( 2) System type (todo: values?? ####)
  'ulPTABPosDelay',0,...            %uint16_t,                              ///< 22: ( 2) PTAb delay ??? TODO: How do we handle this ####
  'lPTABPosX',0,...                 %int32_t,                               ///< 24: ( 4) absolute PTAB position in [µm]
  'lPTABPosY',0,...                 %int32_t,                               ///< 28: ( 4) absolute PTAB position in [µm]
  'lPTABPosZ',0,...                 %int32_t,                               ///< 32: ( 4) absolute PTAB position in [µm]
  'ulReserved1',0,...               %uint32_t,                              ///< 36: ( 4) reserved for future hardware signals
  'aulEvalInfoMask',zeros(1,MDH_NUMBEROFEVALINFOMASK),...       %uint32_t,  ///< 40: ( 8) evaluation info mask field
  'ushSamplesInScan',0,...          %uint16_t,                              ///< 48: ( 2) # of samples acquired in scan
  'ushUsedChannels',0,...           %uint16_t,                              ///< 50: ( 2) # of channels used in scan
  'sLC',sLoopCounter,...            %                                       ///< 52: (28) loop counters
  'sCutOff',sCutOffData,...         %                                       ///< 80: ( 4) cut-off values
  'ushKSpaceCentreColumn',0,...     %uint16_t,                              ///< 84: ( 2) centre of echo
  'ushCoilSelect',0,...             %uint16_t,                              ///< 86: ( 2) Bit 0..3: CoilSelect
  'fReadOutOffcentre',0.0,...       %float,                                 ///< 88: ( 4) ReadOut offcenter value
  'ulTimeSinceLastRF',0,...         %uint32_t,                              ///< 92: ( 4) Sequence time stamp since last RF pulse
  'ushKSpaceCentreLineNo',0,...     %uint16_t,                              ///< 96: ( 2) number of K-space centre line
  'ushKSpaceCentrePartitionNo',0,...%uint16_t,                              ///< 98: ( 2) number of K-space centre partition
  'sSD',sSliceData,...              %                                       ///< 100:(28) Slice Data
  'aushIceProgramPara',zeros(1,MDH_NUMBEROFICEPROGRAMPARA),...   %uint16_t, ///< 128:(48) free parameter for IceProgram
  'aushReservedPara',zeros(1,MDH_RESERVEDHDRPARA),...            %uint16_t, ///< 176:( 8) unused parameter (padding to next 192byte alignment )
  ...                               %                                       ///<          NOTE: These parameters MUST NOT be used by any application (for future use)
  'ushApplicationCounter',0,...     %uint16_t,                              ///< 184 ( 2)
  'ushApplicationMask',0,...        %uint16_t,                          /   ///< 186 ( 2)
  'ulCRC',0 ...                     %uint32_t,                              ///< 188:( 4) CRC 32 checksum
); % New MDH for VD13    

lScanHeaderLength = 192; % bytes; in VB this is 128 as below:

%   'ulDMALength',0,...                                       %% unsigned long  DMA length [bytes] must be                        4 bytes %% first parameter
%   'lMeasUID',0,...                                          %% long           measurement user ID                               4
%   'ulScanCounter',0,...                                     %% unsigned long  scan counter [1...]                               4
%   'ulTimeStamp',0,...                                       %% unsigned long  time stamp [2.5 ms ticks since 00:00]             4
%   'ulPMUTimeStamp',0,...                                    %% unsigned long  PMU time stamp [2.5 ms ticks since last trigger]  4
%   'aulEvalInfoMask',zeros(1,MDH_NUMBEROFEVALINFOMASK),...   %% unsigned long  evaluation info mask field                        8
%   'ushSamplesInScan',0,...                                  %% unsigned short # of samples acquired in scan                     2
%   'ushUsedChannels',0,...                                   %% unsigned short # of channels used in scan                        2   =32
%   'sLC',sLoopCounter,...                                    %% loop counters                                                    28  =60
%   'sCutOff',sCutOffData,...                                 %% cut-off values                                                   4
%   'ushKSpaceCentreColumn',0,...                             %% unsigned short centre of echo                                    2
%   'ushDummy',0,...                                          %% unsigned short for swapping                                      2
%   'fReadOutOffcentre',0.0,...                               %% float          ReadOut offcenter value                           4
%   'ulTimeSinceLastRF',0,...                                 %% unsigned long  Sequence time stamp since last RF pulse           4
%   'ushKSpaceCentreLineNo',0,...                             %% unsigned short number of K-space centre line                     2
%   'ushKSpaceCentrePartitionNo',0,...                        %% unsigned short number of K-space centre partition                2
%   'aushIceProgramPara',zeros(1,MDH_NUMBEROFICEPROGRAMPARA),... %% unsigned short free parameter for IceProgram                  8  =88
%   'aushFreePara',zeros(1,MDH_FREEHDRPARA),...               %% unsigned short free parameter                                    4 * 2 =   8
%   'sSD',sSliceData,...                                      %% Slice Data                                                       28 =124
%   'ulChannelId',0 ...                                       %% unsigned long	 channel Id must be the last parameter            4
% );                                                          %% total length: 32 * 32 Bit (128 Byte)                             128
%% MDH_H %%

%% read header information
tic;
[pa na ex] = fileparts(measfile);
diary_file = fullfile(pa, [na '_diary.txt']);
if exist(diary_file,'file')
    delete(diary_file);
end


diary(diary_file);

fid = fopen(measfile,'r');

% first 8 numbers in the beginning of .dat
% 1           0
% 2           1
% 3        5450
% 4       55097
% 5       10240
% 6           0
% 7     8628288
% 8           0

disp(['This code only reads VD measurement data.']);

ID      = fread(fid,1,'uint32'); % 1:0
NRAID   = fread(fid,1,'uint32'); % 2:1 for VD  max = 64

for i = 1:NRAID % read multi-RAID info 
    aRaidInfo.measID      = fread(fid,1,'uint32'); % 3:5450
    aRaidInfo.fileID      = fread(fid,1,'uint32'); % 4:55097
    aRaidInfo.measOffset  = fread(fid,1,'uint64'); % 5+6:10240
    aRaidInfo.measLength  = fread(fid,1,'uint64'); % 7+8:8628288
    aRaidInfo.patname     = fread(fid,64,'uint8=>char');
    aRaidInfo.protname    = fread(fid,64,'uint8=>char');
    RaidInfo(i) = aRaidInfo;
end % NRAID

disp(['Read the Measurement Data file. This will take long if it is HUGE ... ...']);
disp(['Total Number of RAID Files in this Measurement Data : ' num2str(NRAID)]);

for iRAID = 1:NRAID,


    measOffset = RaidInfo(iRAID).measOffset;
    measLength = RaidInfo(iRAID).measLength;
    fseek(fid,measOffset,'bof'); % this is where the header begins
    hdr_len_RAID  = fread(fid,1,'uint32');
    readout_start = measOffset + hdr_len_RAID; % here is where the read out measurement data start as a scan header
    readout_end   = measOffset + hdr_len_RAID + measLength;
    protname = deblank(RaidInfo(iRAID).protname');

    disp('-------------------------------------------------------')
    disp('-------------------------------------------------------')
    disp(['RAID File # : ' num2str(iRAID) ' of total ' num2str(NRAID)]);
    disp(['Length of the Header (' protname ') : ' num2str(hdr_len_RAID)]);
    disp('-------------------------------------------------------')
    disp('-------------------------------------------------------')

    % save the header for this RAID file
    fseek(fid, measOffset, 'bof'); % rewind to read and save a header file for later use
    a_header = fread(fid, hdr_len_RAID);

    measfile_header = fullfile(pa, [na '_' num2str(iRAID) '_' protname '.header']);

    %PJN edit to accomodate dumb gas exchange name file
    badname = strfind(measfile_header,'<');
    measfile_header(badname) = '_';
    fid_hdr = fopen(measfile_header,'w');
    fwrite(fid_hdr, a_header);
    fclose(fid_hdr);
    disp(['Common Header of this RAID file is saved in : ' measfile_header]);
    
    % init for this RAID file
    rawdata = []; data = []; dimensions = []; loopcounters = []; MrProt = []; Yaps = []; lc_names = []; sCHA = [];

    % find the total size of the file
    [tmp Nbytes] = fread(fid,inf,'uchar'); tmp = []; 
    FileSize = ftell(fid); 
    disp(['Original Measurement Data File Size : ' num2str(FileSize)]);

    next_fid_pos = readout_start;

%   cur_pos = readout_start;

    lcounter = 1;

    LastADCData = 0; KMAX = 0;

    while 1

        fseek(fid, next_fid_pos, 'bof'); %%..maolin modified according to VD data 5/5/2015

        if feof(fid), break; end

        % read the 192 byte miniHeader for each read-out
        sMDH.ulDMALength                                = fread(fid, 1, 'ubit25');          % 0:    ( 4) bytes
        sMDH.ulPackBit                                  = fread(fid, 1, 'ubit1');
        sMDH.ulPCI_rx                                   = fread(fid, 1, 'ubit6');           
        sMDH.lMeasUID                                   = fread(fid, 1, 'int32');           % 4:    ( 4)
        sMDH.ulScanCounter                              = fread(fid, 1, 'uint32');          % 8:    ( 4)
        sMDH.ulTimeStamp                                = fread(fid, 1, 'uint32');          % 12:   ( 4)
        sMDH.ulPMUTimeStamp                             = fread(fid, 1, 'uint32');          % 16:   ( 4)
        sMDH.ushSystemType                              = fread(fid, 1, 'uint16');          % 20:   ( 2)
        sMDH.ulPTABPosDelay = fread(fid, 1, 'uint16');          %uint16_t,               ///< 22:   ( 2) PTAb delay ??? TODO: How do we handle this ####
        sMDH.lPTABPosX      = fread(fid, 1, 'int32');           %int32_t,                ///< 24:   ( 4) absolute PTAB position in [µm]
        sMDH.lPTABPosY      = fread(fid, 1, 'int32');           %int32_t,                ///< 28:   ( 4) absolute PTAB position in [µm]
        sMDH.lPTABPosZ      = fread(fid, 1, 'int32');           %int32_t,                ///< 32:   ( 4) absolute PTAB position in [µm]
        sMDH.ulReserved1    = fread(fid, 1, 'uint32');          %uint32_t,               ///< 36:   ( 4) reserved for future hardware signals

        for i = 1:MDH_NUMBEROFEVALINFOMASK       % 2
            sMDH.aulEvalInfoMask(i)                     = fread(fid, 1, 'uint32');          % 40:   ( 8 = 4 x 2)
        end

        sMDH.ushSamplesInScan                           = fread(fid, 1, 'uint16');          % 48:   ( 2)
        sMDH.ushUsedChannels                            = fread(fid, 1, 'uint16');          % 50:   ( 2)

        K = sMDH.ushUsedChannels;
        if K > MAX_NUMBER_OF_CHANNELS, K = 0; end % some number are set not as CHA ID so I added this to control, max number of channels <= 1024
        KMAX = max(K,KMAX);

        dummy                          = fread(fid,      14, 'uint16');                     % 52: 	(28 = 14 x 2)

        sMDH.sLC.ushLine                             = dummy(1)+1;
        sMDH.sLC.ushAcquisition                      = dummy(2)+1;
        sMDH.sLC.ushSlice                            = dummy(3)+1;
        sMDH.sLC.ushPartition                        = dummy(4)+1;
        sMDH.sLC.ushEcho                             = dummy(5)+1;
        sMDH.sLC.ushPhase                            = dummy(6)+1;
        sMDH.sLC.ushRepetition                       = dummy(7)+1;
        sMDH.sLC.ushSet                              = dummy(8)+1;
        sMDH.sLC.ushSeg                              = dummy(9)+1;
        sMDH.sLC.ushIda                              = dummy(10)+1;
        sMDH.sLC.ushIdb                              = dummy(11)+1;
        sMDH.sLC.ushIdc                              = dummy(12)+1;
        sMDH.sLC.ushIdd                              = dummy(13)+1;
        sMDH.sLC.ushIde                              = dummy(14)+1;                         

        sMDH.sCutOff.ushPre                          = fread(fid, 1, 'uint16');              % 80: 	( 2)
        sMDH.sCutOff.ushPost                         = fread(fid, 1, 'uint16');              % 82:	( 2)

        sMDH.ushKSpaceCentreColumn                   = fread(fid, 1, 'uint16');              % 84:   ( 2)
        sMDH.ushCoilSelect                           = fread(fid, 1, 'uint16');              % 86:   ( 2)
        sMDH.fReadOutOffcentre                       = fread(fid, 1, 'float');               % 88:   ( 4)
        sMDH.ulTimeSinceLastRF                       = fread(fid, 1, 'uint32');              % 92:   ( 4)
        sMDH.ushKSpaceCentreLineNo                   = fread(fid, 1, 'uint16');              % 96:   ( 2)
        sMDH.ushKSpaceCentrePartitionNo              = fread(fid, 1, 'uint16');              % 98:   ( 2)

        sMDH.sSD.sVector.flSag                       = fread(fid, 1, 'float');               % 100:  ( 4)
        sMDH.sSD.sVector.flCor                       = fread(fid, 1, 'float');               % 104:  ( 4)
        sMDH.sSD.sVector.flTra                       = fread(fid, 1, 'float');               % 108:  ( 4)

        for i = 1:4
            sMDH.sSD.aflQuaternion(i)                   = fread(fid, 1, 'float');           % 112:  (16 = 4 x 4)       
        end

        for i = 1:MDH_NUMBEROFICEPROGRAMPARA    % 24
            sMDH.aushIceProgramPara(i)              = fread(fid, 1, 'uint16');              % 128:  (48 = 24 x 2)
        end

        for i = 1:MDH_RESERVEDHDRPARA  % 4
            sMDH.aushReservedPara(i)                = fread(fid, 1, 'uint16');              % 176:  ( 8 = 4 x 2)
        end

        sMDH.ushApplicationCounter                   = fread(fid, 1, 'uint16');              % 184:  ( 2)
        sMDH.ushApplicationMask                      = fread(fid, 1, 'uint16');              % 186:  ( 2)
        sMDH.ulCRC                                   = fread(fid, 1, 'uint32');              % 188:  ( 4) 
                                                                                            % 192:  = total ok!!

        sMDH.ulChannelId                             = K;      %VD? no idea just set to the MAX!

        % Evaluation Information Mask
        eval_info_mask.MDH_ACQEND             = min(bitand(sMDH.aulEvalInfoMask(1), 2^0),1);
        eval_info_mask.MDH_RTFEEDBACK         = min(bitand(sMDH.aulEvalInfoMask(1), 2^1),1);
        eval_info_mask.MDH_HPFEEDBACK         = min(bitand(sMDH.aulEvalInfoMask(1), 2^2),1);
        eval_info_mask.MDH_SYNCDATA           = min(bitand(sMDH.aulEvalInfoMask(1), 2^5), 1);
        eval_info_mask.MDH_RAWDATACORRECTION  = min(bitand(sMDH.aulEvalInfoMask(1), 2^10),1);
        eval_info_mask.MDH_REFPHASESTABSCAN   = min(bitand(sMDH.aulEvalInfoMask(1), 2^14),1);
        eval_info_mask.MDH_PHASESTABSCAN      = min(bitand(sMDH.aulEvalInfoMask(1), 2^15),1);
        eval_info_mask.MDH_SIGNREV            = min(bitand(sMDH.aulEvalInfoMask(1), 2^17),1);
        eval_info_mask.MDH_PHASCOR            = min(bitand(sMDH.aulEvalInfoMask(1), 2^21),1);
        eval_info_mask.MDH_PATREFSCAN         = min(bitand(sMDH.aulEvalInfoMask(1), 2^22),1);
        eval_info_mask.MDH_PATREFANDIMASCAN   = min(bitand(sMDH.aulEvalInfoMask(1), 2^23),1);
        eval_info_mask.MDH_REFLECT            = min(bitand(sMDH.aulEvalInfoMask(1), 2^24),1);
        eval_info_mask.MDH_NOISEADJSCAN       = min(bitand(sMDH.aulEvalInfoMask(1), 2^25),1);
 
% following one miniHeader 3 things could happen: 
%           1) channel header + data for all channels
%           2) end of a RAID file
%           3) sychronization data

        % CASE 1:
        % when it's not sychronizing data, not end of a RAID file, there
        % are channel data following, so jump a fixed number of bytes
        if ~eval_info_mask.MDH_SYNCDATA && ~eval_info_mask.MDH_ACQEND && sMDH.ulDMALength~=0
            DMALen = sMDH.ulDMALength;
            sMDH.ulDMALength = lScanHeaderLength + (2*4*sMDH.ushSamplesInScan + lChannelHeaderLength) * sMDH.ushUsedChannels;
            if (DMALen ~= sMDH.ulDMALength)
                disp(['Readout: ' sprintf('%8.8d',lcounter) ' Modified sMDH.ulDMALength = ' num2str(sMDH.ulDMALength) ' Original = ' num2str(DMALen)]);
            end
        end % if this is TRUE the following 2 cases can not be TRUE so it actually moves to read channel header and channel data
        
        % CASE 2:
        % end of one RAID file and jump to the next
        % align to the next 512-byte block's start
        if eval_info_mask.MDH_ACQEND || sMDH.ulDMALength==0
            disp(['Readout: ' sprintf('%8.8d',lcounter) ' sMDH.MDH_ACQEND = ' num2str(eval_info_mask.MDH_ACQEND) ' sMDH.ulDMALength = ' num2str(sMDH.ulDMALength)]);
            next_fid_pos = next_fid_pos + sMDH.ulDMALength; % just to the end of this acq.
        % align to the next 512-byte block's start
            if mod(next_fid_pos,512)
                pre_next_fid_pos = next_fid_pos;
                next_fid_pos = next_fid_pos + 512 - mod(next_fid_pos,512);
                disp(['Readout: ' sprintf('%8.8d',lcounter) ' next_fid_pos = ' num2str(next_fid_pos) ' Original = ' num2str(pre_next_fid_pos) ' Adjust to next 512-byte block''s start by jumping forward Bytes = ' num2str(next_fid_pos-pre_next_fid_pos)]);
            end
            break; % skip the channel header and channel data reading, jump out of while loop and move to the next RAID file
        end

        % CASE 3:
        % skip SYNCDATA
        if eval_info_mask.MDH_SYNCDATA
            disp(['Readout: ' sprintf('%8.8d',lcounter) ' eval_info_mask.MDH_SYNCDATA = ' num2str(eval_info_mask.MDH_SYNCDATA) ' sMDH.ulDMALength = ' num2str(sMDH.ulDMALength) ' skipped']);
            next_fid_pos = next_fid_pos + sMDH.ulDMALength;
            continue; % skip the channel header and channel data reading, continue to the while looping (same RAID file)
        end


        % read channel header and channel data 
        
        lcha = 0;

        while lcha < K %% DATA = ScanHeader + NCha * (ChannelHeader + ADC(sMDH.ushSamplesInScan*2*4))


            sCHA.ulType           = fread(fid, 1, 'ubit8'); 
            sCHA.ulChannelLength  = fread(fid, 1, 'ubit24'); %uint32_t,         ///< 0: (4) bit  0.. 7: type (0x02 => ChannelHeader)
                                                             %                  ///<        bit  8..31: channel length (header+data) in byte
                                                             %                  ///<        type   := ulTypeAndChannelLength & 0x000000FF
                                                             %                  ///<        length := ulTypeAndChannelLength >> 8
            sCHA.lMeasUID         = fread(fid, 1, 'int32');  %int32_t,          ///< 4: (4) measurement user ID
            sCHA.ulScanCounter    = fread(fid, 1, 'uint32'); %uint32_t,         ///< 8: (4) scan counter [1...]
            sCHA.ulReserved1      = fread(fid, 1, 'uint32'); %uint32_t,         ///< 12:(4) reserved
            sCHA.ulSequenceTime   = fread(fid, 1, 'uint32'); %uint32_t,         ///< 16:(4) Sequence readout starting time bit 31..9 time in [10us]
                                                             %                  ///<                                       bit  8..0 time in [25ns]
            sCHA.ulUnused2        = fread(fid, 1, 'uint32'); %uint32_t,         ///< 20:(4) unused
            sCHA.ulChannelId      = fread(fid, 1, 'uint16'); %uint16_t,         ///< 24:(2) unused
            sCHA.ulUnused3        = fread(fid, 1, 'uint16'); %uint16_t,         ///< 26:(2) unused
            sCHA.ulCRC            = fread(fid, 1, 'uint32'); %uint32_t,         ///< 28:(4) CRC32 checksum of channel header
                                                         % total length:  32 byte

         if lcounter == 1
           sMDH1 = sMDH;
           sCHA1 = sCHA;
           aADC = -999999*ones(1,IDX_DIM+sMDH.ushSamplesInScan*2); % the first IDX_DIM entries are the indices of this ADC in terms of sLoopCounter
           aADC(1)  = sMDH.sLC.ushLine;
           aADC(2)  = sMDH.sLC.ushAcquisition;
           aADC(3)  = sMDH.sLC.ushSlice;
           aADC(4)  = sMDH.sLC.ushPartition;
           aADC(5)  = sMDH.sLC.ushEcho;
           aADC(6)  = sMDH.sLC.ushPhase;
           aADC(7)  = sMDH.sLC.ushRepetition;
           aADC(8)  = sMDH.sLC.ushSet;
           aADC(9)  = sMDH.sLC.ushSeg;
           aADC(10) = sMDH.sLC.ushIda;
           aADC(11) = sMDH.sLC.ushIdb;
           aADC(12) = sMDH.sLC.ushIdc;
           aADC(13) = sMDH.sLC.ushIdd;
           aADC(14) = sMDH.sLC.ushIde;
%          aADC(15) = mod(sCHA.ulChannelId, K)+1;
           aADC(15) = mod(lcha, K)+1; % simply looping cha id
           aADC(16) = sMDH.ushSamplesInScan;
           aADC(17) = sMDH.ushKSpaceCentreColumn;
           aADC(18) = sMDH.fReadOutOffcentre;
           aADC(19:19+MDH_NUMBEROFEVALINFOMASK-1) = sMDH.aulEvalInfoMask;
           aADC(23) = sMDH.ushKSpaceCentreLineNo;
           aADC(24) = sMDH.sCutOff.ushPre;
           aADC(25) = sMDH.sCutOff.ushPost;
           aADC(26) = sMDH.ushKSpaceCentrePartitionNo;
           aADC(27) = sMDH.ulDMALength;
           aADC(28) = sMDH.lMeasUID;
           aADC(29) = sMDH.ulScanCounter;
           aADC(30) = sMDH.ulTimeStamp;
           aADC(31) = sMDH.ulPMUTimeStamp;
           aADC(32) = sMDH.ulTimeSinceLastRF;
        %  aADC(35:35+MDH_NUMBEROFICEPROGRAMPARA-1) = sMDH.aushIceProgramPara;
           aADC(45:48) = sMDH.sSD.aflQuaternion;
           aADC(55:55+MDH_RESERVEDHDRPARA-1) = sMDH.aushReservedPara;
           aADC(60) = sMDH.sSD.sVector.flSag;
           aADC(61) = sMDH.sSD.sVector.flCor;
           aADC(62) = sMDH.sSD.sVector.flTra;
           aADC(64) = sMDH.ushCoilSelect;
           % ... ...
           aADC(IDX_DIM+1:end) = fread(fid, sMDH.ushSamplesInScan*2, 'float'); % 4 bytes in each float
           LADC = 128 + sMDH.ushSamplesInScan*2*4;   % in bytes now
           NADC = floor(Nbytes/LADC*2);
           data = zeros(NADC, size(aADC,2)); % Optimize to speed up!
           data(lcounter, :) = aADC;
       
           
         else
           aADC = -999999*ones(1,IDX_DIM+sMDH.ushSamplesInScan*2); % the first 14 entries are the indices of this ADC in terms of sLoopCounter
           aADC(1)  = sMDH.sLC.ushLine;
           aADC(2)  = sMDH.sLC.ushAcquisition;
           aADC(3)  = sMDH.sLC.ushSlice;
           aADC(4)  = sMDH.sLC.ushPartition;
           aADC(5)  = sMDH.sLC.ushEcho;
           aADC(6)  = sMDH.sLC.ushPhase;
           aADC(7)  = sMDH.sLC.ushRepetition;
           aADC(8)  = sMDH.sLC.ushSet;
           aADC(9)  = sMDH.sLC.ushSeg;
           aADC(10) = sMDH.sLC.ushIda;
           aADC(11) = sMDH.sLC.ushIdb;
           aADC(12) = sMDH.sLC.ushIdc;
           aADC(13) = sMDH.sLC.ushIdd;
           aADC(14) = sMDH.sLC.ushIde;
%          aADC(15) = mod(sCHA.ulChannelId, K)+1;
           aADC(15) = mod(lcha, K)+1; % simply looping cha id
           aADC(16) = sMDH.ushSamplesInScan;
           aADC(17) = sMDH.ushKSpaceCentreColumn;
           aADC(18) = sMDH.fReadOutOffcentre;
           aADC(19:19+MDH_NUMBEROFEVALINFOMASK-1) = sMDH.aulEvalInfoMask;
           aADC(23) = sMDH.ushKSpaceCentreLineNo;
           aADC(24) = sMDH.sCutOff.ushPre;
           aADC(25) = sMDH.sCutOff.ushPost;
           aADC(26) = sMDH.ushKSpaceCentrePartitionNo;
           aADC(27) = sMDH.ulDMALength;
           aADC(28) = sMDH.lMeasUID;
           aADC(29) = sMDH.ulScanCounter;
           aADC(30) = sMDH.ulTimeStamp;
           aADC(31) = sMDH.ulPMUTimeStamp;
           aADC(32) = sMDH.ulTimeSinceLastRF;
           aADC(35:35+MDH_NUMBEROFICEPROGRAMPARA-1) = sMDH.aushIceProgramPara;
           aADC(45:48) = sMDH.sSD.aflQuaternion;
           aADC(55:55+MDH_RESERVEDHDRPARA-1) = sMDH.aushReservedPara;
           aADC(60) = sMDH.sSD.sVector.flSag;
           aADC(61) = sMDH.sSD.sVector.flCor;
           aADC(62) = sMDH.sSD.sVector.flTra;
           aADC(64) = sMDH.ushCoilSelect;
           % ... ...
           aADC(IDX_DIM+1:end) = fread(fid, sMDH.ushSamplesInScan*2, 'float');
           
           
           
           
           if size(aADC,2) == size(data,2) % readouts of the same length
               data(lcounter,:) = aADC;
               
               
           %% CYL commented    
%            else % readouts of a different length
%     %           if size(aADC,2) == IDX_DIM + 32  % I regard this as a STOP sign so discard last N acquisitions of 32 (52-20) floats added by VB
%     %               disp(['WARNING: Discarded scan lcounter = ' num2str(lcounter) ' size(aADC,2) = ' num2str(size(aADC,2)-IDX_DIM) ' that is different from the previously detected ' num2str(size(data,2)-IDX_DIM)]);
%     %  		       break;
%     %           else % adjust the center of this readout by cutting short or padding zeros at the beginning of the readout
%                    data(lcounter,1:IDX_DIM) = aADC(1:IDX_DIM); % loopcounters
%                    delta = size(data,2)-size(aADC,2);
%                    if delta >=0
%                       data(lcounter,delta+IDX_DIM+1:end) = aADC(IDX_DIM+1:end);
%                     else
%                       data(lcounter,IDX_DIM+1:end) = aADC(-delta+IDX_DIM+1:end);
%                    end
%     %           end
%            end
           else %% CYL added for RTfeedback data 
               if size(data,2) > size(aADC,2),
                   data(lcounter,1:size(aADC,2)) = aADC;
               else
                   %% enlarge the data size
                   %PJN, change this line to avoid concatenation error
                  % data = [data zeros(NADC, size(aADC,2)-size(data,2))];
                   data = [data zeros(size(data,1), size(aADC,2)-size(data,2))];
                   data(lcounter,:) = aADC;
               end
           end
                  
         end
         % since for some reasons, several readouts are added at the end let keep the max
         if lcounter == 1
           maxSampleInScan = sMDH.ushSamplesInScan;
         else
           maxSampleInScan = max(maxSampleInScan, sMDH.ushSamplesInScan);
         end
         
         
         curr_pos = ftell(fid);
    %    if curr_pos >= hdr_len + Nbytes, LastADCData = 1; end 
         lcounter = lcounter + 1;
         lcha = lcha + 1;
      end % lcha
      %catch ME
      %    disp(['WARNING : Incomplete Measurement Data!']);
      %    break; % for incomplete data
      %end

      next_fid_pos = next_fid_pos + sMDH.ulDMALength; % beggining of each ADC readout
      % disp([curr_pos next_fid_pos next_fid_pos-curr_pos]);
      % if LastADCData, break; end

%       display('precutoff');
%       sMDH.sCutOff.ushPre
%       display('ushSamplesInScan');
%       sMDH.ushSamplesInScan
      
    end % of while

    
    
    % if feof(fid), break; end
    % put the max maxSampleInScan in sMDH
    % sMDH.ushSamplesInScan = maxSampleInScan;

    lcounter = lcounter - 1;
    data = data(1:lcounter,:);
    loopcounters = data(1:lcounter,1:IDX_DIM); % all start from 1, unused -1
    dimensions = max(loopcounters);

    %PJN Edits
    if numel(dimensions == 1)
        sMDH = sMDH; % return the 1st sMDH
        sMDH.sLC.ushLine = 1  ;
        sMDH.sLC.ushAcquisition = 1  ;
        sMDH.sLC.ushSlice = 1  ;
        sMDH.sLC.ushPartition = 1  ;
        sMDH.sLC.ushEcho = 1  ;
        sMDH.sLC.ushPhase = 1  ;
        sMDH.sLC.ushRepetition = 1  ;
        sMDH.sLC.ushSet = 1  ;
        sMDH.sLC.ushSeg = 1  ;
        sMDH.sLC.ushIda = 1 ;
        sMDH.sLC.ushIdb = 1 ;
        sMDH.sLC.ushIdc = 1 ;
        sMDH.sLC.ushIdd = 1 ;
        sMDH.sLC.ushIde = 1 ;
    else
    % Let sMDH have the maximum values of all loop counters
        sMDH = sMDH; % return the 1st sMDH
        sMDH.sLC.ushLine = dimensions(1)  ;
        sMDH.sLC.ushAcquisition = dimensions(2)  ;
        sMDH.sLC.ushSlice = dimensions(3)  ;
        sMDH.sLC.ushPartition = dimensions(4)  ;
        sMDH.sLC.ushEcho = dimensions(5)  ;
        sMDH.sLC.ushPhase = dimensions(6)  ;
        sMDH.sLC.ushRepetition = dimensions(7)  ;
        sMDH.sLC.ushSet = dimensions(8)  ;
        sMDH.sLC.ushSeg = dimensions(9)  ;
        sMDH.sLC.ushIda = dimensions(10) ;
        sMDH.sLC.ushIdb = dimensions(11) ;
        sMDH.sLC.ushIdc = dimensions(12) ;
        sMDH.sLC.ushIdd = dimensions(13) ;
        sMDH.sLC.ushIde = dimensions(14) ;
    end
    sMDH.sLC.ulChannelId = KMAX; % dimensions(15) ;
    % ... ...
    % let it bear the number of channels used
    sMDH.ulChannelId = KMAX;
    sCHA = sCHA1; % ruturn the 1st sCHA

   
    rawdata = complex(data(1:lcounter,IDX_DIM+1:2:end),data(1:lcounter,IDX_DIM+2:2:end));

    %%%%%%%%%%%%%%%%%%%%%%%% show you some information %%%%%%%%%%%%%%%%%%%%%

    lc_names = {'001:sMDH.sLC.ushLine'
           '002:sMDH.sLC.ushAcquisition'
           '003:sMDH.sLC.ushSlice'
           '004:sMDH.sLC.ushPartition'
           '005:sMDH.sLC.ushEcho'
           '006:sMDH.sLC.ushPhase'
           '007:sMDH.sLC.ushRepetition'
           '008:sMDH.sLC.ushSet'
           '009:sMDH.sLC.ushSeg'
           '010:sMDH.sLC.ushIda'
           '011:sMDH.sLC.ushIdb'
           '012:sMDH.sLC.ushIdc'
           '013:sMDH.sLC.ushIdd'
           '014:sMDH.sLC.ushIde'
           '015:sMDH.sLC.ulChannelId'
           '016:sMDH.ushSamplesInScan'
           '017:sMDH.ushKSpaceCentreColumn'
           '018:sMDH.fReadOutOffcentre'
           '019:sMDH.aulEvalInfoMask(1)'
           '020:sMDH.aulEvalInfoMask(2)' %20
           '021:' %21
           '022:' %22
           '023:sMDH.ushKSpaceCentreLineNo' %23
           '024:sMDH.sCutOff.ushPre' %24
           '025:sMDH.sCutOff.ushPost' %25
           '026:sMDH.ushKSpaceCentrePartitionNo' %26
           '027:sMDH.ulDMALength' %27
           '028:sMDH.lMeasUID' %28
           '029:sMDH.ulScanCounter' %29
           '030:sMDH.ulTimeStamp' %30
           '031:sMDH.ulPMUTimeStamp' %31
           '032:sMDH.ulTimeSinceLastRF' %32
           '033:' %33
           '034:' %34
           '035:' %35
           '036:' %36
           '037:' %37
           '038:' %38
           '039:' %39
           '040:' %40
           '041:' %41
           '042:' %42
           '043:' %43
           '044:' %44
           '045:sMDH.sSD.aflQuaternion(1)' %45
           '046:sMDH.sSD.aflQuaternion(2)' %46
           '047:sMDH.sSD.aflQuaternion(3)' %47
           '048:sMDH.sSD.aflQuaternion(4)' %48
           '049:' %49
           '050:' %50
           '051:' %51
           '052:' %52
           '053:' %53
           '054:' %54
           '055:sMDH.aushReservedPara(1)' %55
           '056:sMDH.aushReservedPara(2)' %56
           '057:sMDH.aushReservedPara(3)' %57
           '058:sMDH.aushReservedPara(4)' %58
           '059:' %59
           '060:sMDH.sSD.sVector.flSag' %60
           '061:sMDH.sSD.sVector.flCor' %61
           '062:sMDH.sSD.sVector.flTra' %62
           '063:' %63
           '064:sMDH.ushCoilSelect'}; %64

    disp('-------------------------- SUMMARY OF RAW DATA -----------------------');
    disp('------ 1st Scan Header Structure : ------');
    disp(sMDH); 
    disp('------ 1st Channel Header Structure : ------');
    disp(sCHA); 
    disp('------ Max Values of Loop Counters : ------');
    loopcounters_max = sMDH.sLC; disp(loopcounters_max);
    disp(['------ Max line index : ' num2str(sMDH.sLC.ushLine) ', Total number of readouts :', num2str(size(loopcounters, 1)) ' ------']);
    for jj = 1:15
      eval(['lcv = ' lc_names{jj,:}(5:end) ';'])
      if lcv > 1
          if jj ==1 %% since there are often many lines, give only the several first and last lines FYI
              disp(['repeats for each line -- loopcounters(' num2str(jj) ') or (' lc_names{jj,:} ') : --']);
              for kk = 1:10
                  disp(['  No. ' num2str(kk) ' = ' num2str(size(find(loopcounters(:, jj)==kk), 1)) ' repeats']);
              end 
              disp('    ... ...');
              for kk = sMDH.sLC.ushLine-10:sMDH.sLC.ushLine
                  disp(['  No. ' num2str(kk) ' = ' num2str(size(find(loopcounters(:, jj)==kk), 1)) ' repeats']);
              end
          else
              disp(['lines for each -- loopcounters(' num2str(jj) ') or (' lc_names{jj,:} ') : --']);
              for kk = 1:lcv
                  disp(['  No. ' num2str(kk) ' = ' num2str(size(find(loopcounters(:, jj)==kk), 1)) ' lines']);
              end
          end
      end
    end
    disp(['All the information masks used in the measurement data : --']);
    [B I] = unique(loopcounters(:, 19));
    for i = 1:size(B, 1);
        idx = find(loopcounters(I(i),20) == loopcounters(:, 20) & B(i) == loopcounters(:, 19));
        occur = size(idx, 1);
        disp(['Mask Values (' num2str([B(i) loopcounters(I(i),20)]) ')    =    '  num2str(occur) ' repeats :']);
        if exist('evalInfoMask.m','file')
          [M T] = evalInfoMask([B(i) loopcounters(I(i),20)]);
          disp(T(1:end,:));
        else
          disp(['M file (evalInfoMask.m) is not found to provide this information']);
        end
    end

    sMDH.sCutOff.ushPre
    aRAID.rawdata = rawdata;
    aRAID.loopcounters = loopcounters;
    aRAID.sMDH = sMDH;
    aRAID.MrProt = MrProt;
    aRAID.Yaps = Yaps;
    aRAID.lc_names = lc_names;
    aRAID.sCHA = sCHA;
    aRAID.RaidInfo = RaidInfo(iRAID);
    MultiRAID(iRAID) = aRAID;
    
end % iRAID


fclose(fid);


disp(['Information is retrieved using Matlab function : ']);
disp(['MultiRAID = ReadSiemensMeasVD13_idea(measfile, dispopt);']);
disp(['For details in MATLAB type help ReadSiemensMeasVD13_idea']);
disp('------------------- END OF DATA SUMMARY -------------------------');
disp(['Time used to read measurement data is : ' num2str(toc) ' seconds']);
disp(['The above information saved in diary file : ' diary_file]);

diary off;

if exist(diary_file,'file')
    delete(diary_file);
    delete(measfile_header);
end

%% try to show some way to rearrange the raw data and reconstrut the images
%% if not correct, following the examples below do it yourself AND ... ...
%% YOU CAN DO IT!!


if strcmp(dispopt, 'on')

    for iRAID = 1:size(MultiRAID, 2) 
        
        clear araw rawdata raw
        
        rawdata         = MultiRAID(iRAID).rawdata;
        loopcounters    = MultiRAID(iRAID).loopcounters;
        sMDH            = MultiRAID(iRAID).sMDH;
        MrProt          = MultiRAID(iRAID).MrProt;
        Yaps            = MultiRAID(iRAID).Yaps;
        lc_names        = MultiRAID(iRAID).lc_names;
        sCHA            = MultiRAID(iRAID).sCHA;
        RaidInfo        = MultiRAID(iRAID).RaidInfo;
        
        protname = deblank(RaidInfo.protname');
        
        what_to_see = 0;

        lin_counter = 0;
        for lc = 1:sMDH.sLC.ushLine

          switch what_to_see
          case 1
            lin_idx = find(loopcounters(:,1)   == lc & loopcounters(:,2)   == 1 & loopcounters(:,3)   == 1 & loopcounters(:,5)   == 1  &  loopcounters(:,15)  == 1); 
          case 2
            lin_idx = find(loopcounters(:,1)  == lc & loopcounters(:, 5)  == 1); % to see Echo = 1
          case 3
            lin_idx = find(loopcounters(:,1)  == lc & loopcounters(:,15)  ==  8); % to see Ch = 8
          otherwise
            lin_idx = find(loopcounters(:,1)  == lc & loopcounters(:,3)  == 1);
          end
          if ~isempty(lin_idx)
            lin_counter = lin_counter + 1;
            araw = rawdata(lin_idx,:);
            if size(araw, 1) == 1
              raw(lin_counter, :) = araw;
            else
              raw(lin_counter, :) = mean(araw, 1);
            end
          end
        end % lc

        raw = raw';
        raw_size = size(raw);
        aimg = fftshift(fft2(fftshift(raw)));
        mag =   abs(aimg); % magnitude image
        pha = angle(aimg); % phase image

        atitle = ['Slice:' num2str(1)];
        figure('Name', ['(' protname ') MAG,PHA,K.REAL,AND K.IM (' atitle ')']); 

        subplot(2, 2, 1), imagesc(mag'); axis image; colormap(gray);
        subplot(2, 2, 2), imagesc(pha'); axis image;
        subplot(2, 2, 3), imagesc(real(raw)'); axis image; 
        subplot(2, 2, 4), imagesc(imag(raw)'); axis image; 

    end %iRAID

end % dispopt

