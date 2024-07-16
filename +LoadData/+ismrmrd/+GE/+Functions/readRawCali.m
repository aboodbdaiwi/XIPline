function [cali_struct] = readRawCali(raw_path)
%Read in twix file or mrd file.
%   readRawCali(X) reads in the twix file or mrd file of the 129Xe
%   calibration scan and returns a structure containing all the variables
%   needed for processing.
% The data structure file contains:
%
% cali_struct.seq_name = sequence name (optional)
% cali_struct.weight = patient weight in pounds %
% cali_struct.te = TE in u-seconds;
% cali_struct.tr = TR in u-seconds;
% cali_struct.dwell_time = dwell time in seconds;
% cali_struct.freq = gas excitation frequency in Hz;
% cali_struct.xeFreqMHz = gas excitation frequency in MHz;
% cali_struct.data = the FID data;
% cali_struct.scan_date = scan date in format (YYYY-MM-DD) optional;
% cali_struct.vref = reference voltage (V) optional;
% cali_struct.rf_excitation_ppm = rf excitation in ppm;

cali_struct = {};
[~, ~, file_extension] = fileparts(raw_path);

% Read in twix or P file and define associated variables
switch file_extension
    case '.dat'
        % Twix file from Siemens
        twix_obj = mapVBVD(raw_path);
        twix_obj.data = squeeze(double(twix_obj.image()));
        % If twix obj contains field sWipMemBlock, change to sWiPMemBlock (capital P)
        if isfield(twix_obj.hdr.MeasYaps, 'sWipMemBlock')
            temp = RenameField(twix_obj.hdr.MeasYaps, 'sWipMemBlock', 'sWiPMemBlock');
            twix_obj.hdr.MeasYaps = temp;
        end

        if isfield(twix_obj.hdr.Phoenix, 'sWipMemBlock')
            temp = RenameField(twix_obj.hdr.Phoenix, 'sWipMemBlock', 'sWiPMemBlock');
            twix_obj.hdr.Phoenix = temp;
        end
        % Extract variables and store in cali_struct
        npts = size(twix_obj.data, 1); % Number of samples per FID
        dwell_time = twix_obj.hdr.MeasYaps.sRXSPEC.alDwellTime{1, 1} * 10^-9; % Receiver bandwidth (kHz); works for both calibration types
        tr = twix_obj.hdr.Config.TR(1) * 1E-6; % Time between each sample
        fids = twix_obj.data;
        % Get the excitation frequency
        if isfield(twix_obj.hdr.Config, 'Frequency')
            % UVA Siemens File
            xeFreqMHz = twix_obj.hdr.Config.Frequency * 1e-6; %34.093484
        elseif isfield(twix_obj.hdr.Meas, 'lFrequency')
            % Duke Siemens File
            xeFreqMHz = twix_obj.hdr.Meas.lFrequency * 1e-6; %34.091516
        end
        % Get the scan date
        scanDate = twix_obj.hdr.Phoenix.tReferenceImage0;
        scanDate = strsplit(scanDate, '.');
        scanDate = scanDate{end};
        scanDateStr = [scanDate(1:4), '-', scanDate(5:6), '-', scanDate(7:8)];
        if isfield(twix_obj.hdr.Phoenix, 'sWiPMemBlock')
            % Duke twix file
            if isfield(twix_obj.hdr.Phoenix.sWiPMemBlock, 'adFree') && length(twix_obj.hdr.Phoenix.sWiPMemBlock.adFree)>3
                VRef = twix_obj.hdr.Phoenix.sWiPMemBlock.adFree{4};
                % seems to be in all sequences
                %rf_amp1 = twix_obj.hdr.Phoenix.sTXSPEC.aRFPULSE{1}.flAmplitude;
                % in cali this is the 3rd flip angle (calibration, usually)
                %rf_amp2 = twix_obj.hdr.Phoenix.sTXSPEC.aRFPULSE{3}.flAmplitude;
                %rf_amp3=1; %dummy override for when using hard pulse
            elseif isfield(twix_obj.hdr.Phoenix.sWiPMemBlock, 'alFree')
                % Duke twix file using UVA sequence
                VRef = twix_obj.hdr.Phoenix.sWiPMemBlock.alFree{1};
                % dissolved phase
                %rf_amp1 = twix_obj.hdr.Phoenix.sTXSPEC.aRFPULSE{2}.flAmplitude;
                % calibration
                %rf_amp2 = twix_obj.hdr.Phoenix.sTXSPEC.aRFPULSE{3}.flAmplitude;
            else
                disp('WARNING: twix file type not supported, cannot determine reference voltage')
            end
        else
            disp('WARNING: twix file type not supported, cannot determine reference voltage')
        end
        % Read RF excitation frequency
        mag_fstrength = twix_obj.hdr.Dicom.flMagneticFieldStrength; % Magnetic Field Strength
        excitation = twix_obj.hdr.Phoenix.sWiPMemBlock.alFree{1, 5}; % Read excitation from twix header,Cali version
        
        % RF excitation will be in ppm, likeley either 218ppm or 208 ppm at Duke
        gyro_ratio = 11.777; % gyromagnetic ratio of 129Xe in MHz/Tesla
        rf_excitation_ppm = round(excitation/(gyro_ratio * mag_fstrength));

        % Save to output struct
        cali_struct.seq_name = twix_obj.hdr.Config.SequenceDescription;
        cali_struct.weight = twix_obj.hdr.Dicom.flUsedPatientWeight; %
        cali_struct.te = twix_obj.hdr.Phoenix.alTE{1};
        cali_struct.tr = twix_obj.hdr.Config.TR(1);
        cali_struct.dwell_time = twix_obj.hdr.MeasYaps.sRXSPEC.alDwellTime{1, 1}* 1E-9;
        cali_struct.freq = xeFreqMHz * 1e6;
        cali_struct.xeFreqMHz = xeFreqMHz;
        cali_struct.data = twix_obj.data;
        cali_struct.scan_date = scanDateStr;
        cali_struct.vref = VRef;
        cali_struct.rf_excitation_ppm = rf_excitation_ppm;

    case '.mrd'
        % Removed '.h5' from the line above
        % mrd file

        % read in mrd file dataset and ismrmrdHeader
        dataset = LoadData.ismrmrd.Dataset(raw_path, 'dataset');
        ismrmrd_header = LoadData.ismrmrd.xml.deserialize(dataset.readxml);

        % convert user parameter fields to maps for easy query
        general_user_params_long = containers.Map();
        data_struct = ismrmrd_header.userParameters.userParameterLong;
        for i = 1:numel(data_struct)
            general_user_params_long(data_struct(i).name) = data_struct(i).value;
        end

        % read in variables
        cali_struct.scan_date = ismrmrd_header.studyInformation.studyDate; % in YYYY-MM-DD
        vendor = ismrmrd_header.acquisitionSystemInformation.systemVendor;
        cali_struct.te = ismrmrd_header.sequenceParameters.TE * 1e3; % in us
        cali_struct.tr = ismrmrd_header.sequenceParameters.TR(2) * 1e3; % dissolved TR in us
        cali_struct.dwell_time = double(dataset.readAcquisition().head.sample_time_us(1)) * 1e-6; % in s
        cali_struct.freq = general_user_params_long("xe_center_frequency"); % in Hz
        cali_struct.xeFreqMHz = cali_struct.freq * 1e-6; % gas excitation frequency in MHz
        field_strength = ismrmrd_header.acquisitionSystemInformation.systemFieldStrength_T; % in T
        freq_dis_excitation_hz = general_user_params_long("xe_dissolved_offset_frequency"); % in Hz

        % calculate rf excitation in ppm
        gyro_ratio = 11.777; % gyromagnetic ratio of 129Xe in MHz/Tesla
        cali_struct.rf_excitation_ppm = round(freq_dis_excitation_hz/(gyro_ratio * field_strength));

        % assign nan to variables not in mrd file
        cali_struct.seq_name = nan;
        cali_struct.weight = nan;
        cali_struct.vref = nan;

        % read k-space data
        npts = size(dataset.readAcquisition(1).data{1},1);
        nfids = dataset.getNumberOfAcquisitions;
        fids_cell = dataset.readAcquisition().data;
        fids = zeros(npts, nfids);
        for i=1:nfids
            fids(:,i) = transpose(double(fids_cell{i}(:,1)));
        end
        
        % if data from GE scanner, take complex conjugate
        if strcmpi(vendor,'ge')
            fids = conj(fids);
        end
        
        % add fid data to struct
        cali_struct.data = fids;
    case '.h5'
        % Removed '.h5' from the line above
        % mrd file

        % read in mrd file dataset and ismrmrdHeader
        dataset = LoadData.ismrmrd.Dataset(raw_path, 'dataset');
        ismrmrd_header = LoadData.ismrmrd.xml.deserialize(dataset.readxml);

        % convert user parameter fields to maps for easy query
        general_user_params_long = containers.Map();
        data_struct = ismrmrd_header.userParameters.userParameterLong;
        for i = 1:numel(data_struct)
            general_user_params_long(data_struct(i).name) = data_struct(i).value;
        end

        % read in variables
        cali_struct.scan_date = ismrmrd_header.studyInformation.studyDate; % in YYYY-MM-DD
        vendor = ismrmrd_header.acquisitionSystemInformation.systemVendor;
        cali_struct.te = ismrmrd_header.sequenceParameters.TE * 1e3; % in us
        cali_struct.tr = ismrmrd_header.sequenceParameters.TR(2) * 1e3; % dissolved TR in us
        cali_struct.dwell_time = double(dataset.readAcquisition().head.sample_time_us(1)) * 1e-6; % in s
        cali_struct.freq = general_user_params_long("xe_center_frequency"); % in Hz
        cali_struct.xeFreqMHz = cali_struct.freq * 1e-6; % gas excitation frequency in MHz
        field_strength = ismrmrd_header.acquisitionSystemInformation.systemFieldStrength_T; % in T
        freq_dis_excitation_hz = general_user_params_long("xe_dissolved_offset_frequency"); % in Hz

        % calculate rf excitation in ppm
        gyro_ratio = 11.777; % gyromagnetic ratio of 129Xe in MHz/Tesla
        cali_struct.rf_excitation_ppm = round(freq_dis_excitation_hz/(gyro_ratio * field_strength));

        % assign nan to variables not in mrd file
        cali_struct.seq_name = nan;
        cali_struct.weight = nan;
        cali_struct.vref = nan;

        % read k-space data
        npts = size(dataset.readAcquisition(1).data{1},1);
        nfids = dataset.getNumberOfAcquisitions;
        fids_cell = dataset.readAcquisition().data;
        fids = zeros(npts, nfids);
        for i=1:nfids
            fids(:,i) = transpose(double(fids_cell{i}(:,1)));
        end
        
        % if data from GE scanner, take complex conjugate
        if strcmpi(vendor,'ge')
            fids = conj(fids);
        end
        
        % add fid data to struct
        cali_struct.data = fids;
    otherwise
        error('Unknown Raw File Type');
end
end
