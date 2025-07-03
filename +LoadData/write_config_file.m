function write_config_file(MainInput)
% write_config_file - Writes config_Xe.txt and/or config_H.txt
% to C:\XIPline\offline_recon, removing old versions if present.

    outputFolder = 'C:\XIPline\offline_recon';

    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
    end

    % Define config file paths
    configXePath = fullfile(outputFolder, 'config_Xe.txt');
    configHPath  = fullfile(outputFolder, 'config_H.txt');

    % -------- Delete old config files if they exist --------
    if exist(configXePath, 'file')
        delete(configXePath);
        fprintf('Deleted old config_Xe.txt\n');
    end
    if exist(configHPath, 'file')
        delete(configHPath);
        fprintf('Deleted old config_H.txt\n');
    end

    % -------- Write Xenon Config --------
    if isfield(MainInput, 'ReconImageMode') && strcmpi(MainInput.ReconImageMode, 'xenon')
        fid = fopen(configXePath, 'w');
        if fid == -1
            error('Unable to open file for writing: %s', configXePath);
        end
        fprintf(fid, 'datalocation = %s\n', MainInput.XeDataLocation);
        fprintf(fid, 'freq_offset = %d\n', MainInput.freqoffset);
        fclose(fid);
        fprintf('Wrote new config_Xe.txt\n');
    end

    % -------- Conditionally Write Proton Config --------
    if isfield(MainInput, 'NoProtonImage')
        hasProton = ...
            (isnumeric(MainInput.NoProtonImage) && MainInput.NoProtonImage == 0) || ...
            (ischar(MainInput.NoProtonImage) && strcmpi(MainInput.NoProtonImage, 'no'));

        if hasProton && isfield(MainInput, 'HDataLocation')
            fid = fopen(configHPath, 'w');
            if fid == -1
                error('Unable to open file for writing: %s', configHPath);
            end
            fprintf(fid, 'datalocation = %s\n', MainInput.HDataLocation);
            fclose(fid);
            fprintf('Wrote new config_H.txt\n');
        else
            fprintf('Proton config not written (either NoProtonImage is false or HDatalocation is missing).\n');
        end
    end
end
