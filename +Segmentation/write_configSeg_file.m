function write_configSeg_file(MainInput)
% write_config_file - Writes config_Xe.txt and/or config_H.txt
% to C:\XIPline\offline_recon, removing old versions if present.

    % RH: was hardcoded to 'C:\XIPline', which doesn't exist on macOS/Linux;
    % branch like the rest of the segmentation pipeline does.
    if ispc
        outputFolder = 'C:\XIPline';
    else
        outputFolder = fullfile(getenv('HOME'), 'XIPline');
    end

    % if ~exist(outputFolder, 'dir')
    %     mkdir(outputFolder);
    % end

    % Define config file paths
    configSegPath = fullfile(outputFolder, 'config_segmentation.txt');

    % -------- Delete old config files if they exist --------
    if exist(configSegPath, 'file')
        delete(configSegPath);
        fprintf('Deleted old configSegPath.txt\n');
    end


    % -------- Write Segmentation Config --------
    if isfield(MainInput, 'SegmentType') 
        fid = fopen(configSegPath, 'w');
        if fid == -1
            error('Unable to open file for writing: %s', configSegPath);
        end
        fprintf(fid, 'SegmentType = %s\n', MainInput.SegmentType);
        fclose(fid);
        fprintf('Wrote new configSegPath.txt\n');
    end

end
