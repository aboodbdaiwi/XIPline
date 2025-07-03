function [Ventilation, MainInput] = philips_XeVent_2DSpiral_recon(MainInput, Ventilation)
% run_xe_reconstruction - Runs Xe image reconstruction using external exe
% and loads the resulting NIfTI file after checking for existence.
%
% Inputs:
%   MainInput - Struct with required fields:
%       .XeDatalocation: path to Xe data
%       .ReconImageMode: 'xenon'
%       .freqoffset
%       .XeFileName: output image filename
%       .XeFullPath: full path to img_ventilation.nii.gz
%   LoadData - Struct with fields:
%       .write_config_file(MainInput): function handle to write config
%       .load_nii(path): load NIfTI
%       .load_untouch_nii(path): fallback NIfTI loader
%
% Output:
%   Ventilation - Struct containing loaded image and metadata

    % Step 1: Write config file
    LoadData.write_config_file(MainInput);

    % Step 2: Run external recon executable
    exePath = 'C:\XIPline\offline_recon\Philips_2DXeSpiralRecon.exe';
    [status, cmdout] = system(['"', exePath, '"']);

    if status ~= 0
        error('Failed to run Philips_2DHSpiralRecon.exe:\n%s', cmdout);
    end

    % Step 3: Wait for img_ventilation.nii.gz
    outputFile = fullfile(MainInput.XeDataLocation, 'img_ventilation.nii.gz');
    waitTimes = [1, 5, 10];  % seconds
    fileFound = false;

    for t = waitTimes
        pause(t);
        if exist(outputFile, 'file')
            fprintf('File found after waiting %d seconds: %s\n', t, outputFile);
            fileFound = true;
            break;
        else
            fprintf('File not found after waiting %d seconds. Retrying...\n', t);
        end
    end

    % Step 4: Load and process image
    if fileFound
        try
            A1 = LoadData.load_nii(outputFile);
        catch
            A1 = LoadData.load_untouch_nii(outputFile);
        end
        A = double(squeeze(A1.img));
        A = permute(A,[1, 3, 2]);
        A = imrotate(A, 90);
        A = flip(A, 2);
        A = flip(A, 3);
        %imslice(A)
        % Assign to output struct
        Ventilation.Image    = A;
        Ventilation.filename = MainInput.XeFileName;
        Ventilation.folder   = MainInput.XeDataLocation;
    else
        error('img_ventilation.nii.gz not found after 150 seconds. Aborting.');
    end

end
