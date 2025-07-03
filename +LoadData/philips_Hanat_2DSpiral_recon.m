function [Proton, MainInput] = philips_Hanat_2DSpiral_recon(MainInput, Proton)
% philips_Hanat_2DSpiral_recon - Runs proton image reconstruction
% using external spiral recon executable and loads the resulting NIfTI.
%
% Inputs:
%   MainInput - Struct with required fields:
%       .HDataLocation: path to proton data
%       .ReconImageMode: 'proton'
%       .NoProtonImage: 0 or 'no' if proton image exists
%
% Output:
%   Proton - Struct with loaded image and metadata

    % Step 1: Write config file
    LoadData.write_config_file(MainInput);

    % Step 2: Run external recon executable
    exePath = 'C:\XIPline\offline_recon\Philips_2DHSpiralRecon.exe';
    [status, cmdout] = system(['"', exePath, '"']);

    if status ~= 0
        error('Failed to run Philips_2DHSpiralRecon.exe:\n%s', cmdout);
    end

    % Step 3: Wait for img_proton.nii.gz
    outputFile = fullfile(MainInput.HDataLocation, 'img_proton.nii.gz');

    if exist(outputFile, 'file')
        fprintf('Proton file already exists: %s\n', outputFile);
    else
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
        if ~fileFound
            error('img_proton.nii.gz not found after waiting. Aborting.');
        end
    end

    % Step 4: Load and process image
    try
        A1 = LoadData.load_nii(outputFile);
    catch
        A1 = LoadData.load_untouch_nii(outputFile);
    end
    A = double(squeeze(A1.img));
    A = permute(A, [1, 3, 2]);
    A = imrotate(A, 90);
    A = flip(A, 2);
    A = flip(A, 3);

    % Assign to output struct
    Proton.Image    = A;
    Proton.filename = 'img_proton.nii.gz';
    Proton.folder   = MainInput.HDataLocation;
end
