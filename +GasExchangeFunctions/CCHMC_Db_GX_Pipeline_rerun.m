function CCHMC_Db_GX_Pipeline_rerun(MainInput, analysispath)
    clc;
    UpdatedImageQuality = MainInput.UpdatedImageQuality;
    UpdatedNote = MainInput.UpdatedNote;
    disp('loading workspace, please wait.....')
    load(fullfile(MainInput.gx_analysis_folder, 'workspace.mat'));
    
    MainInput.UpdatedImageQuality = UpdatedImageQuality;
    MainInput.UpdatedNote = UpdatedNote;

    % mask_file_name = dir(fullfile(MainInput.gx_analysis_folder, 'Mask.nii.gz'));
    % mask_file_name = fullfile(mask_file_name.folder, mask_file_name.name);
    mask_file_name = fullfile(MainInput.gx_analysis_folder, 'Mask.nii.gz');
    [~,~,mask_ext] = fileparts(mask_file_name);
    if strcmp(mask_ext, '.gz') || strcmp(mask_ext, '.nii')
        try
            Mask = LoadData.load_nii(mask_file_name);
        catch
            Mask = LoadData.load_untouch_nii(mask_file_name);
        end
        A = double(Mask.img);
    elseif strcmp(mask_ext, '.dcm')
        A = double(squeeze(dicomread(mask_file_name)));
    end
    GasExchange.LungMask = A;
    GasExchange.Image = double(abs(GasExchange.VentImage));
    [~, GasExchange] = VentilationFunctions.correct_mask_orientation(GasExchange);
    % figure; imslice(GasExchange.LungMask)
    Proton.ProtonMaskRegistred = GasExchange.LungMask;
    Proton.LungMask = GasExchange.LungMask;
    GasExchange.Mask = GasExchange.LungMask;
    
    delete_if_exist = @(pattern) cellfun(@(f) delete(fullfile(MainInput.gx_analysis_folder, f)), ...
        {dir(fullfile(MainInput.gx_analysis_folder, pattern)).name}, 'UniformOutput', false);
    
    % ==== Delete any older matching files ====
    cd(MainInput.gx_analysis_folder)
    delete_if_exist('GasExchange_Report*.pptx*');
    
    if ~isfield(MainInput, 'analysispath') 
         MainInput.analysispath = analysispath;
    end

    [GasExchange] = GasExchangeFunctions.GasExchange_Analysis(GasExchange,Proton,MainInput);
    
    clearvars
end

