function series_description = rename_series(dataPath,common_string)
folder_content=dir(fullfile(dataPath,common_string));
for n_f=1:numel(folder_content)
    series_path = fullfile(dataPath, folder_content(n_f).name);
    dicom_file_list = dir(series_path);
    if ~isempty(dicom_file_list)
        dicom_info = dicominfo(fullfile(series_path,dicom_file_list(end).name));
        fprintf('patient name = %s; scan date=%s\n',dicom_info.PatientName.FamilyName,dicom_info.StudyDate)
        series_description = strrep(dicom_info.SeriesDescription,':','');
        series_description = strrep(series_description,'*','');
        new_series_folder =fullfile(dataPath, series_description);
        if ~isdir(new_series_folder)
            mkdir(new_series_folder)
        end
       movefile(fullfile(series_path,'I*'),new_series_folder,'f');
       rmdir(series_path);
        
    end
end