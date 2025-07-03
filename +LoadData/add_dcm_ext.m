% Define the folder path
folderPath = 'C:\Your\Folder\Path';  % <-- Replace with your folder path
% Get list of all files in the folder
fileList = dir(folderPath);
% Loop through each file
for i = 1:length(fileList)
    % Skip folders
    if fileList(i).isdir
        continue;
    end
    % Get full file name and extension
    oldName = fileList(i).name;
    [~, name, ext] = fileparts(oldName);

    % Skip if already has .dcm extension
    if strcmpi(ext, '.dcm')
        continue;
    end
    % Construct full file paths
    oldFullPath = fullfile(folderPath, oldName);
    newFullPath = fullfile(folderPath, [name ext '.dcm']);
    % Rename file
    movefile(oldFullPath, newFullPath);
end
disp('Renaming completed.');
