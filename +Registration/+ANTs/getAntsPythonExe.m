% RH: new file — locates the ANTsPy-capable Python interpreter for the macOS/Linux
% fallback used by AntsRegistration.m and N4_bias_correction.m
function pythonExe = getAntsPythonExe(registration_path)
%GETANTSPYTHONEXE Locate a Python interpreter with antspyx installed.
%   Used on macOS/Linux, where no native antsRegistration CLI binary is
%   bundled, to run ANTs registration through ANTsPy instead.
%
%   Reads the interpreter path from <registration_path>/ants_python_path.txt
%   (one line, e.g. the python binary inside a conda env that has
%   `pip install antspyx`). Falls back to "python3" on PATH if the file
%   is missing.

    configFile = fullfile(registration_path, 'ants_python_path.txt');
    if exist(configFile, 'file')
        fid = fopen(configFile, 'r');
        pythonExe = strtrim(fgetl(fid));
        fclose(fid);
    else
        pythonExe = 'python3';
    end
end
