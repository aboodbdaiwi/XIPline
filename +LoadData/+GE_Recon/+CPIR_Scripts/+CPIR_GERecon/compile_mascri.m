%% options
name_of_executable = 'matlab_scripter';
matlab_file_name_to_compile = 'matlab_scripter_cchmc.m';
options = '-Cmv';
% C= do not embed deployable archive
% m= generate standalone application
% v= dispaly verbose output

%% compile
v=strsplit(version,'.');
folder_v = [v{1} v{2}];
output_path = ['mascri',folder_v,'_cchmc'];

if exist(output_path,"dir") % backup old versions
    status = movefile(output_path,[output_path,'_old_',datestr(now,'yyyymmdd')]);
end

mcc('-o', name_of_executable,'-A','glnxa64', '-d',output_path, options, matlab_file_name_to_compile)