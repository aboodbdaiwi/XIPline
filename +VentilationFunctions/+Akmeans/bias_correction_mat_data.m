function [corrected_data,bias_field, bias_corrected]=bias_correction_mat_data(raw_3Ddata,path_struct)
%
% ANTs function N4BiasFieldCorrection takes the raw dicom images as input
% and iteratively estimate the B1 inhomogeneity as a slow-variant, B-spline
% function. 
%
% corrected_data = raw_3Ddata./bias_field;
%
% Ref:Tustison et al. N4ITK: Improved N3 bias correction. IEEE Trans Med 
% Imaging. 2010
%
% W. Zha @2014
%%
bias_corrected=false;
[nRows,nCols,nSlices]=size(raw_3Ddata);
%% prepare .mha files for ANTs function call
mha_dest_folder = path_struct.output_path;
subject_string=path_struct.subject;
estimation_file_basename=sprintf('%s_for_fieldestm',subject_string);
write_dot_dat_dot_mha(estimation_file_basename,mha_dest_folder,double(raw_3Ddata));
output_filename=sprintf('%s_biascorrected',subject_string);

%% bias correction
cd(mha_dest_folder)
system(['N4BiasFieldCorrection -d ',num2str(3), ' -i ',...
    strcat(estimation_file_basename,'.mha'),' -s ',num2str(2),...
    ' -o ',strcat(output_filename,'.mha')]);

%% convert back to matrix format
system(['ImageFileConverter ',strcat(output_filename,'.mha')]);

corrected_data = read_in(strcat(output_filename,'.dat'),[nRows nCols nSlices],0,'l','float');
bias_field = corrected_data./raw_3Ddata;
if ~isempty(corrected_data)
    bias_corrected=true;
end