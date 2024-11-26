function defect_mask_kmeans=defect_analysis_script(proton,he,lungsmask,study_id,output_path)
% required inputs: 1.proton --registered to 3He preferablly
%                  2.he
%                  3.lungsmask --
%                        1) If binary, left and right lungs will be
%                           separated when calling identify_lobes. Please note that
%                           it is currently not able to seprate right and left lungs
%                           if connected.
%                        2) this input (lungsmask) could be nonbinary such that each
%                           value representing a lung lobe. The values and matching
%                           lobe labels are:
%                            2 for RUL
%                            3 for RML
%                            4 for RLL
%                            5 for LUL
%                            6 for LLL
%                            and the subsequent output will have defects by lobes;
%
%                        3) If only right and left lungs were separated, the value 2
%                           represents right lung and 5 represents left lung.
%
% optional inputs:
%                  -study_id: a string used as output filename base;
%                  -output_path: directory to save the output;
%
% NOTE: It is important to have a lung mask that matches with proton and
% 3He for estimating vessel mask and segmented defects because every pixel
% inside the lung mask will be considered for defect.
%
%  W. Zha @2015
%% Check inputs
if nargin<4
    path_struct.subject='Unknown study';
    path_struct.output_path=pwd;
elseif nargin<5
    path_struct.subject=study_id;
    path_struct.output_path=pwd;
    
else
    path_struct.subject=study_id;
    path_struct.output_path=output_path;%'/data/users/mummy/sandbox/defect_testing/SARP80-811-020/';
    
end
%% bias correction on 3He, which requires the installation of ANTs at https://github.com/stnava/ANTs.
try
    corrected_he=bias_correction_mat_data(he,path_struct);
    temp_handles.he=corrected_he;
catch
    temp_handles.he=he;
    warning('defect_analysis_script:bias_correction_mat_data','not able to perform intensity correction on 3He.\n');
end

%% set up inputs for defect analysis
temp_handles.proton = proton;

values = unique(lungsmask);
values(values==0)=[];
if isempty(values)
    error('defect_analysis_script:lungsmask','invalid lungsmask.\n');
elseif length(values)>1
    temp_handles.lungsmask=lungsmask;
    temp_handles.lungsmask(temp_handles.lungsmask>0)=1;
    temp_handles.lobemask=lungsmask;
elseif values==1
    temp_handles.lungsmask=lungsmask;
    temp_handles.lobemask=lungsmask;
    
end



temp_handles = check_mask_size(temp_handles);
temp_handles = identify_lobes(temp_handles);
% additional parameters
temp_handles.hist_mu=385; % intensity value corresponds to the 75% landmark on the histogram of proton image ranged in [0,512]
temp_handles.frangi_thresh=0.3;% empirical choice for binarize the vessel mask
temp_handles.vessel_removal = false;%true;
temp_handles.sigmas=2.5; % empirical choice for vesselness filter
temp_handles.plotting_range=get_plotting_range(temp_handles.lungsmask);

%% defect analysis
defect_mask_kmeans = detect_defects_by_kmeans_sarp(temp_handles);
%% visualize defects on 3He

title_string = sprintf('%s defect by Kmeans',strrep(path_struct.subject,'_',' '));
slice_range=80:100;
temp_handles.slice_range=slice_range;
figure,plot_defect_mask(temp_handles.he(:,:,slice_range),defect_mask_kmeans(:,:,slice_range),temp_handles,temp_handles.plotting_range,title_string);
figure,plot_defect_mask(temp_handles.he(:,:,slice_range),temp_handles.lungsmask(:,:,slice_range),temp_handles,temp_handles.plotting_range,title_string);
figure,plot_defect_mask(temp_handles.proton(:,:,slice_range),temp_handles.lungsmask(:,:,slice_range),temp_handles,temp_handles.plotting_range,title_string);