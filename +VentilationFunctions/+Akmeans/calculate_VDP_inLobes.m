studyid = '80-811-011';

data_folder =fullfile('/data/users/mummy/testdata/svdp_debacle/011_rework');
visit = 'mribl';



data_path.he = fullfile(data_folder,visit,'helium');
data_path.proton = fullfile(data_folder,visit,'proton');
data_path.ct = fullfile(data_folder,'ct');

vida_path = fullfile(data_folder,'sublobe');

addpath(genpath('/data/users/mummy/m-files/NIfTI_20140122'));
addpath(genpath('/data/users/mummy/m-files/tools'));
addpath(genpath('/data/users/mummy/m-files/defect_analysis_core'));

sublobe_lookup_init;
%% binary proton lung mask: *.mat
proton_lungmask_matfile = fullfile(data_folder,'lungsmask_011.mat'); % variable name: mask



%% HARD CODED VIDA OUTPUT
vida_output = 'ZUNU_vida-sublobes.img.gz';

ants_folder='/export/home/zha/ants_v210/antsbin/bin';
script_folder='/data/users/zha/antScripts';

output_folder = fullfile(data_folder, '/output');
sublobe_folder = fullfile(data_folder, '/sublobe');

if ~isdir(output_folder)
    mkdir(output_folder)
end
cd(output_folder);
%% READ spreadsheet for matched CT slice list
ct_proton_sheet = sprintf('%s.xlsx',studyid);
ct_slice_matfile = fullfile(data_folder,sprintf('ct_slice%s.mat',studyid));
if exist(fullfile(data_folder,ct_proton_sheet),'file')
[~,~,raw_sheet] = xlsread(fullfile(data_folder,ct_proton_sheet));
proton_sl = cell2mat(raw_sheet(2:end,1));

ct_slice = cell2mat(raw_sheet(2:end,2));
elseif exist(ct_slice_matfile,'file')
    load(ct_slice_matfile); % variable name: ct_slice
end

%% convert data into NIfTI
data_types = fieldnames(data_path);
for n_type = 1:numel(data_types) % DGM 4/5/16 added indexing
    data_type = data_types{n_type};
   % [rawdata.(data_type),image_headers.(data_type)] = load_dicom_series(data_path.(data_type),'IM');
    [rawdata.(data_type),acquisition_order,image_headers.(data_type),IPP] = load_dicom_series(data_path.(data_type),'IM');
   if strcmp(data_type,'ct')
        nii_data.(data_type) = make_nii(double(rawdata.(data_type)(:,:,ct_slice)),[],[],64);
        
    else
        nii_data.(data_type) = make_nii(double(rawdata.(data_type)),[],[],64);
    end
    cd(output_folder)
    nii_files.(data_type) =fullfile(output_folder, sprintf('sarp%s_raw_%s.nii.gz',studyid,data_type));
    save_nii(nii_data.(data_type),nii_files.(data_type));
    if ~strcmp(data_type,'ct')
        % bias correction
        nii_files.cor.(data_type) = strrep(nii_files.(data_type),'raw','biasC');
        system([ants_folder,'/N4BiasFieldCorrection -d 3 -v 1 -i ',...
            nii_files.(data_type),' -s ',num2str(2),' -o ',nii_files.cor.(data_type)]);
        nii_data.cor.(data_type) = load_nii(nii_files.cor.(data_type));
    else
        
    end
end


%% true rigid registration: proton (moving)->he(fixed)
reg_pair.moving = nii_files.proton;
reg_pair.fixed = nii_files.cor.he;
moving_data = 'proton';
    system([script_folder,'/rigid_true ',reg_pair.moving,' ', reg_pair.fixed,' ',moving_data]);
proton_reg_file = strrep(reg_pair.fixed,'he', 'he_fixed_proton_movingWarped');
reg_data = load_nii(proton_reg_file);
title_string = sprintf('%s %s registered',studyid,moving_data);
proton_gas_overlay(reg_data.img,nii_data.cor.he.img,title_string);
proton_gas_overlay(nii_data.proton.img,nii_data.cor.he.img,strrep(title_string,'registered','original'));

%% VIDA results
system(['gunzip ',fullfile(sublobe_folder,vida_output)]);
sublobar_mask = analyze75read(fullfile(sublobe_folder,strrep(vida_output,'.gz','')));
%sublobar_mask_BA = sublobar_mask(:,:,end:-1:1);% wz: Check if it is always needed

%ct_sublobar = sublobar_mask_BA(:, :, ct_slice); % DGM 4/5/16 adding line of code to get sublobar subset
ct_sublobar = sublobar_mask(:, :, ct_slice); % DGM 4/5/16 adding line of code to get sublobar subset

[nRows,nCols,nSlices] = size(ct_sublobar);
lungmask.ct = zeros(nRows,nCols,nSlices);
lungmask.ct(ct_sublobar>0)=1;

 figure,imshow3D(lungmask.ct,[])



lost_regions = setdiff(unique(sublobar_mask),unique(ct_sublobar));
fprintf(' lost sublobar region by selecting matched CT slices: %d\n',lost_regions);
figure,data_mask_overlay(nii_data.ct.img,lungmask.ct,'ct raw') % DGM Don't have this? 

%% Pre-registration processing
% Downsample CT
[resized_ct,lungmask.ct_resz] = resample_CT_to_proton(nii_data.ct.img, image_headers.ct,image_headers.proton,ct_sublobar);
resized_ct(resized_ct<0)=0;
CTmin = min(resized_ct(:));
CTmax = max(resized_ct(:));
nii_data.nlz.ct = make_nii((resized_ct-CTmin)./(CTmax-CTmin),[],[],64);
reg_pair2.moving = strrep(nii_files.ct,'raw','resized_nlz');
% save_nii(nii_data.nlz.ct,reg_pair2.moving);
figure,data_mask_overlay(nii_data.nlz.ct.img,lungmask.ct_resz,'ct resized')


% enhance and normalize proton
enhanced_proton_reg = hist_transformation(reg_data.img,385);
Protonmin = min(double(enhanced_proton_reg(:)));
Protonmax = max(double(enhanced_proton_reg(:)));
nii_data.nlz.proton = make_nii((enhanced_proton_reg-Protonmin)./(Protonmax-Protonmin),[],[],64);
reg_pair2.fixed = fullfile(output_folder,sprintf('sarp%s_reg_nlz_%s.nii.gz',studyid,moving_data));
% save_nii(nii_data.nlz.proton,reg_pair2.fixed);
%% deformable registration: ct(moving)->proton(fixed)
moving_data2 = 'ct';
%ct_reg = zeros(nRows,nCols,nSlices);
%sublobar_reg = zeros(nRows,nCols,nSlices);
%proton_ct_fused = zeros(nRows,nCols,nSlices,3); 

% Temp hack DM 4/5/16
ct_reg = zeros(256, 256,nSlices);
sublobar_reg = zeros(256, 256,nSlices);
proton_ct_fused = zeros(256, 256,nSlices,3); 


tic;
for nsl = 1:nSlices
    clc;
    fprintf('slice %d ...\n', nsl);
    % slice CT
    slice_ct = nii_data.nlz.ct.img(:,:,nsl);
    slice_ct_data = make_nii(slice_ct,[],[],64);
    nii_files.slice.ct = fullfile(output_folder,'slice_ct.nii.gz');
    save_nii(slice_ct_data, nii_files.slice.ct);
    
    % slice proton
    slice_proton = nii_data.nlz.proton.img(:,:,nsl);
    slice_proton_data = make_nii(slice_proton,[],[],64);
    nii_files.slice.proton = fullfile(output_folder,'slice_proton.nii.gz');
    save_nii(slice_proton_data, nii_files.slice.proton);
    
    %slice sublobarmask
        slice_sublobar = lungmask.ct_resz(:,:,nsl);
    slice_sublobar_data = make_nii(slice_sublobar,[],[],64);
    nii_files.slice.sublobar = fullfile(output_folder,'slice_sublobar_mask.nii.gz');
    save_nii(slice_sublobar_data, nii_files.slice.sublobar);

    
    
    system([script_folder,'/ct_proton_def ', nii_files.slice.ct, ' ', nii_files.slice.proton, ' ', moving_data2,' ',nii_files.slice.sublobar]);

    slice_reg_file = fullfile(output_folder,'slice_proton_fixed_ct_movingWarped.nii.gz');
    slice_reg = load_nii(slice_reg_file);
%     figure,imshow([slice_ct,slice_proton,slice_reg.img]);
%      c=imfuse(slice_proton,slice_ct);
%      creg=imfuse(slice_proton,slice_reg.img);
%      figure, imagesc([c,creg]);
    ct_reg(:,:,nsl) = slice_reg.img;
    proton_ct_fused(:,:,nsl,:) = imfuse(nii_data.nlz.proton.img(:,:,nsl),ct_reg(:,:,nsl));
    
    nii_files.slice.sublobar_reg = strrep(slice_reg_file,'Warped','mask_reg');
    sublobar_reg_sl = load_nii(nii_files.slice.sublobar_reg);
    sublobar_reg(:,:,nsl) = sublobar_reg_sl.img;
%     figure,imshow3D([slice_sublobar,sublobar_reg_sl.img],[0,10]);colormap('jet')
%         figure,imshow3D([slice_sublobar,sublobar_reg_sl.img],[100,110]);colormap('jet')
%             figure,imshow3D([slice_sublobar,sublobar_reg_sl.img],[130,170]);colormap('jet')
%                 figure,imshow3D([slice_sublobar,sublobar_reg_sl.img],[230,240]);colormap('jet')



end
total_reg_time = toc;

save registered;

%% registration alternative: 3D data stack with 2D in-plane deformation
%     system([script_folder,'/ct_proton_3D ', reg_pair2.moving, ' ', reg_pair2.fixed, ' ', moving_data2]);

%figure,imshow3D([nii_data.nlz.proton.img, ct_reg])


%% CONVERT BINARY PROTON LUNG MASK TO SUBLOBAR MASK
protonmask = load(proton_lungmask_matfile);


%binary_mask = protonmask.mask(:,:,end:-1:1); % NEED TO MAKE SURE THE SLICE ORDER IS CONSISTENT WITH PROTON
binary_mask = protonmask.lungsmask(:,:,:); % NEED TO MAKE SURE THE SLICE ORDER IS CONSISTENT WITH PROTON
proton_sublobar = merge_ctlobemask_into_lungmask(sublobar_reg,binary_mask);
figure,imshow3D([sublobar_reg,proton_sublobar])

%% defect analysis: measure defects on binary lung mask first simply because 
% defect analysis script currently does not recognize sublobar values.
temp_handles.he = double(nii_data.cor.he.img);
temp_handles.proton = double(nii_data.cor.proton.img);
temp_handles.lungsmask = binary_mask;
temp_handles.lobemask = binary_mask;
temp_handles.plotting_range = get_plotting_range(temp_handles.lungsmask);
temp_handles.hist_mu=385;
temp_handles.frangi_thresh=0.3;
temp_handles.kmeans_threshold = 1;
temp_handles.vessel_removal = true;
temp_handles.low_intn_perc_cutoff=1;
temp_handles.frangi_sigmas = 1;
temp_handles.helium_voxelsize = prod(image_headers.he.PixelSpacing)*image_headers.he.SliceThickness*1e-6;
%temp_handles.frangi_sigmas = 2.5;


%defect_mask_kmeans = detect_defects_by_kmeans_fangi(temp_handles);
defect_mask_kmeans = detect_defects_by_kmeans_sarp(temp_handles);


 [defect_mask_kmeans,~,~,vessel_mask,transformed_proton,low_intn_perc,ventilation_mask] = ventilation_defect_by_kmeans(temp_handles);

figure,plot_defect_mask(temp_handles.he,defect_mask_kmeans,temp_handles,temp_handles.plotting_range,studyid)

defect_sublobar = proton_sublobar.*defect_mask_kmeans;


% BASELINE DEFECTS
total_count = zeros(234, 1);
defect_count = zeros(234, 1);

cd(data_folder);
for ii = 1:nSlices
    
   
   defect_temp = defect_mask_kmeans(:, :, ii);
    
    % Get current sublobe mask
    sublobe_temp = proton_sublobar(:, :, ii);
    
    % Create combination of the two
    overlay_temp = zeros([256, 256]);
    overlay_temp(defect_temp > 0) = sublobe_temp(defect_temp > 0); 
    
    for jj = 1:256;
   
        for kk = 1:256;
        if(sublobe_temp(jj, kk) ~= 0) 
            total_count(sublobe_temp(jj, kk)) = total_count(sublobe_temp(jj, kk)) + 1;
        end
    
        if(overlay_temp(jj, kk) ~= 0)
            defect_count(overlay_temp(jj, kk)) = defect_count(overlay_temp(jj, kk)) + 1;
        end
        end   
    end
    
end

percent_defect = -1* ones(234, 1);
percent_defect(total_count ~= 0) = defect_count(total_count ~= 0) ./ total_count(total_count ~= 0);

fid = fopen([output_folder, '_results_baseline.csv'], 'w');
sublobe_lookup_init;
emptyCells = cellfun('isempty', sublobe_lookup);
fprintf(fid, 'studyid, region, def_voxel, tot_voxel, vdp\n');

LUL_defect = 0;
LUL_total = 0;
LLL_defect = 0;
LLL_total = 0;
RUL_defect = 0;
RUL_total = 0;
RML_defect = 0;
RML_total = 0;
RLL_defect = 0;
RLL_total = 0;
RL_defect = 0;
RL_total = 0;
LL_defect = 0;
LL_total = 0;
whole_defect = 0;
whole_total = 0;



for kk = 1:234
   if(emptyCells(kk) == 0 && percent_defect(kk) ~= -1)
      disp(['Region ', sublobe_lookup{kk}, ' is ', num2str(defect_count(kk)) ' / ', num2str(total_count(kk)), ' or ' num2str(percent_defect(kk)*100) , ' percent baseline defected']); 
      fprintf(fid, '%s, %s, %s, %s, %s\n', studyid, sublobe_lookup{kk}, num2str(defect_count(kk)), num2str(total_count(kk)), num2str(percent_defect(kk)*100));

      defect_count_temp = defect_count(kk);
      total_count_temp = total_count(kk); 
      %LUL
      if(ismember(kk, [1 2 3 4 5 6 7 8]));
         LUL_defect = LUL_defect + defect_count_temp;
         LUL_total = LUL_total + total_count_temp;
         LL_defect = LL_defect + defect_count_temp;
         LL_total = LL_total + total_count_temp;
      end
      
      %LLL
      if(ismember(kk, [102 104 105 106 16]));
         LLL_defect = LLL_defect + defect_count_temp;
         LLL_total = LLL_total + total_count_temp;
         LL_defect = LL_defect + defect_count_temp;
         LL_total = LL_total + total_count_temp;   
      end
      
      %RUL
      if(ismember(kk, [129 130 131 32]));
         RUL_defect = RUL_defect + defect_count_temp;
         RUL_total = RUL_total + total_count_temp;
         RL_defect = RL_defect + defect_count_temp;
         RL_total = RL_total + total_count_temp;   
      end
      
      %RML
      if(ismember(kk, [164 165 64]));
         RML_defect = RML_defect + defect_count_temp;
         RML_total = RML_total + total_count_temp;
         RL_defect = RL_defect + defect_count_temp;
         RL_total = RL_total + total_count_temp;   
      end
   
      %RLL
      if(ismember(kk, [230 231 232 233 234 128]));
         RLL_defect = RLL_defect + defect_count_temp;
         RLL_total = RLL_total + total_count_temp;
         RL_defect = RL_defect + defect_count_temp;
         RL_total = RL_total + total_count_temp;   
      end;
      
      % Whole lung
      whole_defect = whole_defect + defect_count_temp;
      whole_total = whole_total + total_count_temp;
   end
end

fprintf(fid, '%s, LUL, %s, %s, %s\n', studyid, num2str(LUL_defect), num2str(LUL_total), num2str(LUL_defect/LUL_total*100));
fprintf(fid, '%s, LLL, %s, %s, %s\n', studyid, num2str(LLL_defect), num2str(LLL_total), num2str(LLL_defect/LLL_total*100));
fprintf(fid, '%s, RUL, %s, %s, %s\n', studyid, num2str(RUL_defect), num2str(RUL_total), num2str(RUL_defect/RUL_total*100));
fprintf(fid, '%s, RML, %s, %s, %s\n', studyid, num2str(RML_defect), num2str(RML_total), num2str(RML_defect/RML_total*100));
fprintf(fid, '%s, RLL, %s, %s, %s\n', studyid, num2str(RLL_defect), num2str(RLL_total), num2str(RLL_defect/RLL_total*100));
fprintf(fid, '%s, LL, %s, %s, %s\n', studyid, num2str(LL_defect), num2str(LL_total), num2str(LL_defect/LL_total*100));
fprintf(fid, '%s, RL, %s, %s, %s\n', studyid, num2str(RL_defect), num2str(RL_total), num2str(RL_defect/RL_total*100));
fprintf(fid, '%s, whole, %s, %s, %s\n', studyid, num2str(whole_defect), num2str(whole_total), num2str(whole_defect/whole_total*100));


fclose(fid);


cd(data_folder);
save save_point_final; 

