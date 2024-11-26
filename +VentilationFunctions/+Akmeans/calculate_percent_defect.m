 function volume_struct = calculate_percent_defect(defect_mask,temp_handles)
[nTran,nSag,nCor]=size(temp_handles.he);
% helium_voxel_size=str2double(temp_handles.xres)*str2double(temp_handles.yres)*str2double(temp_handles.zres)*1e-6;

lungsmask =temp_handles.lungsmask(:,:,1:nCor);

[temp_handles,lobe_label_values,lobe_label_list]=VentilationFunctions.Akmeans.identify_lobes(temp_handles);
nLobes = length(lobe_label_values);
num_mask_slices = size(temp_handles.lobemask,3);
if num_mask_slices>=nCor
    lung_lobe_mask = double(temp_handles.lobemask(:,:,1:nCor));
else
    lung_lobe_mask=zeros(nTran,nSag,nCor);
    lung_lobe_mask(:,:,1:num_mask_slices)=double(temp_handles.lobemask);
end

%%
defect_mask=defect_mask(:,:,1:nCor).*lungsmask;
defect_mask(defect_mask>0)=1;
volume_struct.labels={'whole lung'};
volume_struct.labels(1+(1:nLobes))=lobe_label_list;
if ~isfield(temp_handles,'helium_voxelsize')
    temp_handles.helium_voxelsize = 1;
end
volume_struct.absolute.lung(1)=sum(lungsmask(:))*temp_handles.helium_voxelsize;
volume_struct.absolute.defect(1)=sum(defect_mask(:))*temp_handles.helium_voxelsize;
volume_struct.absolute.ventilation(1)=(sum(lungsmask(:))-sum(defect_mask(:)))*temp_handles.helium_voxelsize;

volume_struct.percent.defect(1) =100*sum(defect_mask(:))/sum(lungsmask(:));
volume_struct.percent.ventilation(1) = 100-volume_struct.percent.defect(1);
dmask_values=unique(defect_mask);
if max(dmask_values)==1
    
    defect_mask_lobe_labeled=defect_mask.*lung_lobe_mask;
else
    
    defect_mask_lobe_labeled=defect_mask;
end
if isfield(temp_handles,'vessel_mask')
    vessel_mask=temp_handles.vessel_mask.*lung_lobe_mask;
   volume_struct.absolute.vessel(1)=sum(temp_handles.vessel_mask(:))*temp_handles.helium_voxelsize; 
end

for n_lobe = 1:nLobes
    
    this_lobe_value=lobe_label_values(n_lobe);
    volume_struct.absolute.lung(1+n_lobe)=sum(lung_lobe_mask(:)==this_lobe_value)*temp_handles.helium_voxelsize;
    volume_struct.absolute.defect(1+n_lobe) =sum(defect_mask_lobe_labeled(:)==this_lobe_value)*temp_handles.helium_voxelsize;
    volume_struct.absolute.ventilation(1+n_lobe)=(sum(lung_lobe_mask(:)==this_lobe_value)- ...
        sum(defect_mask_lobe_labeled(:)==this_lobe_value))* ...
        temp_handles.helium_voxelsize;

    volume_struct.percent.ventilation(1+n_lobe) = 100* volume_struct.absolute.ventilation(1+n_lobe)/volume_struct.absolute.lung(1+n_lobe);
    volume_struct.percent.defect(1+n_lobe) = 100* volume_struct.absolute.defect(1+n_lobe)/volume_struct.absolute.lung(1+n_lobe);
     if isfield(temp_handles,'vessel_mask')
       
    volume_struct.absolute.vessel(1+n_lobe) =sum(vessel_mask(:)==this_lobe_value)*temp_handles.helium_voxelsize;
    volume_struct.percent.vessel(1+n_lobe) = 100* volume_struct.absolute.vessel(1+n_lobe)/volume_struct.absolute.lung(1+n_lobe);

    end
    
    fprintf('The percent defected in the %s is %6.2f%%.\n',lobe_label_list{n_lobe},volume_struct.percent.defect(1+n_lobe));
    
end

fprintf('The percent defected in the whole lung is %6.2f%%.\n',volume_struct.percent.defect(1) );

% lungmask_ = permute(lungsmask,[3,2,1]);
%% regional assessment: apical to base
% regional_assessment.tran = partition_volume_by_regions(lungsmask,defect_mask,'Axial');

%% regional assessment: anterior to posterior
% regional_assessment.cor = partition_volume_by_regions(lungsmask,defect_mask,'Coronal');
%% -------------------------------------------------------------
% function regional_assessment = partition_volume_by_regions(lungsmask,defectmask,direction)
% [nTran,nSag,nCor]=size(lungsmask);
% switch direction
%     case 'Coronal'
%         lungsmask_reformat=lungsmask;
%         defect_mask_reformat=defectmask;
%         nSlices=nCor;
%         region_mask=ones(nTran,nSag,nCor);
%         region_labels={'anterior','middle','posterior'};
%     case {'Axial','Transverse'}
%         lungsmask_reformat =reshape_coronal_view_to_axial(lungsmask);
%         defect_mask_reformat = reshape_coronal_view_to_axial(defectmask);
%         region_mask=ones(nSag,nCor,nTran);
%         region_labels={'apex','mid','base'};
%         
%         nSlices=nTran;
%     case 'sagital'
%         % no need so far
%     otherwise
%         error(calculate_percent_defect:partition_volume_by_regions,'Unknown direction: %s',direction);
% end
% nRegions=numel(region_labels);
% nLungPixels = sum(lungsmask(:));
% one_third = round(nLungPixels/3);
% two_thirds = nLungPixels-one_third;
% 
% slice_pixel_list = squeeze(sum(sum(lungsmask_reformat,1),2));
% slice_sum=zeros(nSlices,1);
% for n_sl=1:nSlices
%     slice_sum(n_sl)=sum(slice_pixel_list(1:n_sl));
% end
% [~,one_third_slice]=min(abs(slice_sum-one_third));
% [~,two_third_slice]=min(abs(slice_sum-two_thirds));
% region_mask(:,:,(one_third_slice+1):two_third_slice)=2*ones(size(region_mask,1),size(region_mask,2),two_third_slice-one_third_slice);
% region_mask(:,:,(two_third_slice+1):nSlices)=3*ones(size(region_mask,1),size(region_mask,2),nSlices-two_third_slice);
% 
% lungsmask_tran_value_list=lungsmask_reformat(:).*region_mask(:);
% defect_mask_value_list=defect_mask_reformat.*region_mask;
% labeled_defect_value_list=defect_mask_value_list(:);
% 
% ventilation_reformat=(lungsmask_reformat-defect_mask_value_list).*region_mask;
% labled_vent_value_list = ventilation_reformat(:);
% regional_assessment.labels{1}='whole';
% regional_assessment.labels((1:numel(region_labels))+1) =region_labels;
% regional_assessment.defect=zeros(nRegions+1,1);
% regional_assessment.ventilation=zeros(nRegions+1,1);
% regional_assessment.defect(1)=100*sum(labeled_defect_value_list>0)/sum(lungsmask_tran_value_list>0);
% regional_assessment.ventilation(1)=100*sum(labled_vent_value_list>0)/sum(lungsmask_tran_value_list>0);
% for n_reg=1:nRegions
%     regional_assessment.defect(n_reg+1)=100*sum(labeled_defect_value_list==n_reg)/sum(lungsmask_tran_value_list==n_reg);
%     regional_assessment.ventilation(n_reg+1)=100*sum(labled_vent_value_list==n_reg)/sum(lungsmask_tran_value_list==n_reg);
% end
%        
% %% -------------------------------------------------------------
% function axial_mask = reshape_coronal_view_to_axial(coronal_view)
% [nTran,nSag,nCor]=size(coronal_view);
% axial_mask=zeros(nSag,nCor,nTran);
% for n_tran=1:nTran
%     axial_mask(:,:,n_tran)=coronal_view(n_tran,:,:);
% end
    