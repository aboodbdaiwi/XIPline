function  [cleaned_mask,cluster_stack,proton_stack,...
    truncated_vessel_mask,transformed_data,low_intn_percent,ventilation_mask] = detect_defects_by_kmeans_fangi(temp_handles)
%%
% Inputs: temp_handles.lungmask2
%                     .he
%                     .proton
%                     .hist_mu              -used for stadardizing transform
%                     .frangi_thresh        -used for binarizing vessel mask
%                     .Kmeans_trehsold      -used for binarizing defect mask
%                     .vessel_removal       -used to include or exclude vessels
%                                            from defect mask
%                     .low_intn_perc_cutoff -used to swtich 2 clustering
%                                            options
%                     .num_clusters(=4)     -number of Kmeans clusters used
%                                            in the first round of clustering
% Features:
% 1) Vascular mask was estimated using 2D Frangi filter
% Frangi et al. Model-Based Quantitation of 3-D Magnetic Resonance
% Angiographic Images. IEEE Trans Med Imaging. 1999.
%
% 2) Adaptive defect calling: if the low-intensity signal over lungs is
% larger than the empirical cutoff value, option 2 (1st cluster was
% labeled as defect); otherwise, the 1st cluster was further clustered into
% 4 subclusters and the first 2 subclusters were determined as defect. The
% rationale is that highly defected lungs tend to have less signal-intensity
% variations, whereas less defected lungs tend to have more subtle signal
% heterogeneity requiring refined sub-classification.
%
% 3) Partial volume correction: the lung mask segmented from proton was
% eroded by 2 pixels to get the rought defect mask. If boundary defects
% detected on the eroded mask, such defects were dilated back towards the
% boundary to retore the defect.
%
% 4) Lobe separationg was performed on the defect mask if lobes were
% labeled.
%
%
% W. Zha @2015
%%
debug =0;
defect_mean_intn = [];
defect_std_intn=[];

[nRows,nCols,nSlices]=size(temp_handles.he);
mid_slice = ceil(nSlices/2);
vasculature_slices = mid_slice+(-3:3);
lungsmask =temp_handles.lungsmask(:,:,1:nSlices);
valid_lung_slice_list = find(squeeze(sum(squeeze(sum(lungsmask,1)),1))>1);

proton_stack = zeros(nRows,nCols,nSlices);
se5=strel('disk',5);

se=strel('disk',2);
se4=strel('disk',4);
se6=strel('disk',6);
if isfield(temp_handles,'registered_proton')
    transformed_data=hist_transformation(temp_handles.registered_proton.proton_reg,temp_handles.hist_mu);
else
    transformed_data=hist_transformation(temp_handles.proton,temp_handles.hist_mu);
end
if isfield(temp_handles,'proton_transform') && ~temp_handles.proton_trasform
    frangi_proton=temp_handles.proton; %fiesta
else
    frangi_proton=transformed_data;%ssfse


end
num_clusters = 4;
%% 2D Frangifilter
vessel_stack=zeros(nRows,nCols,nSlices);
options.BlackWhite=false;
options.sigmas=temp_handles.frangi_sigmas;
shrunk_lungmask = zeros(nRows,nCols,nSlices);
outter_boundary = split_lung_boundaries(temp_handles,1:nSlices);
dilated_lung=zeros(nRows,nCols,nSlices);


for sl=1:nSlices
    
%     [~,~,~,filtered_at_sigmas]=FrangiFilter2D(transformed_data(:,:,sl),options);
    [~,~,~,filtered_at_sigmas]=FrangiFilter2D(frangi_proton(:,:,sl),options);

    vessel_stack(:,:,sl)=filtered_at_sigmas;
    shrunk_lungmask(:,:,sl)=imerode(imdilate(temp_handles.lungsmask(:,:,sl),se4),se6);
        dilated_lung(:,:,sl)=imdilate(edge(outter_boundary(:,:,sl)),se5);

end

vessel_stack_binarized=vessel_stack;
vessel_stack_binarized(vessel_stack_binarized<temp_handles.frangi_thresh)=0;
% 
vessel_stack_binarized(vessel_stack_binarized>0)=1;
vessel_stack_binarized=vessel_stack_binarized.*lungsmask;


%%
vessel_mask_cleaned=vessel_stack_binarized-dilated_lung;
vessel_mask_cleaned(vessel_mask_cleaned<0)=0;

% vessel_mask_cleaned = vessel_mask_cleanup(temp_handles,vessel_stack_binarized,1:nSlices);
% truncated_vessel_mask=zeros(nRows,nCols,nSlices);
% truncated_vessel_mask(:,:,vasculature_slices)=vessel_mask_cleaned(:,:,vasculature_slices);
truncated_vessel_mask=vessel_mask_cleaned;
truncated_vessel_mask= bwareaopen(truncated_vessel_mask,temp_handles.vessel_thresh);

%%
% lungsmask_corrected=lungsmask-truncated_vessel_mask;
lungsmask_corrected=lungsmask;
masked_he_stack=temp_handles.he.*lungsmask_corrected;
pixel_list=masked_he_stack(lungsmask==1);
bin_count=hist(pixel_list);
low_intn_percent=bin_count(1)/sum(bin_count);
if low_intn_percent<.04 % the relative number of low intensity pixels decides the size of the 1st cluster of Kmeans
    cluster_option=1;
    num_clusters=5;
elseif low_intn_percent<.1 && low_intn_percent>=0.04
    cluster_option=3;
elseif low_intn_percent>=0.1
    cluster_option=2;%equivalent to call kmeans with 4 clusters
end
cluster_stack=zeros(nRows,nCols,nSlices);
for n_sl=1:nSlices
    cluster_stack(:,:,n_sl)=ventCluster(temp_handles.he(:,:,n_sl),lungsmask_corrected(:,:,n_sl),num_clusters,cluster_option);
end
% cluster_stack=ventCluster2(temp_handles.he,lungsmask_corrected,num_clusters,cluster_option);
kmeans_cutoff=1;
VVmask=zeros(nRows,nCols,nSlices);
VVmask(cluster_stack(:)>kmeans_cutoff)=10;
VVmask(cluster_stack(:)<=kmeans_cutoff & cluster_stack(:)>0)=1;

defect_mask=zeros(size(cluster_stack));
defect_mask(VVmask(:)==1)=1;
cleaned_mask = defect_mask;

if debug
    temp_handles.plotting_range=get_plotting_range(temp_handles.lungsmask);
    title_string = 'raw helium';
    figure,plot_defect_mask(temp_handles.he,temp_handles.lungsmask,temp_handles,temp_handles.plotting_range,title_string);
    title_string = 'raw proton';
    figure,plot_defect_mask(temp_handles.registered_proton.proton_reg,zeros(nRows,nCols,nSlices),temp_handles,temp_handles.plotting_range,title_string);
    figure,plot_defect_mask(temp_handles.proton,temp_handles.lungsmask,temp_handles,temp_handles.plotting_range,title_string);
    title_string= 'contrast enhanced helium';
    figure,plot_defect_mask(he_stack_eq,temp_handles.lungsmask,temp_handles,temp_handles.plotting_range,title_string);
    
    title_string= 'clustered proton';
    figure,plot_defect_mask(proton_stack,temp_handles.lungsmask,temp_handles,temp_handles.plotting_range,title_string);
    
    title_string= 'clustered helium';
    figure,plot_defect_mask(cluster_stack,temp_handles.lungsmask,temp_handles,temp_handles.plotting_range,title_string);
    title_string= '  vessel mask';
    figure,plot_defect_mask(vessel_stack_binarized+vessel_mask_cleaned,zeros(nRows,nCols,nSlices),temp_handles,temp_handles.plotting_range,title_string);
    figure,plot_defect_mask(temp_handles.proton,vessel_mask_cleaned,temp_handles,temp_handles.plotting_range,title_string);
        figure,plot_defect_mask(transformed_data,truncated_vessel_mask,temp_handles,temp_handles.plotting_range,title_string);
                figure,plot_defect_mask(transformed_data,vessel_stack_binarized,temp_handles,temp_handles.plotting_range,title_string);
figure,plot_defect_mask(temp_handles.he,vessel_mask_cleaned,temp_handles,temp_handles.plotting_range,title_string);
title_string='standardized proton';
    figure,plot_defect_mask(transformed_data,zeros(nRows,nCols,nSlices),temp_handles,temp_handles.plotting_range,title_string);
    title_string= 'defect mask';
    figure,plot_defect_mask(cluster_stack,defect_mask,temp_handles,temp_handles.plotting_range,title_string);
end
defect_slice_flag = zeros(nSlices,1);
defect_stat=struct;
defect_count = 0;
for sl=valid_lung_slice_list
    this_lung_edge = edge(imfill(lungsmask(:,:,sl),'holes'));
    dilated_lung_mask = imdilate(this_lung_edge,se);
        dilated_outter_boundary=imdilate(outter_boundary(:,:,sl),se4).*dilated_lung_mask;

    this_lobe_mask = temp_handles.lobemask(:,:,sl);
    RUL = sum(this_lobe_mask(:)==2);
    LUL = sum(this_lobe_mask(:)==5);
    
    this_defect_kmeans = cleaned_mask(:,:,sl);
    if temp_handles.vessel_removal
    defect_mask_shrunk = this_defect_kmeans-dilated_lung_mask-truncated_vessel_mask(:,:,sl);
    else
        defect_mask_shrunk = this_defect_kmeans-dilated_lung_mask;
    end
    defect_mask_shrunk(defect_mask_shrunk<0)=0;
    defect_mask_shrunk_cleaned=defect_mask_shrunk;
    
    cc=bwconncomp(defect_mask_shrunk_cleaned,8);
    if cc.NumObjects==0
        this_labeled_mask=zeros(nRows,nCols);
    else
        this_labeled_mask=double(labelmatrix(cc));
        major_length = regionprops(cc,'MajorAxisLength');
        minor_length = regionprops(cc,'MinorAxisLength');
        axisLength=[cat(1,major_length.MajorAxisLength),cat(1,minor_length.MinorAxisLength)];
        axisRatio = axisLength(:,1)./axisLength(:,2);
        eccentricity_struct = regionprops(cc,'Eccentricity');
        eccentricity=cat(1,eccentricity_struct.Eccentricity);
        num = length(unique(this_labeled_mask))-1;
        
        for n_ob=1:num
            object_mask = zeros(nRows,nCols);
            object_mask(this_labeled_mask==n_ob)=n_ob;
            dilated_object_mask=imdilate(object_mask,se);
             if sum(vasculature_slices==sl)==1
                lost_object_edge_pixels=find((dilated_outter_boundary+dilated_object_mask)==(n_ob+1));
            else
                lost_object_edge_pixels=find((dilated_lung_mask+dilated_object_mask)==(n_ob+1));
             end
            if ~isempty(lost_object_edge_pixels)
                object_mask(lost_object_edge_pixels)=n_ob;
            end
            eroded_defect = imerode(object_mask,se);
            eroded_pixels = sum(eroded_defect(:)>0);
%             lost_object_edge_pixels=find((dilated_object_mask+dilated_lung_mask)==(n_ob+1));
%             object_mask_tmp=object_mask;
%             if ~isempty(lost_object_edge_pixels) && sum(vasculature_slices==sl)==1%2/10/2015
%                 object_mask_tmp(lost_object_edge_pixels)=n_ob;
%                 eroded_defect = imerode(object_mask_tmp,se);
%             else
%                 eroded_defect = imerode(object_mask_tmp,se);
%                 
%             end
%             eroded_pixels = sum(eroded_defect(:)>0);
            
            
            if (eccentricity(n_ob)>.85 && axisRatio(n_ob)>3 && isempty(lost_object_edge_pixels))...
                    || eroded_pixels==0
                
                this_labeled_mask(this_labeled_mask==n_ob)=0;
                lost_object_edge_pixels=[];
                
            end
            
            
            
            this_labeled_mask(lost_object_edge_pixels)=n_ob;

            
            if sum(this_labeled_mask(:)==n_ob)~=0
                defect_count=defect_count+1;
                defect_stat(defect_count).slice_no=sl;
                defect_stat(defect_count).defect_mask = zeros(nRows,nCols);
                defect_stat(defect_count).defect_mask(this_labeled_mask==n_ob) = 1;
                this_defect=bwconncomp(defect_stat(defect_count).defect_mask,8);
                property_list = {'Area','Eccentricity','MajorAxisLength','MinorAxisLength','Extent'};
                
                for n_p=1:numel(property_list)
                    this_property = regionprops(this_defect,property_list{n_p});
                    
                    if n_p==1 && this_defect.NumObjects>1
                        labeled_mask=bwlabel(defect_stat(defect_count).defect_mask);
                        defect_stat(defect_count).Area=cat(1,this_property.Area);
                        [~,max_ind]=max(defect_stat(defect_count).Area);
                        defect_stat(defect_count).defect_mask(labeled_mask~=max_ind)=0;
                        this_labeled_mask(this_labeled_mask==n_ob)=0;
                        this_labeled_mask(labeled_mask==max_ind)=n_ob;
                        this_defect=bwconncomp(defect_stat(defect_count).defect_mask,8);
                        this_property = regionprops(this_defect,property_list{n_p});
                        
                    end
                    defect_stat(defect_count).(property_list{n_p})=cat(1,this_property.(property_list{n_p}));
                end
                
                this_property = regionprops(this_defect,temp_handles.he(:,:,sl),'MeanIntensity');
                defect_stat(defect_count).MeanIntensity = cat(1,this_property.MeanIntensity);
                this_property = regionprops(this_defect,vessel_mask_cleaned(:,:,sl),'MeanIntensity');
                defect_stat(defect_count).Mean_vessel_cluster = cat(1,this_property.MeanIntensity);
            end
        end
    end
    this_labeled_mask(this_labeled_mask>0)=1;
    cleaned_mask(:,:,sl)=this_labeled_mask;
    if sum(this_labeled_mask(:))~=0
        defect_slice_flag(sl) = 1;
    end
end
% ventilation map
ventilation_mask = cluster_stack;
ventilation_mask(truncated_vessel_mask==1)=0;% vessels
ventilation_mask(isnan(ventilation_mask))=0;
vent_mask = ventilation_mask>1;

vent_stack = kmeans_clustering(vent_mask.*temp_handles.he,5);
vent_stack = vent_stack+1;
vent_stack(vent_stack == 1) = 0 ;
ventilation_mask(ventilation_mask>1)=vent_stack(vent_stack>0);
%% lobe defect separation
temp_handles=identify_lobes(temp_handles);
% for n_lobe=1:length(valid_lobe_labels)
%     this_lobe = sprintf('lobemask%s',upper(lobe_label_list{n_lobe}));
%     if isfield(temp_handles,this_lobe)
%         nSlice_lobe =size(temp_handles.(this_lobe),3);
%         if nSlice_lobe>=nSlices
%             this_lobe_mask=temp_handles.(this_lobe)(:,:,1:nSlices);
%         else
%             this_lobe_mask=zeros(size(defect_mask));
%             this_lobe_mask(:,:,1:nSlice_lobe)=temp_handles.(this_lobe);
%         end
%         defect_mask(defect_mask(:)==1 &this_lobe_mask(:)==1)=lobe_label_values(n_lobe);
%         cleaned_mask(cleaned_mask(:)>0 &this_lobe_mask(:)==1)=lobe_label_values(n_lobe);
%     end
% end
cleaned_mask = cleaned_mask.*temp_handles.lobemask;