function  [cleaned_mask,vessel_stack_binarized,low_intn_percent] = detect_defects_by_kmeans_sarp(temp_handles)
%
% Inputs: temp_handles.lungmask
%                     .he                   -acquired in axial
%                     .proton               
%                     .hist_mu              -used for stadardizing transform
%                     .frangi_thresh        -used for binarizing vessel mask
%                     .Kmeans_trehsold      -used for binarizing defect mask
%                     .vessel_removal       -used to include or exclude vessels
%                                            from defect mask
%                     
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
% Copyright W. Zha @2015
%%
debug =0;


[nRows,nCols,nSlices]=size(temp_handles.he);
mid_slice = ceil(nSlices/2);

lungsmask =temp_handles.lungsmask(:,:,1:nSlices);
valid_lung_slice_list = find(squeeze(sum(squeeze(sum(lungsmask,1)),1))>1);

slice_lung_voxels= squeeze(sum(sum(lungsmask,1),2));
slice_count_diff=find(slice_lung_voxels(1:end-1)-slice_lung_voxels(2:end)>1000);

slice_lung_voxels= squeeze(sum(sum(lungsmask,1),2));
diaphragm_slices=slice_count_diff(1):nSlices;
diaphragm_slices(diaphragm_slices<mid_slice)=[];
vasculature_slices=find(slice_lung_voxels>=.6*slice_lung_voxels(mid_slice));

se=strel('disk',2);
se4=strel('disk',4);
se2=strel('disk',6);
se3=strel('disk',3);

transformed_data=hist_transformation(temp_handles.proton,temp_handles.hist_mu);

%% 2D Frangifilter
vessel_stack=zeros(nRows,nCols,nSlices);
options.BlackWhite=false;
options.sigmas=temp_handles.sigmas;
shrunk_lungmask = zeros(nRows,nCols,nSlices);
dilated_lung=zeros(nRows,nCols,nSlices);
[outer_boundary,inner_boundary] = split_lung_boundaries_axial(temp_handles,1:nSlices);

for sl=valid_lung_slice_list%1:nSlices
    [~,~,~,filtered_at_sigmas]=FrangiFilter2D(transformed_data(:,:,sl),options);
    vessel_stack(:,:,sl)=filtered_at_sigmas;
    shrunk_lungmask(:,:,sl)=imerode(imdilate(temp_handles.lungsmask(:,:,sl),se4),se2);
    dilated_lung(:,:,sl)=imdilate(edge(outer_boundary(:,:,sl)),se3);
    
end

vessel_stack_binarized=vessel_stack;
vessel_stack_binarized(vessel_stack_binarized<temp_handles.frangi_thresh)=0;
vessel_stack_binarized(vessel_stack_binarized>0)=1;
vessel_stack_binarized=vessel_stack_binarized.*lungsmask;

vessel_mask_cleaned=vessel_stack_binarized-dilated_lung;
vessel_mask_cleaned(vessel_mask_cleaned<0)=0;
truncated_vessel_mask=vessel_mask_cleaned;


%% calculate percent low-intesity signals and decide the cluster option

masked_he_stack=temp_handles.he.*lungsmask;
pixel_list=masked_he_stack(lungsmask==1);
bin_count=hist(pixel_list);
low_intn_percent=bin_count(1)/sum(bin_count);


if low_intn_percent<.04
    cluster_option=1;
    num_clusters=5;
elseif low_intn_percent<0.1 && low_intn_percent>=0.04
    cluster_option=3;
    num_clusters=4;
elseif low_intn_percent>=0.1
    
    cluster_option=2;
    num_clusters=4;
end

%% cluster he3
cluster_stack=zeros(nRows,nCols,nSlices);
for n_def=valid_lung_slice_list
    cluster_stack(:,:,n_def)=ventCluster(temp_handles.he(:,:,n_def),lungsmask(:,:,n_def),num_clusters,cluster_option);
end
kmeans_cutoff=1;
VVmask=zeros(nRows,nCols,nSlices);
VVmask(cluster_stack(:)>kmeans_cutoff)=10;
VVmask(cluster_stack(:)<=kmeans_cutoff & cluster_stack(:)>0)=1;

defect_mask=zeros(size(cluster_stack));
defect_mask(VVmask(:)==1)=1;
cleaned_mask = defect_mask;

%%
if debug
    title_string = 'raw helium'; %#ok<UNRCH>
    figure,plot_defect_mask(temp_handles.he,temp_handles.lungsmask,temp_handles,temp_handles.plotting_range,title_string);

    title_string= 'contrast enhanced helium';
    figure,plot_defect_mask(transformed_data,temp_handles.lungsmask,temp_handles,temp_handles.plotting_range,title_string);
      
    title_string= 'clustered helium';
    figure,plot_defect_mask(cluster_stack,temp_handles.lungsmask,temp_handles,temp_handles.plotting_range,title_string);
  
    title_string= 'vessel mask';
    figure,plot_defect_mask(temp_handles.he,vessel_mask_cleaned,temp_handles,temp_handles.plotting_range,title_string);

    title_string= 'defect mask';
    figure,plot_defect_mask(temp_handles.he,defect_mask,temp_handles,temp_handles.plotting_range,title_string);
end
%% Morphometric Correction 
defect_slice_flag = zeros(nSlices,1);
defect_stat=struct;
defect_count = 0;
for sl=valid_lung_slice_list
    filled_lungsmask = imfill(lungsmask(:,:,sl),'holes');
        this_lung_edge = edge(filled_lungsmask);
        
        dilated_lung_mask = filled_lungsmask-imerode(filled_lungsmask,se);
    dilated_inner_boundary=imdilate(inner_boundary(:,:,sl),se3);
    if sum(diaphragm_slices==sl)==1
        dilated_lung_mask=dilated_lung_mask+dilated_inner_boundary;
        dilated_lung_mask(dilated_lung_mask>0)=1;
    end
    
        dilated_outer_boundary=imdilate(outer_boundary(:,:,sl),se4).*dilated_lung_mask;
    
    this_lobe_mask = double(temp_handles.lobemask(:,:,sl));
    this_lobe_mask(this_lobe_mask(:)==3 | this_lobe_mask(:)==4)=2;
    this_lobe_mask(this_lobe_mask(:)==6)=5;
    
        dilated_lung_mask_3pixels = imdilate(this_lung_edge,se3);
    this_defect_kmeans = cleaned_mask(:,:,sl);
    if ~temp_handles.vessel_removal
        defect_mask_shrunk = this_defect_kmeans-dilated_lung_mask;
    else
    defect_mask_shrunk = this_defect_kmeans-dilated_lung_mask-truncated_vessel_mask(:,:,sl);
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
        extent_struct = regionprops(cc,'Extent');
        extent = cat(1,extent_struct.Extent);
        num = length(unique(this_labeled_mask))-1;
        lung_bounary_struct = regionprops(cc,dilated_lung_mask_3pixels,'MeanIntensity');
        lung_boundary_intn = cat(1,lung_bounary_struct.MeanIntensity);
        
        for n_ob=1:num
            object_mask = zeros(nRows,nCols);
            object_mask(this_labeled_mask==n_ob)=n_ob;
            dilated_object_mask=imdilate(object_mask,se);
            if sum(vasculature_slices==sl)==1
                lost_object_edge_pixels=find((dilated_outer_boundary+dilated_object_mask)==(n_ob+1));
            else
                lost_object_edge_pixels=find((dilated_lung_mask+dilated_object_mask)==(n_ob+1));
            end
            if ~isempty(lost_object_edge_pixels)
                object_mask(lost_object_edge_pixels)=n_ob;
            end
            eroded_defect = imerode(object_mask,se);
            eroded_pixels = sum(eroded_defect(:)>0);
            %% penalize more on 1) inner lung boundaries where vessels may be attached to
            % and 2) tubular shaped vessels that were not caught up by
            % the vessel mask
            if sum(vasculature_slices==sl)==1 %vessels
                
                if (eccentricity(n_ob)>.85 && axisRatio(n_ob)>3 || eroded_pixels<15) && isempty(lost_object_edge_pixels)
                    this_labeled_mask(this_labeled_mask==n_ob)=0;
                    lost_object_edge_pixels=[];
                
                end
            end
            %% Correct for partial voluming-remove defects aligned along the lung edges.
            if (lung_boundary_intn(n_ob)>0.7 && extent(n_ob)<0.3)|| eroded_pixels==0
                this_labeled_mask(this_labeled_mask==n_ob)=0;
                lost_object_edge_pixels=[];
            end
            
            
            %% Clean up the false defects that occurred due to mismatch proton
            % and 3He due to the appearance of diaphragm in AXIAL view
            if sum(diaphragm_slices==sl)==1 && sum(this_labeled_mask(:)==n_ob)~=0
                object_mask(lost_object_edge_pixels)=n_ob;
                outer_boundary_pixels= (dilated_outer_boundary+dilated_object_mask)==(n_ob+1);
                if isempty(outer_boundary_pixels)
                    lobe_defect = double(this_lobe_mask).*object_mask;
                    
                    if sum(lobe_defect(:))~=0
                        overlapped_lobevalues=unique(lobe_defect(:))/n_ob;
                        
                        overlapped_lobevalues(overlapped_lobevalues==0)=[];
                        defect_lobe_ratio=zeros(length(overlapped_lobevalues),1);
                        
                        for n_v=1:length(overlapped_lobevalues)
                            lobe_no=overlapped_lobevalues(n_v);
                            defect_lobe_ratio(n_v)=sum(lobe_defect(:)==lobe_no*n_ob)/sum(this_lobe_mask(:)== lobe_no);
                        end
                        
                        if max(defect_lobe_ratio)>.15 && sum(object_mask(:))>400
                            this_labeled_mask(this_labeled_mask==n_ob & cleaned_mask(:,:,sl-1)==0)=0;
                            
                            
                          
                        end
                    end
                end
            end
           %% Compute defect statistics 
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


%% lobe defect separation
[temp_handles,lobe_label_values,lobe_label_list,valid_lobe_labels]=identify_lobes(temp_handles);
for n_lobe=1:length(valid_lobe_labels)
    this_lobe = sprintf('lobemask%s',upper(lobe_label_list{n_lobe}));
    if isfield(temp_handles,this_lobe)
        nSlice_lobe =size(temp_handles.(this_lobe),3);
        if nSlice_lobe>=nSlices
            this_lobe_mask=temp_handles.(this_lobe)(:,:,1:nSlices);
        else
            this_lobe_mask=zeros(size(defect_mask));
            this_lobe_mask(:,:,1:nSlice_lobe)=temp_handles.(this_lobe);
        end
        defect_mask(defect_mask(:)==1 &this_lobe_mask(:)==1)=lobe_label_values(n_lobe);
        cleaned_mask(cleaned_mask(:)>0 &this_lobe_mask(:)==1)=lobe_label_values(n_lobe);
    end
end
