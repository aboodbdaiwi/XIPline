function  [cleaned_mask,cluster_stack,proton_stack,...
    truncated_vessel_mask,transformed_data,low_intn_percent,ventilation_mask,lungsmask_corrected,frangi_thresh] = ventilation_defect_by_kmeans(temp_handles)
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
he_ori = temp_handles.he;
he = temp_handles.he;

[nRows,nCols,nSlices] = size(he);
he_flat=zeros(nRows,nCols,nSlices);
for n_sl=1:nSlices
    slice_he = he(:,:,n_sl);
    maxhe = max(slice_he(:));
    minhe = min(slice_he(:));
    he_flat(:,:,n_sl) = (slice_he-minhe)*255/(maxhe-minhe);
end
lungsmask = temp_handles.lungsmask(:,:,1:nSlices);
contoured_slices = find(squeeze(sum(sum(lungsmask,1),2))>0).';
if ~isfield(temp_handles,'view')
    temp_handles.view = 'axial';
end
switch temp_handles.view
    case 'axial'
        diaphragm_sl = locate_diaphragm_slice(lungsmask);
        vasculature_slices = contoured_slices(2): min(diaphragm_sl+2,contoured_slices(end));
        diaphragm_slices = diaphragm_sl:contoured_slices(end);
    case 'coronal'
        
        vasculature_slices = contoured_slices(3:(end-2));
        diaphragm_slices = 0;
        lungsmask = fill_lungsmask(lungsmask,temp_handles.view);
    otherwise
        error('ventilation_defect_by_kmeans:view','%s is an unknown view type.',view);
        
end
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
    frangi_proton = temp_handles.proton; %fiesta
else
    frangi_proton = transformed_data;%ssfse
    
    
end
num_clusters = 4;
%% 2D Frangifilter
vessel_stack=zeros(nRows,nCols,nSlices);
options.BlackWhite=false;
options.sigmas=temp_handles.frangi_sigmas;
% options.FrangiScaleRange=[1 1];
% options.FrangiBeta =0.4;
% options.FrangiAlpha = 0.5;
% options.FrangiC = 1500;
% options.FrangiScaleRatio =1;
shrunk_lungmask = zeros(nRows,nCols,nSlices);
if strcmp(temp_handles.view,'coronal')
    outter_boundary = split_lung_boundaries_coronal(temp_handles.lungsmask,1:nSlices);
else
    outter_boundary = split_lung_boundaries(temp_handles,1:nSlices);
end

dilated_lung=zeros(nRows,nCols,nSlices);
%% FRANGI 2d
if ~isfield(temp_handles,'vessel_mask')
    for sl = vasculature_slices
        
        try
            [~,~,~,filtered_at_sigmas]=FrangiFilter2D(frangi_proton(:,:,sl),options);
        catch
            continue;
        end
        vessel_stack(:,:,sl)=filtered_at_sigmas;
        shrunk_lungmask(:,:,sl)=imerode(imdilate(lungsmask(:,:,sl),se4),se6);
        dilated_lung(:,:,sl)=imdilate(edge(outter_boundary(:,:,sl)),se5);
        
    end
    %%
    if isfield(temp_handles,'frangi_thresh')
        thresh_range = temp_handles.frangi_thresh;
    else
    thresh_range = 0.02:0.01:max(vessel_stack(:));
    end
    vessel_percent = zeros(length(thresh_range),1);
    truncated_vessel_mask_tmp = zeros(nRows, nCols,nSlices,length(thresh_range));
    proton = temp_handles.proton;
    vessel_thresh=temp_handles.vessel_thresh;
    for n_thresh = 1:length(thresh_range)
        vessel_stack_binarized=vessel_stack;
        
        this_thresh = thresh_range(n_thresh);
        %     vessel_stack_binarized(vessel_stack_binarized<temp_handles.frangi_thresh)=0;
        vessel_stack_binarized(vessel_stack_binarized<this_thresh)=0;
        %
        vessel_stack_binarized(vessel_stack_binarized>0)=1;
        % vessel_stack_binarized=vessel_stack_binarized.*lungsmask;
        vessel_stack_binarized=vessel_stack_binarized.*shrunk_lungmask;
        % vessel_stack_binarized=VesselnessFilter_mu(lungsmask,transformed_data,1);
        
        %
        vessel_mask_cleaned=vessel_stack_binarized-dilated_lung;
        vessel_mask_cleaned(vessel_mask_cleaned<0)=0;
        
        truncated_vessel_mask_tmp(:,:,:,n_thresh)=vessel_mask_cleaned;
        clustered_proton = kmeans_clustering(proton,4);
        clustered_proton(clustered_proton>1)=0;
        clustered_proton = clustered_proton.*lungsmask;
        truncated_vessel_mask_tmp(:,:,:,n_thresh) = truncated_vessel_mask_tmp(:,:,:,n_thresh).*clustered_proton;
        threshed_vessel =truncated_vessel_mask_tmp(:,:,:,n_thresh);
        
        for  sl = vasculature_slices
            threshed_vessel(:,:,sl)= bwareaopen(threshed_vessel(:,:,sl),vessel_thresh);
        end
        truncated_vessel_mask_tmp(:,:,:,n_thresh)=threshed_vessel;
        vessel_percent(n_thresh) = sum(threshed_vessel(:))*100/sum(lungsmask(:));
        %         if vessel_percent(n_thresh)>5 && vessel_percent(n_thresh)<6
        %             truncated_vessel_mask = truncated_vessel_mask_tmp(:,:,:,n_thresh);
        %             fprintf('frangi threshold = %6.2f.\n',this_thresh);
        %
        % %             break;
        %
        %         end
    end
    %     if ~exist('truncated_vessel_mask','var')
    [~,min_index] = min(abs(vessel_percent-5.4));
    
%     [sorted_percent,sorted_order] = sort(abs(vessel_percent-5.4),'descend');
%     sorted_percent(sorted_percent<0)=Inf;
%     [~,min_index] = min(sorted_percent);
    frangi_thresh = thresh_range(min_index);
    fprintf('frangi threshold = %6.2f.\n',frangi_thresh);
%         fprintf('vessel percent = %6.2f.\n',min_percent);

    truncated_vessel_mask =  truncated_vessel_mask_tmp(:,:,:,min_index);
    %     end
    %%
else
    truncated_vessel_mask = temp_handles.vessel_mask;
end

%%
lowerhalf = zeros(nRows,nCols,nSlices);
for n =  vasculature_slices(3:end)%vasculature_slices(3:(end-2))
    
    slicerange = get_plotting_range(lungsmask(:,:,n),-15);
    if ~isempty(slicerange.x)
        regionmask = zeros(nRows,nCols);
        %     regionmask(slicerange.y,slicerange.x)=1;
        %     regionmask(slicerange.y,round(mean(slicerange.x)-10):end)=1;
        regionmask(slicerange.y(15:(end-10)),slicerange.x)=1;
        
        
        lowerhalf(:,:,n) = regionmask;
    end
end

%% ADAPTIVE K-MEANS: determine clustering option
% lungsmask_corrected=lungsmask;
if ~isfield(temp_handles,'smooth_mask_flag')
    temp_handles.smooth_mask_flag = false;
end

if temp_handles.smooth_mask_flag
    lungsmask_corrected = smooth_lungmask(lungsmask,'coronal');
else
    lungsmask_corrected = lungsmask;
end
masked_he_stack=he.*lungsmask_corrected;
pixel_list=masked_he_stack(lungsmask_corrected==1);
bin_count = hist(pixel_list);
low_intn_percent=bin_count(1)/sum(bin_count);
lower_cutoff = 0.035;%0.04;%
upper_cutoff = 0.08;%0.1;
if low_intn_percent<lower_cutoff % the relative number of low intensity pixels decides the size of the 1st cluster of Kmeans
    cluster_option=1;
    num_clusters=5;
elseif low_intn_percent<upper_cutoff && low_intn_percent>=lower_cutoff
    cluster_option=3;
elseif low_intn_percent>=upper_cutoff
    cluster_option=2;%equivalent to call kmeans with 4 clusters
end
% cluster_stack = ventCluster2(he,lungsmask_corrected,num_clusters,cluster_option);
cluster_stack=zeros(nRows,nCols,nSlices);
%     cluster_stack = ventCluster2(he,lungsmask_corrected,num_clusters,cluster_option);


for n_sl=1:nSlices
%     cluster_stack(:,:,n_sl) = ventCluster2(he(:,:,n_sl),lungsmask_corrected(:,:,n_sl),num_clusters,cluster_option);
    cluster_stack(:,:,n_sl) = ventCluster_masked(he(:,:,n_sl),lungsmask_corrected(:,:,n_sl),num_clusters,cluster_option);
    %     cluster_stack(:,:,n_sl) = adaptiveKmeans(he(:,:,n_sl),lungsmask_corrected(:,:,n_sl),num_clusters,cluster_option);
    
end
kmeans_cutoff=1;
VVmask=zeros(nRows,nCols,nSlices);
VVmask(cluster_stack(:)>kmeans_cutoff)=10;
VVmask(cluster_stack(:)<=kmeans_cutoff & cluster_stack(:)>0)=1;

defect_mask=zeros(size(cluster_stack));
defect_mask(VVmask(:)==1)=1;

cleaned_mask = defect_mask;

if debug
    %     temp_handles.plotting_range=get_plotting_range(lungsmask);
    %     title_string = 'raw helium';
    figure,plot_defect_mask(he,lungsmask,temp_handles,temp_handles.plotting_range,title_string);
    title_string = 'raw proton';
    figure,plot_defect_mask(temp_handles.registered_proton.proton_reg,zeros(nRows,nCols,nSlices),temp_handles,temp_handles.plotting_range,title_string);
    figure,plot_defect_mask(temp_handles.proton,lungsmask,temp_handles,temp_handles.plotting_range,title_string);
    title_string= 'contrast enhanced helium';
    figure,plot_defect_mask(he_stack_eq,lungsmask,temp_handles,temp_handles.plotting_range,title_string);
    
    title_string= 'clustered proton';
    figure,plot_defect_mask(proton_stack,lungsmask,temp_handles,temp_handles.plotting_range,title_string);
    
    title_string= 'clustered helium';
    figure,plot_defect_mask(cluster_stack,lungsmask,temp_handles,temp_handles.plotting_range,title_string);
    title_string= '  vessel mask';
    figure,plot_defect_mask(vessel_stack_binarized+vessel_mask_cleaned,zeros(nRows,nCols,nSlices),temp_handles,temp_handles.plotting_range,title_string);
    figure,plot_defect_mask(temp_handles.proton,vessel_mask_cleaned,temp_handles,temp_handles.plotting_range,title_string);
    figure,plot_defect_mask(transformed_data,truncated_vessel_mask,temp_handles,temp_handles.plotting_range,title_string);
    figure,plot_defect_mask(transformed_data,vessel_stack_binarized,temp_handles,temp_handles.plotting_range,title_string);
    figure,plot_defect_mask(he,vessel_mask_cleaned,temp_handles,temp_handles.plotting_range,title_string);
    title_string='standardized proton';
    figure,plot_defect_mask(transformed_data,zeros(nRows,nCols,nSlices),temp_handles,temp_handles.plotting_range,title_string);
    title_string= 'defect mask';
    figure,plot_defect_mask(cluster_stack,defect_mask,temp_handles,temp_handles.plotting_range,title_string);
end
defect_slice_flag = zeros(nSlices,1);
se3 = strel('disk',3);
se2=strel('disk',2);
%% ADAPTIVE K-MEANS: morphometric correction of partial voluming

for sl=valid_lung_slice_list
    %     this_lung_edge = lungsmask(:,:,sl)-imerode(lungsmask(:,:,sl),se1);
    %     dilated_lung_mask = imdilate(this_lung_edge,se);
    dilated_lung_mask = lungsmask(:,:,sl)-imerode(lungsmask(:,:,sl),se2);
    dilated_outter_boundary=imdilate(outter_boundary(:,:,sl),se4).*dilated_lung_mask;
    
    
    this_lobe_mask = double(temp_handles.lobemask(:,:,sl));
    
        dilated_lung_mask_3pixels = lungsmask(:,:,sl)-imerode(lungsmask(:,:,sl),se3);%imdilate(this_lungT_edge,se3);
    dilated_lung_mask_2pixels = lungsmask(:,:,sl)-imerode(lungsmask(:,:,sl),se2);%imdilate(this_lung_edge,se3);
    
    this_defect_kmeans = cleaned_mask(:,:,sl);
    
    if temp_handles.vessel_removal
        defect_mask_shrunk = this_defect_kmeans-dilated_lung_mask-truncated_vessel_mask(:,:,sl);
    else
        % For obstructive lung disease with persistent focal defect (e.g.
        % CF), vessels are likely merged into large defects. vessel removal
        % may not be necessary.
        defect_mask_shrunk = this_defect_kmeans-dilated_lung_mask;
        truncated_vessel_mask = zeros(nRows,nCols,nSlices);
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
                lost_object_edge_pixels=find((dilated_outter_boundary+dilated_object_mask)==(n_ob+1));
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
                
                if (eccentricity(n_ob)>.85 && axisRatio(n_ob)>3 && eroded_pixels<15) && isempty(lost_object_edge_pixels)% changed 20 back to 15 of 0803
                    
                    this_labeled_mask(this_labeled_mask==n_ob)=0;
                    lost_object_edge_pixels=[];
                    
                end
            end
            %% Correct for partial voluming-remove defects aligned along the lung edges.
            %             boundary_overlap = dilated_lung_mask_3pixels.*object_mask;
            boundary_overlap = dilated_lung_mask_2pixels.*object_mask;
            vessel_overlap = object_mask.*lowerhalf(:,:,sl);
            
            if eroded_pixels==0 || (sum(boundary_overlap(:))/sum(object_mask(:))>0.6 && eroded_pixels<10)...
                    || (sum(vessel_overlap(:))>0.9*sum(object_mask(:)) )% deleted 0.9* sum(object_mask(:)) 1/5/17
                %             if (lung_boundary_intn(n_ob)>0.7 && extent(n_ob)<0.3)|| eroded_pixels==0
                this_labeled_mask(this_labeled_mask==n_ob)=0;
                lost_object_edge_pixels=[];
            end
            %% Clean up the false defects that occurred due to mismatch proton
            % and 3He due to the appearance of diaphragm
            if ~isempty(lost_object_edge_pixels)
                dilated_defect = imdilate(object_mask,se3);
                lost_object_edge_pixels=unique([lost_object_edge_pixels;find((this_defect_kmeans+dilated_defect)==(n_ob+1))]);
            end
            this_labeled_mask(lost_object_edge_pixels)=n_ob;
            if sum(diaphragm_slices==sl)==1
                object_mask(lost_object_edge_pixels)=n_ob;
                lobe_defect = double(this_lobe_mask).*object_mask;
                overlapped_lobevalues=unique(lobe_defect)/n_ob;
                overlapped_lobevalues(overlapped_lobevalues==0)=[];
                defect_lobe_ratio=zeros(length(overlapped_lobevalues),1);
                for n_v=1:length(overlapped_lobevalues)
                    lobe_no=overlapped_lobevalues(n_v);
                    defect_lobe_ratio(n_v)=sum(lobe_defect(:)==lobe_no*n_ob)/sum(this_lobe_mask(:)== lobe_no);
                end
                
                %                 if ~isempty(defect_lobe_ratio) && max(defect_lobe_ratio)>.15
                %                     this_labeled_mask(this_labeled_mask==n_ob & cleaned_mask(:,:,sl-1)==0)=0;
                %
                %
                %                 end
            end
            
            
        end
    end
    this_labeled_mask(this_labeled_mask>0)=1;
    cleaned_mask(:,:,sl)=this_labeled_mask;
    if sum(this_labeled_mask(:))~=0
        defect_slice_flag(sl) = 1;
    end
end
cleaned_mask = cleaned_mask.*logical(temp_handles.lungsmask);
%% VENTILATION MAP
fprintf('low_intn_perc = %6.3f\n',low_intn_percent);
% cluster_stack2 = ventCluster2(he_ori,lungsmask_corrected,num_clusters,cluster_option);

% vent_mask = lungsmask;

if isfield(temp_handles,'ventmaskbinary')
    bounded_mask = logical(temp_handles.ventmaskbinary);
    
    vent_mask = temp_handles.ventmaskbinary;
else
    vent_mask = temp_handles.lungsmask;
    bounded_mask = logical(temp_handles.lungsmask);
    
end

vent_mask(truncated_vessel_mask==1) = 0;% vessels
vent_mask(cleaned_mask==1) = 0; %defect
% vent_mask = vent_mask>1;
% lung_intensities = he_ori(lungsmask(:)==1);

he_nlz = (he_ori-min(he_ori(:)))./(max(he_ori(:))-min(he_ori(:)));
VentLevels = 5;
vent_stack = kmeans_clustering(vent_mask.*he_nlz,VentLevels);

cl4 = zeros(nRows,nCols,nSlices);
cl4(vent_stack==4)=1;
cl4_stack = kmeans_clustering(cl4.*he_nlz,VentLevels);


cl2 = zeros(nRows,nCols,nSlices);
cl2(vent_stack==2)=1;
% cl2(vent_stack==1)=1; %tcon
cl2_stack = kmeans_clustering(cl2.*he_nlz,VentLevels);



%% re-classify vent_stack into 3 ventilated bins
vent_stack = vent_stack+1; % SAVE CLUSTER # 1 FOR DEFECT
vent_stack(vent_stack == 1) = 0;
% vent_stack(cl2_stack>=1 & cl2_stack<=2)=2; % LOW VENTILATED REGION
% 7.21.2016 removed
% if low_intn_percent<0.1 %tcon
%     vent_stack(cl2_stack>2) = 3;
% end
if low_intn_percent<upper_cutoff
    vent_stack(cl2_stack>0 & cl2_stack<=2) =2;
else
    vent_stack(vent_stack==3)=2;
end

vent_stack(vent_stack == 4) = 3;
if low_intn_percent>=upper_cutoff
    
    vent_stack(vent_stack ==5) = 3;% MODERATE VENTILATED REGION
    %
    
    vent_stack(vent_stack==6)=4;%10/22
else
    %     vent_stack(vent_stack==5)=4; %tcon
    vent_stack(vent_stack == 5) = 3;
    vent_stack(cl4_stack>=4)=4;
    vent_stack(vent_stack == 6) = 4; % HYPER VENTILATED REGION
    
end

ventilation_mask = vent_stack.*bounded_mask;
ventilation_mask (cleaned_mask==1 & bounded_mask==1)= 1;
% if sum(truncated_vessel_mask(:))>0
%     ventilation_mask (cluster_stack==1 & truncated_vessel_mask==0 & cleaned_mask==0)= 2;
% else
%     ventilation_mask (cluster_stack==1 & cleaned_mask==0)= 2;
% end
ventilation_mask=ventilation_mask.*bounded_mask;

%  vent_stack(cl4_stack<4 & cl4_stack>=1)=4;
% if isfield(temp_handles,'ventmaskbinary')
%     ventilation_mask = vent_stack.*temp_handles.ventmaskbinary;
%                 ventilation_mask (cleaned_mask==1 & temp_handles.ventmaskbinary==1)= 1;
%
% else
%     ventilation_mask = vent_stack.*lungsmask;
% ventilation_mask (cleaned_mask==1)= 1;
%
% end
%
% ventilation_mask (cluster_stack==1 & truncated_vessel_mask==0 & cleaned_mask==0)= 2;
% %% lobe defect separation
% temp_handles=identify_lobes(temp_handles);
%
cleaned_mask = cleaned_mask.*double(temp_handles.lobemask);