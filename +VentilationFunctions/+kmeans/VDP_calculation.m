function VDP_calculation(filename_LungMask, filename_He, subj_ID, wholelung_VDPfile, regional_VDPfile)

    % load the deformed lung/lobe mask. 
    LungMask = metaImageRead(filename_LungMask);

    % k-means clustering using our previous code.
    clusters_3D = clustering_3D(filename_He);

    % get the number of He image clusters. Usually, one integral corresponds to one cluster.
    clusters_3D_labels = unique(clusters_3D);
    clusters_3D_label_length = length (clusters_3D_labels);
    
    % get the number of lobes. Usually, one lobel corresponds to one lobe.
    lobemask_labels = unique(LungMask);
    lobemask_labels = lobemask_labels(lobemask_labels>0);
    lobemask_label_length = length (lobemask_labels);
    
    % create a 2D array to save the regional VDP.
    % row represents lobes in the order of left upper---left lower---right upper---right
    % middle---right lower 
    lobe_wise_VDP = zeros(lobemask_label_length,clusters_3D_label_length + 1);
    
    % interate through each lobe% open and close file by file ID.
    for i = 1 : lobemask_label_length
        
        % find the current lobe i.
        lobe = (LungMask==lobemask_labels(i));
        lobe_wise_VDP(i,1) = lobemask_labels(i);
        % interate through each He cluster
        for j = 1 : clusters_3D_label_length
            % find the current cluster.
            cluster = (clusters_3D==clusters_3D_labels(j));
            % find cluster j in the current lobe i. for example, we want to
            % see cluster 0 in lobe 1.=, similarly, we can see cluster 1 in lobe 1
            cluster_i_in_lobe_j = lobe.*cluster;
            % calculate the region of the current cluster, in the units of
            % voxel.
            lobe_wise_VDP(i,j+1) = sum(cluster_i_in_lobe_j(:));
        end
        % normalize the region of current cluster to the current lobe for
        % regional VDP calculation.
        lobe_wise_VDP(i,2:clusters_3D_label_length + 1) = lobe_wise_VDP(i,2:clusters_3D_label_length + 1)/sum(lobe(:))*100;
    end
        
    % whole lung VDP calculation, we actually calculate the ratio of each
    % He image cluster over the whole lung.
    wholelung_VDP = zeros (clusters_3D_label_length,1);
    % get the whole lung volume.
    total_lung_volume = (LungMask~=0);
    CTlung_volume = sum (total_lung_volume(:));
    % whole lung volume size (derived from CT) to compare with H MRI manaul
    % lung volume size. 
    CTlung_volume = CTlung_volume*1.5625^3/1000000;
    for i = 1 : clusters_3D_label_length 
            cluster_within_lung = (clusters_3D==clusters_3D_labels(i)).*total_lung_volume;
            wholelung_VDP(i) = sum(cluster_within_lung(:))/sum(total_lung_volume(:))*100;
    end

    % write the regional VDP results to the text file with a name specified
    % by "regional_VDPfile"
    % open and close file by file ID.
    fid = fopen(regional_VDPfile,'a');
    fprintf(fid,'%s\t', subj_ID);

    for i = 1 : size(lobe_wise_VDP,1)
        for j = 1: size(lobe_wise_VDP,2)
            fprintf(fid,'%f\t',lobe_wise_VDP(i,j));
        end
        fprintf(fid,'\n');
        fprintf(fid,'\t\t');
    end
    fprintf(fid,'\b\r');
    fclose(fid);

     
    % write the regional VDP results as well as the CT-derived lung size to  
    % the text file with a name specified by "wholelung_VDPfile"
    % open and close file by file ID.
    fid = fopen(wholelung_VDPfile,'a');
    fprintf(fid,'%s\t\t', subj_ID);
    
    for i = 1: size(wholelung_VDP,1)
            fprintf(fid,'%f\t',wholelung_VDP(i));
    end
    fprintf(fid,'%s\t','Lung volume (L):');
    fprintf(fid,'%f\t',CTlung_volume);
    fprintf(fid,'\n');
    fclose(fid);
end
