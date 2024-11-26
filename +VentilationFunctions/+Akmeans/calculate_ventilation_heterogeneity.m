function [heterogeneity_index,CoVstd,CoVmap,localScale,gas_fraction] = calculate_ventilation_heterogeneity(gas,lungsmask,ventilated_mask,gas_volsize)

gas_lung = gas.*lungsmask;
[~,noise_mean] = estimate_image_noise(gas,lungsmask);
gas_lung = gas_lung - noise_mean;
gas_lung(gas_lung <0) = 0;

% gas_lung = gas_lung.*ventilated_mask;

localScale = get_local_scale(lungsmask);
% localScale =5;
C_intToVol = 1/sum(gas_lung(ventilated_mask>0));
gas_fraction = C_intToVol*gas_lung.*lungsmask/gas_volsize;

% gas_fraction = C_intToVol*gas_lung.*ventilated_mask/gas_volsize;
gas_fraction = gas_fraction./max(gas_fraction(:));
% gas_fraction(lungsmask==0) = -1;

mea_kernel = strel('square',localScale);
%%
contoured_slices =  find(squeeze(sum(sum(ventilated_mask,1),2)));
nSlices = length(contoured_slices);
[nRows, nCols, ~] = size(ventilated_mask);
% CoVvalues = [];
CoVmap = zeros(size(ventilated_mask));
for sl = 1:nSlices
    sliceVent = ventilated_mask(:,:,contoured_slices(sl));
    sliceGasFrac  = gas_fraction(:,:,contoured_slices(sl));

%         sliceGas = gas(:,:,contoured_slices(sl));
    pixel_list = find(sliceVent(:)==1);
    [row_list,col_list] = ind2sub([nRows, nCols],pixel_list);
    CoV_list = zeros(length(pixel_list),1);
    CoVSlicemap = zeros(nRows,nCols);
    for n = 1:length(pixel_list)
        this_row = row_list(n);
        this_col = col_list(n);
        pixel_neighbor = zeros(nRows, nCols);
        pixel_neighbor(this_row,this_col) = 1;

pixel_neighbor = imdilate(pixel_neighbor,mea_kernel);
        gasfrac_neighbor = pixel_neighbor.*sliceGasFrac;
%         neighbor_Vfrac  = gasfrac_neighbor(pixel_neighbor==1);

           neighbor_Vfrac  = gasfrac_neighbor(gasfrac_neighbor>0);
     
        CoV_list(n) = std(neighbor_Vfrac)/mean(neighbor_Vfrac);

    end
    CoVSlicemap(pixel_list) = CoV_list;
    CoVmap(:,:,contoured_slices(sl)) = CoVSlicemap;
%     CoVvalues = [CoVvalues;CoV_list]; %#ok<AGROW>
    
end
CoVvalues = CoVmap(ventilated_mask==1);
CoVvalues(isnan(CoVvalues))=[];
heterogeneity_index = mean(CoVvalues);
CoVstd = std(CoVvalues);