function [transformed_data,standardized_hist,bin_counts]=hist_transformation(V,mu_s)
%
% step 2: map to the new intensity scale [1,512] with 75% max intensity 
% landmark maped to the 25% max intensity for a contrast enhancement on
% vessels inside the lungs.
% 
% Ref:Udupa JK. On Standardizing the MR Image Intensity Scale. Magn Reson 
% Med. 1999;42(6):1072-1081. 
% W. Zha @2015

[nRows,nCols,nSlices]=size(V);
all_pixels = V(:);
max_intn=max(all_pixels);
[bin_counts,intn_bins]=hist(V(:),max(max_intn,50));

mu_i = .25*max_intn;


pixel_list1=all_pixels(all_pixels<=mu_i);
pixel_list2=all_pixels(all_pixels>mu_i);

pc1=0.001;
pc2=0.998;
s1=1;
s2=512;
p1= pc1*max_intn;
p2= pc2*max_intn;


standarized_list1 = round(mu_s+(pixel_list1-mu_i)*(s1-mu_s)/(p1-mu_i));
standarized_list2 = round(mu_s+(pixel_list2-mu_i)*(s2-mu_s)/(p2-mu_i));

standarized_values = zeros(length(all_pixels),1);
standarized_values(all_pixels<=mu_i)=standarized_list1;
standarized_values(all_pixels>mu_i)=standarized_list2;
transformed_data=reshape(standarized_values,nRows,nCols,nSlices);


standardized_hist(intn_bins<=mu_i)= mu_s+(intn_bins(intn_bins<=mu_i)-mu_i)*(s1-mu_s)/(p1-mu_i);
standardized_hist(intn_bins>mu_i)= mu_s+(intn_bins(intn_bins>mu_i)-mu_i)*(s1-mu_s)/(p1-mu_i);