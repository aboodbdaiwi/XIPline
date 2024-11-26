function [mu_i,mu_s,mapped_intn_bins,bin_counts] = estimate_mu_from_training_data(V,params)
%
% step 1: Estimate the intensity value corresponding to the 25% max 
% intensity for a contrast enhancement inside the lungs.
% 
% Ref:Udupa JK. On Standardizing the MR Image Intensity Scale. Magn Reson 
% Med. 1999;42(6):1072-1081. 
%
% Copyright W. Zha @2015
pc1=0;
pc2=0.998;
s1=1;
s2=512;
max_intn=max(V(:));
[bin_counts,intn_bins]=hist(V(:),100);
p1= pc1*max_intn;
p2= pc2*max_intn;


mu_i = .75*max_intn;

mapped_intn_bins=s1+(intn_bins-p1)*(s2-s1)/(p2-p1);
mu_s=s1+(mu_i-p1)*(s2-s1)/(p2-p1);
