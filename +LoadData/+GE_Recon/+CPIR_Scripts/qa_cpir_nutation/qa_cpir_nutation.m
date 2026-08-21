% Nutation, CPIR QA
% In theory, can just use default spectroscopy fidall2 mode and use custom 
% rf flip angles

% Set params
init_flip=10;
final_flip=180;
n_flips=18;
flips = linspace(init_flip,final_flip,n_flips);

filename = sprintf('qa_cpir_nutation_nexc%d_start%d_end%d_flip.fdl',n_flips,init_flip,final_flip);
write_fdl(filename,flips,'flip');

