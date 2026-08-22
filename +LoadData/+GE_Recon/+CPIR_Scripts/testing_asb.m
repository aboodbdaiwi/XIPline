
% Calibration   | 01_cal_cpir_dissolved
clc; clear; close all; 
d = 'D:\OneDrive - cchmc\Lab\GE_Acquisition_And_Recon-main\CPIR_Data\01_cal_cpir_dissolved\ScanArchive_513636T1MR_20260728_142337356_Cartesian2D.h5';
h = [];
% h = 'D:\OneDrive - cchmc\Lab\GE_Acquisition_And_Recon-main\CPIR_Data\01_cal_cpir_dissolved\ScanArchive_513636T1MR_20260728_142337356_Cartesian2D.hdr';
sv = 'cal_cpir_dissolved';
CPIR_GERecon.recon_calibration_cchmc(d,h,sv)

%% Diffusion |   07_diff_cpir_XXCart4bs2d

clc; clear; close all; 
d = 'D:\OneDrive - cchmc\Lab\GE_Acquisition_And_Recon-main\CPIR_Data\07_diff_cpir_XXCart4bs2d\Series11_bag\ScanArchive_513636T1MR_20260318_123238925_Cartesian2D.h5';
h = 'D:\OneDrive - cchmc\Lab\GE_Acquisition_And_Recon-main\CPIR_Data\07_diff_cpir_XXCart4bs2d\Series11_bag\ScanArchive_513636T1MR_20260318_123238925_Cartesian2D.hdr';
wfn = 'D:\OneDrive - cchmc\Lab\GE_Acquisition_And_Recon-main\CPIR_Scripts\diff_cpir_2dCart4bs\diffcart2D_129Xe_fov400_mtx80_nexc320_kdt52_gmax29_smax98_dur14p9_co.mat';
[bb,bbabs] = CPIR_GERecon.recon_cart_diffusion(d,h,wfn);



%% GX   |  02_vent_cpir_floret
clc; clear; close all; 

d = 'D:\OneDrive - cchmc\Lab\GE_Acquisition_And_Recon-main\CPIR_Data\02_vent_cpir_floret\ScanArchive_513636T1MR_20260728_144503172_Cartesian2D.h5';
h = 'D:\OneDrive - cchmc\Lab\GE_Acquisition_And_Recon-main\CPIR_Data\02_vent_cpir_floret\ScanArchive_513636T1MR_20260728_144503172_Cartesian2D.hdr';
wfn = 'D:\OneDrive - cchmc\Lab\GE_Acquisition_And_Recon-main\CPIR_Scripts\vent_cpir_floret\floret_xe129_fov375_mtx97_arms7p1_hub3_0p4_intlv1312_bonus10_kdt8_gmax21_smax139_dur7p7_M21.mat';
mtx = [192,192,192];
[bb,bbabs] = LoadData.GE_Recon.CPIR_Scripts.CPIR_GERecon.recon_grid_interp(d,h,wfn,mtx);


figure; imslice(abs(bb));












