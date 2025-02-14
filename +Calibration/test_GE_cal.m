

d = 'D:\OneDrive - cchmc\Lab\Xe_App\testtingData\Test_Data_forXIPline\Gasexchange\GE_ROB0046-029\Testing\Cal\P17920.7';
h = [];
wfnpath = 'C:\XIPline\GE\waveforms\xe_calibration'; % files inside the wfnpath folder = {'*.mat', '*freq.fdl', '*flip.fdl'};
[meanRbc2barrier,te90,targetAX,targetTG] = Calibration.recon_calibration(d,h,wfnpath);

%% 
MainInput.CalFullPath = 'D:\OneDrive - cchmc\Lab\Xe_App\testtingData\Test_Data_forXIPline\Gasexchange\GE_ROB0046-029\Testing\Cal\P17920.7';
[GasExResults, CalResults] = Calibration.GE_calibration(MainInput);

%% 

d = 'D:\OneDrive - cchmc\Lab\Xe_App\testtingData\Test_Data_forXIPline\Gasexchange\GE_ROB0046-029\Testing\DP\raw\P18944.7';
h = [];
wfnpath = 'C:\XIPline\GE\waveforms\xe_dissolved';
wfn_freq = '';
fname = d;
LoadData.ismrmrd.GE.prep_dixon_gx(d,h,wfnpath,wfn_freq,fname)
