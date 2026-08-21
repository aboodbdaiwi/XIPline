%% CartesianRecon
publish_options.codeToEvaluate=[ ...
'CartesianRecon(''../../../../recongold1/2D/256x128-1acq-4slice-1ch/P02048.7'');' char(10) ...
];
publish('CartesianRecon.m',publish_options);
 
clear mex;
clear all;
close all;

%% AssetRecon
publish_options.codeToEvaluate=[ ...
'AssetRecon(''../../../../recongold1/Asset/3D-256x128-1acq-28slice-8ch-Asset2x1/P94208.7'', ''../../../../recongold1/Calibration2D/32x32-1acq-25slice-8ch-Brain/P89088.7'');' char(10) ...
];
publish('AssetRecon.m',publish_options);
 
clear mex;
clear all;
close all;

%% EpiComputeAndApplyPhaseCorrection
publish_options.codeToEvaluate=[ ...
'EpiComputeAndApplyPhaseCorrection(''../../../../recongold1/EpiReferenceScan/MultiPhase-ASSET-128x128/P02048.7'');' char(10) ...
];
publish('EpiComputeAndApplyPhaseCorrection.m',publish_options);
 
clear mex;
clear all;
close all;

%% EpiRecon
publish_options.codeToEvaluate=[ ...
'EpiRecon(''../../../../recongold1/Epi/GE-EPI-TEMinFull-128x128-1Shot-1Slice-1ch-FastRampSamp/P23552.7'');' char(10) ...
];
publish('EpiRecon.m',publish_options);

clear mex;
clear all;
close all;

%% EpiMultiPhaseRecon
publish_options.codeToEvaluate=[ ...
'EpiMultiPhaseRecon(''../../../../recongold1/EpiMultiPhase/fMRI-Asset/P41472.7'');' char(10) ...
];
publish('EpiMultiPhaseRecon.m',publish_options);
 
clear mex;
clear all;
close all;

%% EpiDiffusionRecon
publish_options.codeToEvaluate=[ ...
'EpiDiffusionRecon(''../../../../recongold1/EpiDiffusion/1T2-1BVal-1DifDir-1NEX-Zerofill/P10752.7'');' char(10) ...
];
publish('EpiDiffusionRecon.m',publish_options);

clear mex;
clear all;
close all;

%% SpectroMultiChannelMultiVoxel
publish_options.codeToEvaluate=[ ...
'SpectroMultiChannelMultiVoxel(''../../../../recongold1/SpectroMCSI/8x8x1-8HRBrain-WithLocalizerAndCalData/P08704.7'', ''../../../../recongold1/SpectroMCSI/8x8x1-8HRBrain-WithLocalizerAndCalData/McsiCalibration.h5'');' char(10) ...
];
publish('SpectroMultiChannelMultiVoxel.m',publish_options);

clear mex;
clear all;
close all;

%% SpectroPerChannelMultiVoxel
publish_options.codeToEvaluate=[ ...
'SpectroPerChannelMultiVoxel(''../../../../recongold1/SpectroMultiVoxel/8x8x1-SDK-DemoData/P11264.7'');' char(10) ...
];
publish('SpectroPerChannelMultiVoxel.m',publish_options);

clear mex;
clear all;
close all;

%% SpectroSingleVoxelRecon
publish_options.codeToEvaluate=[ ...
'SpectroSingleVoxelRecon(''../../../../recongold1/SpectroSingleVoxel/ProbeP-8ch-SDKDemo/P12288.7'');' char(10) ...
];
publish('SpectroSingleVoxelRecon.m',publish_options);

clear mex;
clear all;
close all;

%% FlexRecon
publish_options.codeToEvaluate=[ ...
'FlexRecon(''../../../../recongold1/Flex/3D-Sag-256x224x52-8ch-BiPolar-ARC-FiltB-WatFatOp/P11776.7'', ''../../../../recongold1/Flex/3D-Sag-256x224x52-8ch-BiPolar-ARC-FiltB-WatFatOp/kacq_yz.txt.906151242'');' char(10) ...
];
publish('FlexRecon.m',publish_options);
 
clear mex;
clear all;
close all;

%% EpiScanArchiveRecon
publish_options.codeToEvaluate=[ ...
'EpiScanArchiveRecon(''../../../../recongold1/UnitTestData/SDK/Examples/Data/EpiScanArchive/ScanArchive.h5'');' char(10) ...
];
publish('EpiScanArchiveRecon.m',publish_options);
 
clear mex;
clear all;
close all;

%% Spiral2DRecon
publish_options.codeToEvaluate=[ ...
'Spiral2DRecon(''../../../../recongold1/UnitTestData/SDK/Examples/Data/Spiral2D/P11264.7'');' char(10) ...
];
publish('Spiral2DRecon.m',publish_options);
 
clear mex;
clear all;
close all;

%% ASL3dScanArchiveRecon
publish_options.codeToEvaluate=[ ...
'ASL3dScanArchiveRecon(''../../../../recongold1/UnitTestData/SDK/Examples/Data/3DASL/ScanArchive.h5'');' char(10) ...
];
publish('ASL3dScanArchiveRecon.m',publish_options);
 
clear mex;
clear all;
close all;