function outStruct = PK_3T_1C13_pyr()
%% .M file to assemble the bounds, priorKnowledge and initialValues structs for the matlab implementation of AMARES
%
%Each of B, PK and IV is a 1xN struct, where N is the number of peaks. Note
%multiplets are counted as one peak.
%The fields are as follows:
%bounds           initialValues          priorKnowledge
%
%peakName         peakName               peakName
%chemShift        chemShift              multiplet
%linewidth        linewidth              chemShiftDelta
%amplitude        amplitude              amplitudeRatio
%phase            phase                  G_linewidth
%chemShiftDelta                          G_amplitude
%amplitudeRatio                          G_phase
%                                        G_chemShiftDelta
%                                        refPeak
%
%% 3T file
% V1 Mary McLean
% V2 Alixander Khan
% V3 Rolf Schulte


%% Bounds
fields.Bounds = {
    'peakName', 'chemShift', 'linewidth', 'amplitude', 'phase', 'chemShiftDelta', 'amplitudeRatio'};
values.boundsCellArray = {...
    'Pyruvate', [-0.2,0.2],  [2,15],      [0,inf],     [-180,180], [],            [];
    'Lactate',  [11.5,12.8], [2,15],      [0,inf],     [-180,180], [],            [];
    'Bicarb',   [-10.3,-9.9],[2,15],      [0,inf],     [-180,180], [],            [];
    'Alanine',  [5.3,5.9],   [2,15],      [0,inf],     [-180,180], [],            [];
    'Pyr-H'     [7,9],       [2,20],      [0,inf],     [-180,180], [],            [];
    };


%% initialValues
fields.IV = {
    'peakName', 'chemShift', 'linewidth', 'amplitude', 'phase'};
values.IVCellArray = {...
    'Pyruvate',  0,           10,          10,          0;
    'Lactate',  12.2,         10,          5,           0;
    'Bicarb',   -10.1,        10,          2,           0;
    'Alanine',  5.5,          10,          1,           0;
    'Pyr-H',    8.3,          10,          1,           0;
    };


%%
fields.PK = {
    'peakName', 'multiplet', 'chemShiftDelta', 'amplitudeRatio', 'G_linewidth', 'G_amplitude', 'G_phase' ,'RelPhase', 'G_chemShiftDelta', 'refPeak'};
values.PKCellArray = {...
    'Pyruvate', [],          [],               [],               1,             [],            1,        [],          [],                 1;
    'Lactate',  [],          [],               [],               1,             [],            1,        [],          [],                 [];
    'Bicarb',   [],          [],               [],               1,             [],            1,        [],          [],                 [];
    'Alanine',  [],          [],               [],               [],            [],            1,        [],          [],                 [];
    'Pyr-H',    [],          [],               [],               1,             [],            1,        [],          [],                 [];
    };


%% Pass to the function which assembles the constraints into structs and saves them
outStruct = AMARES.priorKnowledge.preparePriorKnowledge(fields,values);
outStruct.svnVersion = '$Rev: 3 $';
outStruct.svnHeader = '$Header: RFS $';
