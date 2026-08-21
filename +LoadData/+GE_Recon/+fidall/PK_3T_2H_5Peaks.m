function outStruct = PK_3T_2H_5Peaks(ppm_offset)
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
%% 3T Deuterium PK file
% V1. Created by Liam Young. University of Oxford.
% V2: Rolf Schulte, GEHC
% V3: Michael Vaeggemose, GEHC -- remove sgn as input
% V4: Alixander Khan & Michael Vaeggemose adding lipids and adjusting lactate
% V5: Alixander Khan & Rolf Schulte; tuning lipds/lactate; optional ppm_offset

if ~exist('ppm_offset','var'), ppm_offset = []; end


%% Bounds
fields.Bounds = {
    'peakName', 'chemShift', 'linewidth', 'amplitude', 'phase', 'chemShiftDelta', 'amplitudeRatio'};
values.boundsCellArray = {...
    'D2O',      [4.5,4.9],   [0,inf],     [0,inf],     [0,360], [],               [];
    'Glucose',  [3.6,4.0],   [9,11],      [0,inf],     [0,360], [],               [];
    'Glx',      [2.3,2.5],   [0,100],     [0,inf],     [0,360], [],               [];
    'Lactate',  [1.15,1.45], [0,100],     [0,inf],     [0,360], [],               [];
    'Lipid',    [0.7,1.0],   [0,100],     [0,inf],     [0,360], [],               [];
    };


%% initialValues
fields.IV = {
    'peakName', 'chemShift', 'linewidth', 'amplitude', 'phase'};
values.IVCellArray = {...
    'D2O',       4.7,         8,           10,          0;
    'Glucose',   3.8,         10,          3,           0;
    'Glx',       2.3,         8,           1,           0;
    'Lactate',   1.35,        8,           0.5,         0;
    'Lipid',     0.9,         8,           0.5,         0;
    };


%%
fields.PK = {
    'peakName', 'multiplet', 'chemShiftDelta', 'amplitudeRatio', 'G_linewidth', 'G_amplitude', 'G_phase', 'RelPhase', 'G_chemShiftDelta', 'refPeak'};
values.PKCellArray = {...
    'D2O',       [],          [],               [],               [],            [],            1,         [],         [],                 1;
    'Glucose',   [],          [],               [],               [],            [],            1,         [],         [],                 0;
    'Glx',       [],          [],               [],               [],            [],            1,         [],         [],                 0;
    'Lactate',   [],          [],               [],               [],            [],            1,         [],         [],                 0;
    'Lipid',     [],          [],               [],               [],            [],            [],        [],         [],                 0;
    };


%% Correct for freq offset in spectro-prescan (default: ppm_offset = 0)
if ~isempty(ppm_offset)
    for l=1:size(values.boundsCellArray,1)
        values.boundsCellArray{l,2} = values.boundsCellArray{l,2} + ppm_offset;
        values.IVCellArray{l,2} = values.IVCellArray{l,2} + ppm_offset;
    end
end


%% Pass to the function which assembles the constraints into structs and saves them
outStruct = AMARES.priorKnowledge.preparePriorKnowledge(fields,values);
outStruct.svnVersion = '$Rev: 5 $';
outStruct.svnHeader = '$Header: RFS $';
