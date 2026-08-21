function outStruct = PK_3T_129Xe_3Peaks(ppm_offset)
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
%% 3T Xenon PK file
% version 1. Created by Michael Vaeggemose, GE HC
% version 2. Created by Michael Vaeggemose and Jack J. Miller - add G_sigma & G_phase

if ~exist('ppm_offset','var'), ppm_offset = []; end


%% Bounds
fields.Bounds = {
    'peakName', 'chemShift', 'linewidth', 'amplitude', 'phase', 'chemShiftDelta', 'amplitudeRatio', 'sigma'};
values.boundsCellArray = {...
    'Gas',      [-25,25],    [0,inf],     [0,inf],     [0,360], [],               [],               [0,inf];
    'Tissue',   [192,205],   [0,100],     [0,inf],     [0,360], [],               [],               [0,150];
    'Blood',    [205,215],   [0,100],     [0,inf],     [0,360], [],               [],               [0,150];
    };


%% initialValues
fields.IV = {
    'peakName', 'chemShift', 'linewidth', 'amplitude', 'phase', 'sigma'};
values.IVCellArray = {...
    'Gas',       0,           40,          1,           0,       100;
    'Tissue',    197,         10,          0.6,         0,       10;
    'Blood',     210,         10,          0.2,         0,       10;
    };


%%
fields.PK = {
    'peakName', 'multiplet', 'chemShiftDelta', 'amplitudeRatio', 'G_linewidth', 'G_amplitude', 'G_phase', 'RelPhase', 'G_chemShiftDelta', 'refPeak', 'G_sigma'};
values.PKCellArray = {...
    'Gas',       [],          [],               [],               [],            [],            1,         [],         [],                 1, [];
    'Tissue',    [],          [],               [],               [],            [],            2,         [],         [],                 0, [];
    'Blood',     [],          [],               [],               [],            [],            2,         [],         [],                 0, [];
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
outStruct.svnVersion = '$Rev: 8034 $';
outStruct.svnHeader = '$Header: RFS $';
