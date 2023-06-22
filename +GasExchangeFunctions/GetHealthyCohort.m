function [HealthyCohort, CohortTag] = GetHealthyCohort()
%GetHealthyCohort - Returns a healthy cohort structure for use in GasExchange processing
%   Either uncomment the cohort desired or create a new one
%   Add any new healthy subject to all applicable cohorts
%   If changing/adding a cohort update its Name/Date as well

%% Current Healthy Subjects Acquired
%%Any Notes following data and subject for folder name will not be used
%%even if in the cohort
%    'IRC186H-281',...%7F
%    'IRC186H-274',...%9M - Broken Coil; unused
%    'IRC186H-273',...%11M - Broken Coil; unused
%    'IRC186H-327',...%11M
%    'IRC186H-283',...%23F
%    'IRC186H-326',...%20F
%    'IRC186H-328',...%40M
%    'IRC186H-332',...%8M - noticeable ventilation defect
%    'IRC186H-336',...%28F
%    'IRC186H-338',...%5M - poorish breath hold
%    'IRC186H-342',...%39F
%    'IRC186H-348',...%29F
%    'IRC186H-349',...%20F
%    'IRC186H-1022',...%28M
%    'IRC186H-1022',...%28M - 3TR Dual Loop; unused since only want 1 of each subject
%    'IRC186H-1022',...%28M - 3TR Polarean; unused since only want 1 of each subject
%    'IRC186H-1019',...%28M - noticeable ventilation defect
%    'IRC186H-1026',...%30M
%    'IRC186H-1010',...%36M
%    'IRC186H-1010',...%36M - exercise; unused since only want 1 of each subject



%% Potential Cohorts

%Full Cohort
HealthyCohort = {
    'IRC186H-281',...%7F
    'IRC186H-274',...%9M - Broken Coil; unused
    'IRC186H-273',...%11M - Broken Coil; unused
    'IRC186H-327',...%11M
    'IRC186H-283',...%23F
    'IRC186H-326',...%20F
    'IRC186H-328',...%40M
    'IRC186H-332',...%8M - noticeable ventilation defect
    'IRC186H-336',...%28F
    'IRC186H-338',...%5M - poorish breath hold
    'IRC186H-342',...%39F
    'IRC186H-348',...%29F
    'IRC186H-349',...%18F
    'IRC186H-1022',...%28M - 3 acqs with only the 3T-T1 Polarean acq used
    'IRC186H-1019',...%28M
    'IRC186H-1026',...%30M
    'IRC186H-1010',...%36M - 2 acqs with only the baseline (non-exercise) used
    }';
CohortTag = '01142021_FullHealthyCohort';

% %Pediatric Cohort <=18
% HealthyCohort = {
%     'IRC186H-338',...%5M - poorish breath hold
%     'IRC186H-281',...%7F
%     'IRC186H-332',...%8M - noticeable ventilation defect
%     'IRC186H-274',...%9M - Broken Coil; unused
%     'IRC186H-273',...%11M - Broken Coil; unused
%     'IRC186H-327',...%11M
%     'IRC186H-349',...%18F
%     }';
% CohortTag = '11182020_PediatricHealthyCohort';

% %Adult Cohort >18
% HealthyCohort = {
%     'IRC186H-326',...%20F
%     'IRC186H-283',...%23F
%     'IRC186H-336',...%28F
%     'IRC186H-1022',...%28M - 3 acqs with only the 3T-T1 Polarean acq used
%     'IRC186H-1019',...%28M
%     'IRC186H-348',...%29F
%     'IRC186H-1026',...%30M
%     'IRC186H-1010',...%36M - 2 acqs with only the baseline (non-exercise) used
%     'IRC186H-342',...%39F
%     'IRC186H-328',...%40M
%     }';
% CohortTag = '01142021_AdultHealthyCohort';

end