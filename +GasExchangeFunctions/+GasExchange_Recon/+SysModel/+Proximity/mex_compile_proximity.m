% Get the starting directory
rootDistribDir = pwd();

% Change to the mex directory
cd([rootDistribDir filesep() 'GasExchangeFunctions' filesep()...
    '+GasExchange_Recon' filesep() '+SysModel' filesep() '+Proximity' filesep()]);

% Compile mex code
mex -g sparse_gridding_distance_mex.c;

% Change back to starting directory 
cd(rootDistribDir);