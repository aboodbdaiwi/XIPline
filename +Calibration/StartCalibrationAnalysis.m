
function [GasExResults, CalResults] = StartCalibrationAnalysis(MainInput)
%   Inputs:
%      
%
%   Outputs:
%                   
%
%   Package: https://github.com/aboodbdaiwi/HP129Xe_Analysis_App
%
%
%   Author: Abdullah S. Bdaiwi
%   Work email: abdullah.bdaiwi@cchmc.org
%   Personal email: abdaiwi89@gmail.com
%   Website: https://www.cincinnatichildrens.org/research/divisions/c/cpir

% Initialize output variables
%Results to Pass to Gas Exchange Recon
GasExResults = [];
%Results to Pass to Calibration information
CalResults= [];

    % call functions 
    if strcmp(MainInput.XeDataext,'.data') == 1 && (strcmp(MainInput.Institute, 'XeCTC')==1 ||...
            strcmp(MainInput.Institute, 'CCHMC') == 1 )                 
        [GasExResults, CalResults] = Calibration.XeCTC_Calibration(MainInput); 
    elseif strcmp(MainInput.XeDataext,'.dat') == 1 && (strcmp(MainInput.Institute, 'Duke') == 1 ||...
            strcmp(MainInput.Institute, 'UVA') == 1)
        [GasExResults, CalResults] = Calibration.Xe_duke_UVA_calibration(MainInput); 
        
    %-------------------- add new function here -------------------------------
    % elseif strcmp(DataType,'add DataType') == 1
        %  [GasExResults, CalResults] = add your function (DataLocation); 
    end 
end