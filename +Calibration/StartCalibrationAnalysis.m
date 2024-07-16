
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
    if strcmp(MainInput.CalDataext,'.data') && (strcmp(MainInput.Institute, 'XeCTC') ||...
            strcmp(MainInput.Scanner, 'Philips'))                 
        [GasExResults, CalResults] = Calibration.XeCTC_Calibration(MainInput); 
    elseif strcmp(MainInput.CalDataext,'.dat') && (strcmp(MainInput.Institute, 'XeCTC') ||...
             strcmp(MainInput.Scanner, 'Siemens'))
        [GasExResults, CalResults] = Calibration.Xe_duke_UVA_calibration(MainInput); 
    elseif strcmp(MainInput.CalDataext,'.7') && strcmp(MainInput.Scanner, 'GE')   
            LoadData.ismrmrd.GE.calibration_to_ismrmrd([MainInput.Cal_name '.h5'],MainInput);
            MainInput.XeFileName = [MainInput.Cal_name '.h5'];
        [GasExResults, CalResults] = Calibration.XeCTC_Calibration_GEMRD(MainInput);     
    elseif (strcmp(MainInput.CalDataext,'.h5') || strcmp(MainInput.CalDataext,'.MRD') ||...
            strcmp(MainInput.CalDataext,'.mrd')) && strcmp(MainInput.Scanner, 'GE')   
        [GasExResults, CalResults] = Calibration.XeCTC_Calibration_GEMRD(MainInput);         
    elseif strcmp(MainInput.CalDataext,'.h5') || strcmp(MainInput.CalDataext,'.MRD') || strcmp(MainInput.CalDataext,'.mrd') 
        [GasExResults, CalResults] = Calibration.XeCTC_Calibration_MRD(MainInput);         

        %-------------------- add new function here -------------------------------
    % elseif strcmp(DataType,'add DataType') == 1
        %  [GasExResults, CalResults] = add your function (DataLocation); 
    end 
end