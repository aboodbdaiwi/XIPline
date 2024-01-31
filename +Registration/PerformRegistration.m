
function [Proton,GasExchange,MainInput] = PerformRegistration(Proton,Ventilation,GasExchange,MainInput)
%   Inputs:
%          
%   Outputs:
%      LungMask
%
%   Example: 
%   Package: 
%
%   Author: Abdullah Bdaiwi 
%   Work email: abdullah.bdaiwi@cchmc.org
%   Personal email: abdaiwi89@gmail.com
%   Website: https://www.cincinnatichildrens.org/research/divisions/c/cpir

%% Initialize output variables

if strcmp(MainInput.AnalysisType,'Ventilation') == 1                 
    [Proton,MainInput] = Registration.GeneralRegisterProton_to_Xenon(...
    Proton,...
    Ventilation,...
    MainInput);     

elseif strcmp(MainInput.AnalysisType,'GasExchange') == 1 
    if (strcmp(MainInput.Institute,'CCHMC') == 1 || strcmp(MainInput.Institute,'XeCTC') == 1) && strcmp(MainInput.SequenceType, '3D Radial') == 1 && MainInput.NoProtonImage == 0  
        [Proton,GasExchange] = Registration.GasExchange_RegisterProton_to_Xenon(Proton,GasExchange,MainInput);
    else
        [Proton] = Registration.GeneralRegisterProton_to_Xenon(...
        Proton.Image,...
        GasExchange.VentImage,...
        MainInput);  
    end 
end 

end
