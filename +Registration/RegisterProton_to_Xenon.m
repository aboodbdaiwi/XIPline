
function [Proton,GasExchange,MainInput] = RegisterProton_to_Xenon(Proton,Ventilation,GasExchange,MainInput)
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
%   Website: https://cpir.cchmc.org/

%% Initialize output variables
% ImageRegistration = [];


%% 

if strcmp(MainInput.AnalysisType,'Ventilation') == 1                 
    [Proton,MainInput] = Registration.GeneralRegisterProton_to_Xenon(...
    Proton,...
    Ventilation,...
    MainInput);     

elseif strcmp(MainInput.AnalysisType,'GasExchange') == 1 
    if strcmp(MainInput.Institute,'CCHMC') == 0  && strcmp(MainInput.Institute,'XeCTC') == 0             
        [Proton] = Registration.GeneralRegisterProton_to_Xenon(...
        Proton.Image,...
        GasExchange.VentImage,...
        MainInput);     

    elseif (strcmp(MainInput.Institute,'CCHMC') == 1 || strcmp(MainInput.Institute,'XeCTC') == 1)&& strcmp(MainInput.SequenceType, '3D Radial') == 1  &&...
           strcmp(MainInput.XeDataext,'.data') == 1  && strcmp(MainInput.Scanner, 'Philips') == 1    
        [Proton,GasExchange] = Registration.GasExchange_RegisterProton_to_Xenon(Proton,GasExchange,MainInput);
    end 
end 

end
