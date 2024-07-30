
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

if strcmp(MainInput.AnalysisType,'Ventilation')                 
    [Proton,MainInput,GasExchange] = Registration.GeneralRegisterProton_to_Xenon(...
    Proton,...
    Ventilation,...
    MainInput,GasExchange);     

elseif strcmp(MainInput.AnalysisType,'GasExchange') 
        [Proton,MainInput,GasExchange] = Registration.GeneralRegisterProton_to_Xenon(...
        Proton,...
        GasExchange.VentImage,...
        MainInput, GasExchange); 
    % if (strcmp(MainInput.Institute,'CCHMC') == 1 || strcmp(MainInput.Institute,'XeCTC')) &&...
    %         strcmp(MainInput.SequenceType, '3D Radial') && MainInput.NoProtonImage == 0 && (strcmp(MainInput.Scanner, 'Philips') ||...
    %         strcmp(MainInput.Scanner,'Siemens') ) 
    %     [Proton,GasExchange] = Registration.GasExchange_RegisterProton_to_Xenon(Proton,GasExchange,MainInput);
    % else
    %     [Proton] = Registration.GeneralRegisterProton_to_Xenon(...
    %     Proton,...
    %     GasExchange.VentImage,...
    %     MainInput);  
    % end 
end 

end
