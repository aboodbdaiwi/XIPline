
function [Proton, Ventilation, GasExchange, MainInput] = PerformRegistration(Proton,Ventilation,GasExchange,MainInput)
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
    switch MainInput.RegistrationType
        case 'ANTs'
            [MainInput, Proton, Ventilation, GasExchange] = Registration.ANTs.AntsRegistration(MainInput, Proton, Ventilation, GasExchange);
        case 'Multimodal'
            [Proton, Ventilation, GasExchange, MainInput] = Registration.GeneralRegisterProton_to_Xenon(...
            Proton,...
            Ventilation,...
            MainInput,GasExchange);   
        case 'Manual'
            [MainInput, Proton, Ventilation, GasExchange] = Registration.ControlPointsRegistration(MainInput, Proton, Ventilation, GasExchange);
        case 'Pretrained DL Model'
            % will be added later 
    end

elseif strcmp(MainInput.AnalysisType,'GasExchange') 
        [Proton,GasExchange] = Registration.GasExchange_RegisterProton_to_Xenon(Proton,GasExchange,MainInput);

    % switch MainInput.RegistrationType
    %     case 'ANTs'
    %         %[MainInput, Proton,  Ventilation, GasExchange] = Registration.ANTs.AntsRegistration(MainInput, Proton, Ventilation, GasExchange);
    %     case 'Multimodal'
    %         [Proton, Ventilation, GasExchange, MainInput] = Registration.GeneralRegisterProton_to_Xenon(...
    %         Proton,...
    %         GasExchange.VentImage,...
    %         MainInput, GasExchange); 
    %     case 'Manual'
    % 
    %     case 'Pretrained DL Model'
    % end 
    % if (strcmp(MainInput.Institute,'CCHMC') == 1 || strcmp(MainInput.Institute,'XeCTC')) &&...
    %         strcmp(MainInput.SequenceType, '3D Radial') && MainInput.NoProtonImage == 0 && (strcmp(MainInput.Scanner, 'Philips') ||...
    %         strcmp(MainInput.Scanner,'Siemens') ) 
    %     [Proton,GasExchange] = Registration.GasExchange_RegisterProton_to_Xenon(Proton,GasExchange,MainInput);
    % else
    %     [Proton, Ventilation, GasExchange, MainInput] = Registration.GeneralRegisterProton_to_Xenon(...
    %     Proton,...
    %     GasExchange.VentImage,...
    %     MainInput, GasExchange);  
    % end 
end 

end
