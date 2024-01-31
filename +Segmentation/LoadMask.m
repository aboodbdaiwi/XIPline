function [Proton,Ventilation,Diffusion,GasExchange] = LoadMask(Proton,Ventilation,Diffusion,GasExchange,MainInput)
%   Inputs:
%      
%   Outputs:
%                   
%   Package: https://github.com/aboodbdaiwi/HP129Xe_Analysis_App
%
%   Author: Abdullah S. Bdaiwi
%   Work email: abdullah.bdaiwi@cchmc.org
%   Personal email: abdaiwi89@gmail.com
%   Website: https://www.cincinnatichildrens.org/research/divisions/c/cpir
%
%   Please add updates at the end. Ex: 3/10/24 - ASB: update .... 

        [filename, path] = uigetfile({'*.mat';'*.dcm'},'Select File');
        [~,~,ext] = fileparts(filename);
        if  strcmp(ext, '.mat') 
            LungMask = load([path,filename]);
        elseif strcmp(ext, '.dcm') 
            LungMask = double(dicomread(filename));
        end
        
        % store lung mask
        switch MainInput.AnalysisType
            case 'Ventilation'
                Ventilation.LungMask = LungMask;
            case 'Diffusion'
                Diffusion.LungMask = LungMask;
            case 'GasExchange'
                GasExchange.LungMask = LungMask;
        end  

end