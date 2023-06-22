
function [Proton,Ventilation,Diffusion,GasExchange] = LoadMask(Proton,Ventilation,Diffusion,GasExchange,MainInput)
        
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