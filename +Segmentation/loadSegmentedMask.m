
function [Proton,Ventilation,Diffusion, GasExchange] = loadSegmentedMask(Proton,Ventilation,Diffusion,GasExchange,MainInput)
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

cd(MainInput.XeDataLocation);
mask_existing1 = exist('mask','var');
mask_existing2 = exist('Mask','var');
mask_existing3 = exist('MASK','var');

% load mask
if mask_existing1 == 1 || mask_existing2 == 1 || mask_existing3 == 1
    mask_existing = 1;
    if mask_existing == 1
        Mask = double(mask);
    elseif mask_existing2 == 1
        Mask = double(Mask);
    elseif mask_existing3 == 1
        Mask = double(MASK);
    end
elseif isfile('mask.mat')                                
        Mask = load([MainInput.XeDataLocation,'\mask.mat']);
        Mask = double(Mask.labels);
        mask_existing = 1;
elseif isfile('Mask.mat')                                
        Mask = load([MainInput.XeDataLocation,'\Mask.mat']);
        Mask = double(Mask.labels);
        mask_existing = 1;
elseif isfile('MASK.mat')                                
        Mask = load([MainInput.XeDataLocation,'\MASK.mat']);
        Mask = double(Mask.labels);
        mask_existing = 1;
else
    mask_existing = 0;
end
        
% store mask
if mask_existing == 1
    % assign lung and airway mask
    Mask(isinf(Mask)|isnan(Mask)) = 0;
    lungmask = zeros(size(Mask));
    lungmask(Mask == 1)= 1;
    airwaymask = zeros(size(Mask));
    airwaymask(Mask == 2)= 1;  
    
    % store mask
    switch MainInput.AnalysisType
        case 'Ventilation'
           Ventilation.Mask = Mask;
           Ventilation.LungMask = lungmask;
           Ventilation.AirwayMask = airwaymask;                       
        case 'Diffusion'
           Diffusion.Mask = Mask;
           Diffusion.LungMask = lungmask;
           Diffusion.AirwayMask = airwaymask;
        case 'GasExchange'
           GasExchange.Mask = Mask;
           GasExchange.LungMask = lungmask;
           Proton.LungMask = lungmask;
           Proton.ProtonMaskRegistred = lungmask;
           GasExchange.AirwayMask = airwaymask;                           
    end    
    
    % Save images and maskes
    mkdir([MainInput.XeDataLocation, '\Mask']);
    % Mask 
    niftiwrite(abs(fliplr(rot90(Mask,-1))),[MainInput.XeDataLocation,'\Mask\Mask.nii'],'Compressed',true);
    info = niftiinfo([MainInput.XeDataLocation,'\Mask\Mask.nii.gz']);
    info.Description = strcat('Package Version: ', 'Version1');
    niftiwrite(abs(fliplr(rot90(Mask,-1))),[MainInput.XeDataLocation,'\Mask.nii'],info,'Compressed',true);
    % lung mask
    niftiwrite(abs(fliplr(rot90(lungmask,-1))),[MainInput.XeDataLocation,'\Mask\LungMask.nii'],'Compressed',true);
    info = niftiinfo([MainInput.XeDataLocation,'\Mask\LungMask.nii.gz']);
    info.Description = strcat('Package Version: ', 'Version1');
    niftiwrite(abs(fliplr(rot90(lungmask,-1))),[MainInput.XeDataLocation,'\Mask\LungMask.nii'],info,'Compressed',true);
    % airway mask
    niftiwrite(abs(fliplr(rot90(airwaymask,-1))),[MainInput.XeDataLocation,'\Mask\AirwayMask.nii'],'Compressed',true);
    info = niftiinfo([MainInput.XeDataLocation,'\Mask\AirwayMask.nii.gz']);
    info.Description = strcat('Package Version: ', 'Version1');
    niftiwrite(abs(fliplr(rot90(airwaymask,-1))),[MainInput.XeDataLocation,'\Mask\AirwayMask.nii'],info,'Compressed',true);
                    
    % update app progress message statement           
    MainInput.Status = 'Lung Mask Segmentation completed....';              
else
    % update app progress message statement           
    MainInput.Status = 'Lung Mask does not exist...';                                          
end

end

