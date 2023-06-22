
function [Proton,Ventilation,Diffusion, GasExchange] = PerformSegmention(Proton,Ventilation,Diffusion,GasExchange,MainInput)
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

%% 
% choose image
switch MainInput.AnalysisType
    case 'Ventilation'
        Image = Ventilation.Image;
    case 'Diffusion'
        Image = Diffusion.Image;
    case 'GasExchange'
        Image = GasExchange.VentImage;
end
% perform segmention
switch MainInput.Segment
    case 'parenchyma'
        switch MainInput.SegmentMethod
            case 'threshold'
                LungMask = Segmentation.SegmrntLungthresh(Image,MainInput.SE,MainInput.thresholdlevel);
            case 'manual'
                switch MainInput.SegmentManual
                    case 'basic'
                        LungMask = Segmentation.SegmentLungParenchyma(Image);
                    case 'manual'
                        volumeSegmenter(Image);
                        if isfile('lungmask.mat')
                            LungMask = load([MainInput.XeDataLocation,'\lungmask.mat']);
                            LungMask = double(LungMask);
                        else
                             disp('File does not exist.')
                        end
                end
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
    case 'airway'
        AirwayMask = Segmentation.SegmentAirway(Image);
        % store airway mask
        switch MainInput.AnalysisType
            case 'Ventilation'
                Ventilation.AirwayMask = AirwayMask;
            case 'Diffusion'
                Diffusion.AirwayMask = AirwayMask;
            case 'GasExchange'
                GasExchange.AirwayMask = AirwayMask;
        end              
end

  
end

