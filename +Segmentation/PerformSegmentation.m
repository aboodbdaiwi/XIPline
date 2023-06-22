
function [Proton,Ventilation,Diffusion, GasExchange] = PerformSegmentation(Proton,Ventilation,Diffusion,GasExchange,MainInput)
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
        LungMask = [];
        switch MainInput.SegmentMethod
            case 'threshold'
                LungMask = Segmentation.SegmentLungthresh(Image,MainInput.SE,MainInput.thresholdlevel);
            case 'manual'
                switch MainInput.SegmentManual
                    case 'basic'
                        LungMask = Segmentation.SegmentLungParenchyma(Image);
                    case 'advance'
%                         VolumSegmt = volumeSegmenter(Image);
                        volumeSegmenter(Image);
%                         uiwait(VolumSegmt)
                        cd(MainInput.XeDataLocation)
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
                Ventilation.LungMask = double(LungMask);
            case 'Diffusion'
                Diffusion.LungMask = double(LungMask);
            case 'GasExchange'
                GasExchange.LungMask = double(LungMask);
        end                
    case 'airway'
        AirwayMask = [];
        AirwayMask = double(Segmentation.SegmentAirway(Image));
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

