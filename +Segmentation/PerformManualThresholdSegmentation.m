function [Proton,Ventilation,Diffusion, GasExchange] = PerformManualThresholdSegmentation(Proton,Ventilation,Diffusion,GasExchange,MainInput)
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
        switch MainInput.Imagestosegment
            case 'Xe & Proton Registered'
               Image_to_Segment = Proton.ProtonRegisteredColored;
            case 'Xenon'
               Image_to_Segment = Ventilation.Image;               
            case 'Registered Proton'
               Image_to_Segment = Proton.ProtonRegistered;  
        end 

    case 'Diffusion'
        switch MainInput.Imagestosegment
            case 'Xe & Proton Registered'
               Image_to_Segment = Diffusion.Image; 
            case 'Xenon'
               Image_to_Segment = Diffusion.Image;             
            case 'Registered Proton'
               Image_to_Segment = Diffusion.Image; 
        end

    case 'GasExchange'
        Image = GasExchange.VentImage;
        switch MainInput.Imagestosegment
            case 'Xe & Proton Registered'
               Image_to_Segment = GasExchange.ProtonRegisteredColored;
            case 'Xenon'
               Image_to_Segment = GasExchange.VentImage;               
            case 'Registered Proton'
               Image_to_Segment = GasExchange.ProtonRegistered;  
        end
end

% perform segmention
switch MainInput.SegmentAnatomy
    case 'Parenchyma'
        LungMask = [];
        switch MainInput.SegmentationMethod
            case 'Threshold'
                LungMask = Segmentation.SegmentLungthresh(Image_to_Segment,MainInput.SE,MainInput.thresholdlevel);
            case 'Manual'
                switch MainInput.SegmentManual
                    case 'Freehand'
                        LungMask = Segmentation.SegmentLungParenchyma(Image_to_Segment);
                    case 'AppSegmenter'
%                         VolumSegmt = volumeSegmenter(Image);
                        volumeSegmenter(Image_to_Segment);
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
    case 'Airway'
        AirwayMask = [];
        AirwayMask = double(Segmentation.SegmentAirway(Image_to_Segment));
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