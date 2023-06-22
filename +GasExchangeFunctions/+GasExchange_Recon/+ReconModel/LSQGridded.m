classdef LSQGridded < GasExchangeFunctions.GasExchange_Recon.ReconModel.GriddedReconModel
	properties
		D;
	end
	methods
		% Constructor 
		function [obj] = LSQGridded(systemModel, D, verbose)
			% Call super constructor to build recon obj
			obj = obj@GasExchangeFunctions.GasExchange_Recon.ReconModel.GriddedReconModel(systemModel, verbose);
		
            % Store properties
            obj.D = D;
			obj.unique_string = ['grid_' obj.system.unique_string '_' obj.D.unique_string];
        end
        
        function gridVol = grid(obj, data)
            if(isa(obj.system,'TemporalSystemModel'))
                 switch(obj.D.space)
                    case 'gridspace'
                        error('Not implemented');
                    case 'dataspace'
                        [nSamples nCoils nTimes] = size(data);
                        dcfWeightedData = zeros([nSamples nCoils nTimes]);
                        for iTime=1:nTimes
                            for iCoil=1:nCoils
                                dcfWeightedData(:,iCoil,iTime) = data(:,iCoil,iTime).*obj.D.dcf{iTime};
                            end
                        end
                        gridVol = obj.system' * dcfWeightedData;
                    otherwise
                        error('DCF style not recognized');
                end 
            else
                switch(obj.D.space)
                    case 'gridspace'
                        gridVol = (obj.system' * data) .* obj.D;
                    case 'dataspace'
                        gridVol = obj.system' * (obj.D .* data);
                    otherwise
                        error('DCF style not recognized');
                end
            end
		end
	end
end
