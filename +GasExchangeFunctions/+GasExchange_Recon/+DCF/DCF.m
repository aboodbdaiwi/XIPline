%% DCF
%
%   An abstract class describes how to calculate a density compensation
%   function (DCF) as well as store the DCF coefficients. 
%
% Author: Scott Haile Robertson
% Website: www.ScottHaileRobertson.com
classdef (Abstract) DCF
	properties
		dcf;
		unique_string;
		space;
		verbose;
    end
    methods
        % Constructor
        function obj = DCF(verbosity)           
			% Store properties of DCF
			obj.verbose = verbosity;
        end
        
        % c = a.*b
        function c = times(a,b)
            if(isa(a, 'GasExchangeFunctions.GasExchange_Recon.DCF.DCF'))
                c = a.dcf.*b;
            elseif(isa(b, 'GasExchangeFunctions.GasExchange_Recon.DCF.DCF'))
                c = b.dcf.*a;
            else
                error('Using times directly is not supported.');
            end
        end
    end
end
