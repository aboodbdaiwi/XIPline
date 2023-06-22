%% UNITY
%
%   An unity DCF class treats all samples as having the same (unit) density
%
%   Author: Scott Haile Robertson
%   Website: www.ScottHaileRobertson.com
%
classdef Unity < GasExchange_Recon.DCF.DCF
	methods
		% Constructor
		function obj = Unity(traj, verbosity)      
            % Call super constructor to build obj
            obj = obj@GasExchange_Recon.DCF.DCF(verbosity);
            
			% Store properties of DCF
			obj.unique_string = 'unityDcf';
			obj.space = 'dataspace';
			
			obj.dcf = ones(size(traj,1),1);
		end
	end
end
