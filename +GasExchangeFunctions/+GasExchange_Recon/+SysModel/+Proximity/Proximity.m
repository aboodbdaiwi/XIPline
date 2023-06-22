%% PROXIMITY
% 
%   An abstract class defining a metric of proximity. There are different 
%   ways to measure "nearness" in k-space, so this serves as a generic
%   object that can have a specific implementation that the user doesn't
%   need to know about
%
% Author: Scott Haile Robertson
% Website: www.ScottHaileRobertson.com
%
classdef (Abstract) Proximity
	properties
		kernel;
        verbose;
		unique_string;
	end
	
	methods
        % Constructor
		function obj = Proximity(kernObj, verbosity)
			% Save properties to object
			obj.kernel = kernObj;
            obj.verbose = verbosity;
		end
	end
	
	methods (Abstract)
		[kernel_vals idxOutOfBounds] = evaluateKernel(obj, traj, matrixSize);
	end
end