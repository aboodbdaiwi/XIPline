%% KERNEL
% 
%   An abstract class of gridding kernel. 
%
% Author: Scott Haile Robertson
% Website: www.ScottHaileRobertson.com
%
classdef (Abstract) Kernel
	properties
		verbose; 
        extent;
		unique_string;
	end
	
	methods
        % Constructor
		function obj = Kernel(kernelExtent, verbosity)
			% Save properties of object
			obj.verbose = verbosity;
            
            obj.extent = kernelExtent;
		end
	end
	
	methods (Abstract)
		 kernel_vals = evaluate(obj, distances);	
	end
end