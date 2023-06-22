%% TRIANGLE
%
%   A triangle class of gridding kernel. 
%   inputs: kernelExtent    - The nonzero range of the kernel in units of 
%                             pre-overgridded k-space voxels
%           verbose         - If 1, it will verbosely print information
%
%   Author: Scott Haile Robertson
%   Website: www.ScottHaileRobertson.com
%
classdef Triangle < GasExchange_Recon.SysModel.Kernel.Kernel
	methods
		% Constructor
		function obj = Triangle(kernelExtent, verbose)
			% Call super constructor to build obj
			obj = obj@GasExchange_Recon.SysModel.Kernel.Kernel(kernelExtent, verbose);
			
			% Fill in unique string
			obj.unique_string = ['Tri_e' num2str(obj.extent)];
		end
		
		function [kernel_vals] = evaluate(obj, kdistance_preovergrid)
			% Calculate Normalized Gaussian Function
            kernel_vals = 1-(kdistance_preovergrid/obj.extent);
		end
	end
end