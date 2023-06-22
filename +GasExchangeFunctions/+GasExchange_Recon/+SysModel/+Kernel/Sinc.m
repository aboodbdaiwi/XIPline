%% Sinc
%   A sinc class of gridding kernel.
%
%   inputs: firstZero       - location of the first sinc zero
%           kernelExtent    - The nonzero range of the kernel in units of 
%                             pre-overgridded k-space voxels
%           verbose         - If 1, it will verbosely print information
%
% Author: Scott Haile Robertson
% Website: www.ScottHaileRobertson.com
%
classdef Sinc < GasExchange_Recon.SysModel.Kernel.Kernel
	properties
        firstzero;
	end
	
	methods
		% Constructor
		function obj = Sinc(firstZero, kernelExtent, verbose)
			% Call super constructor to build obj
			obj = obj@GasExchange_Recon.SysModel.Kernel.Kernel(kernelExtent, verbose);
			
			% Store properties
			obj.firstzero = firstZero;
						
			% Fill in unique string
			obj.unique_string = ['sinc_e' num2str(obj.extent) ...
				'_s' num2str(obj.firstzero)];
		end
		
		function [kernel_vals] = evaluate(obj, distances)
			% Calculate sinc Function
			kernel_vals = sin(2*pi*distances/obj.firstzero)./(2*pi*distances/obj.firstzero);
			kernel_vals(distances==0)=1;
		end
	end
end