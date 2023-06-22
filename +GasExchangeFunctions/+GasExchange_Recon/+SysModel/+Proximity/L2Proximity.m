%% L2PROXIMITY
%   An L2 proximity class defining distance in an L2 sense also known as 
%   the Euclidean/pythagorean distance.
%
%   inputs: kernelObj       - A kernel object for ecaluating the kernel
%           verbose         - If 1, it will verbosely print information
%
% Author: Scott Haile Robertson
% Website: www.ScottHaileRobertson.com
%
classdef L2Proximity < GasExchangeFunctions.GasExchange_Recon.SysModel.Proximity.Proximity
	methods
		function obj = L2Proximity(kernelObj, verbosity)
			% Call super constructor to build obj
			obj = obj@GasExchangeFunctions.GasExchange_Recon.SysModel.Proximity.Proximity(kernelObj, verbosity);
			
			% Save properties of proximity object and contained kernel
			obj.unique_string = ['L2_'  kernelObj.unique_string];
		end

		function [sample_idx,voxel_idx,kernel_vals] = ...
                evaluateKernel(obj, traj, overgridFactor, matrixSize)
            
			% Calculate gridding in pre-overgridding distances
			if(obj.verbose)
				disp('Calculating L2 distances...');
            end
            
			[sample_idx,voxel_idx,preovergrid_distances] = ...
				GasExchangeFunctions.GasExchange_Recon.SysModel.Proximity.sparse_gridding_distance_mex(traj',...
				obj.kernel.extent*overgridFactor,...
				uint32(matrixSize'), -1); % Calculate accross all dimensions
            preovergrid_distances = preovergrid_distances/overgridFactor;
			if(obj.verbose)
				disp('Finished calculating L2 distances.');
            end
			
			% Look for any values that are still out of bounds
			if(obj.verbose)
				disp('Applying L2 Bounds...');
            end
            keepValues =  (sample_idx > 0) & (voxel_idx > 0);
			sample_idx = sample_idx(keepValues);
			voxel_idx = voxel_idx(keepValues);
			preovergrid_distances = preovergrid_distances(keepValues);
			clear keepValues;
			if(obj.verbose)
				disp('Finished applying L2 Bounds...');
			end
			
			if(obj.verbose)
				disp('Applying Kernel...');
			end
			% Calculate kernel values
			kernel_vals = obj.kernel.evaluate(preovergrid_distances);
			clear preovergrid_distances;
			if(obj.verbose)
				disp('Finished applying Kernel...');
			end
		end
	end
end