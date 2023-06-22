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
classdef L1Proximity < GasExchange_Recon.SysModel.Proximity.Proximity
	methods
		function obj = L1Proximity(kernelObj, verbosity)
			% Call super constructor to build obj
			obj = obj@GasExchange_Recon.SysModel.Proximity.Proximity(kernelObj, verbosity);
			
			% Save properties of proximity object and contained kernel
			obj.unique_string = ['L1_og' kernelObj.unique_string];
        end
		
        function [sample_idx,voxel_idx,kernel_vals] = evaluateKernel(...
                obj, traj, overgridFactor, matrixSize)
			
            % Calculate gridding in pre-overgridding distances
			if(obj.verbose)
				disp('Calculating L1 distances...');
            end
			[nSamps nDims] = size(traj);
			keepValues = true(1,1);
			kernel_vals = 1;
			for iDim = 1:nDims
				if(obj.verbose)
					disp(['  Dim ' num2str(iDim) '...']);
                end
							
                % Only calculate for this dimension (c uses 0-indexing!)
				[sample_idx,voxel_idx,distDim] = ...
					GasExchange_Recon.SysModel.Proximity.sparse_gridding_distance_mex(traj',...
					obj.kernel.extent*overgridFactor,...
					uint32(matrixSize'), iDim-1);
				
				% Look for any values that are still out of bounds
				keepValues = keepValues & (sample_idx > 0) & (voxel_idx > 0);
				
				% Accumulate kernel values (saves memory this way for >2D)
                if(obj.verbose)
					disp(['  Evaluating kernel for dim ' num2str(iDim) '...']);
                end
				kernel_vals = kernel_vals.* obj.kernel.evaluate(distDim);
            end
			if(obj.verbose)
				disp('Finished calculating L1 distances.');
            end
			
            % Look for any values that are still out of bounds
			if(obj.verbose)
				disp('Applying L1 Bounds...');
            end
			sample_idx = sample_idx(keepValues);
			voxel_idx = voxel_idx(keepValues);
			kernel_vals = kernel_vals(keepValues);
			clear keepValues;
			if(obj.verbose)
				disp('Finished applying L1 Bounds...');
			end
		end
	end
end