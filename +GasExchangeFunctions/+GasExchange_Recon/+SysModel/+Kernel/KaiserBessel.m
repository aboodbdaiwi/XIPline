%% KAISERBESSEL
%   A Kaiser-Bessel class of gridding kernel.
%   creates a Kaiser Bessel function from 0-kernel_extent_oversamp. 

%   Jackson et al. showed that the Kaiser Bessel function is a good 
%   approximation to the prolate spheroidal wave function:
%       Selection of a Convolution Function for Fourier Inversion using 
%       Gridding. Jackson et al. 1991.
%
%   Choice of Beta from: 
%       Rapid Gridding Reconstruction with a minimal Oversampling ratio" 
%       Beaty et all IEEE 2005.
%
%   inputs: beta            - The sharpness of the Kaiser-Bessel function
%           kernelExtent    - The nonzero range of the kernel in units of 
%                             pre-overgridded k-space voxels
%           verbose         - If 1, it will verbosely print information
%
% Author: Scott Haile Robertson
% Website: www.ScottHaileRobertson.com
%
classdef KaiserBessel < GasExchange_Recon.SysModel.Kernel.Kernel
	properties
		beta;
	end
	
	methods
		% Constructor
		function obj = KaiserBessel(kernelBeta, kernelExtent, verbose)
			% Call super constructor to build obj
			obj = obj@GasExchange_Recon.SysModel.Kernel.Kernel(kernelExtent, verbose);
			
			% Store properties
			obj.extent = kernelExtent;
			
			%Calculate Beta value if not given
			if(isempty(kernelBeta))
                % use value from:
                % Rapid Gridding Reconstruction With a Minimal Oversampling 
                % Ratio. Beatty et al. 2005.)
				obj.beta = pi*sqrt( (0.5*obj.extent)^2-0.8 );
			else
				obj.beta = kernelBeta;
            end
            
			% Fill in unique string
			obj.unique_string = ['kb_e' num2str(obj.extent) ...
				'_s' num2str(obj.beta)];
		end
		
		function [kernel_vals] = evaluate(obj, distances)	
			% Calculate Kaiser Bessel Function
			kernel_vals = besseli(0,obj.beta*...
				sqrt(1 - (2*distances/obj.extent).^2))./obj.extent;
			
			% Its tempting to just use some fixed size (like the matrix
			% size), however the kaiser bessel function is created to be
			% zero at the kernel size, thus its best to just use the kernel
			% size. This makes optimization a pain, but theres not much we
			% can do...
			kernel_vals = besseli(0,obj.beta*...
				sqrt(1 - (2*distances/obj.extent).^2))./obj.extent;

			%Normalize - according to maximum of function
            max_val =  besseli(0,obj.beta)/obj.extent;
			kernel_vals = kernel_vals/max_val;
		end
	end
end