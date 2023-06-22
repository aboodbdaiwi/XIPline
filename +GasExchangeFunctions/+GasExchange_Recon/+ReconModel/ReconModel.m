classdef ReconModel
	properties
		verbose;
		unique_string;
        system;
        deapodize = true;
        crop = true;
        PixelShift = [0, 0 ,0]; %MMW
	end
	
	methods
		function obj = ReconModel(system_model, verbosity)
			obj.verbose = verbosity;
			obj.system = system_model;
        end
    end
    
    methods (Abstract)
		% Reconstructs an image volume using the given data
		reconVol = reconstruct(obj,data, traj);
	end
end
