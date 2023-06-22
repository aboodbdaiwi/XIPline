%% ANALYTICAL3DRADIAL
%
%   An analytical DCF class based off the assumption that the radial rays
%   perfectly diverge according to kr^-2
%
%   Author: Scott Haile Robertson
%   Website: www.ScottHaileRobertson.com
%
classdef Analytical3dRadial < GasExchange_Recon.DCF.DCF
	methods
		function obj = Analytical3dRadial(traj, verbose)
            % obj = Analytical3dRadial(traj, verbosity)
            %   traj = k-space trajectory (normalized to -0.5 and +0.5 for
            %          upper and lower Nyquist limits, respectively
            %   verbose = Show verbose print outs?
            
            % Call super constructor to build obj
            obj = obj@GasExchange_Recon.DCF.DCF(verbose);
            
			% Store properties of DCF
			obj.unique_string = 'analyt';
			obj.space = 'dataspace';
			
            % Calculate squared radial position
			obj.dcf = sum(traj.^2,2);			
		end
	end
end
