classdef ConjugateGradient < GasExchange_Recon.ReconModel.GriddedReconModel
    properties
        iterations;
        weights;
    end
    methods
        function obj = ConjugateGradient(system_model, iter, verbose)
            % Call super constructor to build recon obj
            obj = obj@GasExchange_Recon.ReconModel.GriddedReconModel(system_model, verbose);
            
            % Save properties
            obj.deapodize = 0;
            obj.iterations = iter;
            obj.unique_string = ['cgiter' num2str(obj.iterations) '_' obj.system.unique_string];
        end
       
        function gridVol = grid(obj, data, traj, varargin)
            % Reconstructs an image volume from the given data using the
            % conjugate gradient algorithm
            
            %Parse inputs
            if(nargin<3 || isempty(varargin{1}))
                % If initial guess isnt given, assume zeros
                gridVol = zeros(obj.system.fullsize);
            else
                % first input is initial guess
                gridVol = varargin{1};
            end
            
            % Vectorize initial guess
            gridVol = gridVol(:);
            
            % Be clever and only calculate non-sparse k-space voxels by
            % using a sparser system model. This way you dont grid or
            % ungrid zeros unnecessarily
            if(isa(obj.system,'GasExchange_Recon.SysModel.MatrixSystemModel') )
                [obj.system gridVol] = obj.system.makeSuperSparse(gridVol);
            end
            
            if(obj.verbose)
                disp('   Starting Conjugate Gradient Reconstruction...');
            end    
            [gridVol,flag,relres,iter,resvec] = lsqr(@aFun, data, 0, ...
                obj.iterations,[],[],gridVol);
            if(obj.verbose)
                disp('   Finished Conjugate Gradient Reconstruction.');
            end   
            
            % Be clever and only calculate non-sparse k-space voxels by
            % using a sparser system model. This way you dont grid or
            % ungrid zeros unnecessarily
            if(isa(obj.system,'GasExchange_Recon.SysModel.MatrixSystemModel'))
                [obj.system gridVol] = obj.system.revertSparseness(gridVol);
            end 
                        
            function y = aFun(x_,transp_flag)
                if strcmp(transp_flag,'transp')      % x = A'*y
                    y = (obj.system'*x_);
                elseif strcmp(transp_flag,'notransp') % y = A*x
                    y = (obj.system*x_);
                end
            end
            
        end
    end
end
