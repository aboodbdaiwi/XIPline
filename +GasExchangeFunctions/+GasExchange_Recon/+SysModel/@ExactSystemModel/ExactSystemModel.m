%% EXACTSYSTEMMODEL
%
%   An exact system model class that calculates the DFT of the data. Note
%   that this can be an extremely slow computation, but it uses available
%   parallel resources, and serves as a nice gold standard.
%
% Author: Scott Haile Robertson
% Website: www.ScottHaileRobertson.com
%
classdef ExactSystemModel < GasExchange_Recon.SysModel.SystemModel
    properties
        vox;
        traj;
        isDFT = true;
        loopVoxels = false;
    end
    
    methods
        % Constructor
        function obj = ExactSystemModel(traj, overgridFactor, imageSize, verbosity)
            % Call super constructor to build obj
            obj = obj@GasExchange_Recon.SysModel.SystemModel(overgridFactor, imageSize, [], verbosity);
            
            % Fill in unique string
            obj.unique_string = ['ExactMod_o' num2str(overgridFactor) 'x_sx' num2str(imageSize(1))];
            
            % Make voxel index vector
            lsp_x = [1:obj.fullsize(1)] - 0.5*obj.fullsize(1);
            lsp_y = [1:obj.fullsize(2)] - 0.5*obj.fullsize(2);
            lsp_z = [1:obj.fullsize(3)] - 0.5*obj.fullsize(3);
            [x y z] = meshgrid(lsp_x, lsp_y, lsp_z);
            clear lsp_x lsp_y lsp_z;
            obj.vox = [x(:) y(:) z(:)];
            obj.traj = traj;
            clear x y z;
            
            if(size(obj.traj,1) > size(obj.vox,1))
                obj.loopVoxels = true;
            end
        end
        
        % c = a*b
        function c = mtimes(obj,b)
            if ~isa(obj, 'GasExchange_Recon.SysModel.ExactSystemModel')
                error('Using mtimes directly is not supported.');
            end
            
            if(obj.isDFT)
                % Perform dft - calculate kspace
                % Initialize output
                c = zeros([size(obj.traj,1) 1]);
                
                % Precalculate what we can
                mult = -2*pi*1i;
                
                if(obj.loopVoxels)
                    nSamp = size(obj.vox,1);
%                     p = Progress(nSamp,1);
                    parfor iS = 1:nSamp
                        c = c + b(iS).*exp(mult*(...
                            obj.traj(:,1)*obj.vox(iS,1) + ...
                            obj.traj(:,2)*obj.vox(iS,2) + ...
                            obj.traj(:,3)*obj.vox(iS,3)));
%                         p.progress();
                    end
                else
                    nSamp = size(obj.traj,1);
%                     p = Progress(nSamp,1);
                    parfor iS = 1:nSamp
                        c(iS) = c(iS) + sum(b.*exp(mult*(...
                            obj.traj(iS,1)*obj.vox(:,1) + ...
                            obj.traj(iS,2)*obj.vox(:,2) + ...
                            obj.traj(iS,3)*obj.vox(:,3))));
%                         p.progress();
                    end
                end
%                 p.close();
%                 delete(p);
            else
                % perform idft
                % Initialize output
                c = zeros([size(obj.vox,1) 1]);
                
                % Precalculate what we can
                mult = 2*pi*1i;

                if(obj.loopVoxels)
                    nSamp = size(obj.vox,1);
%                     p = Progress(nSamp,1);
                    parfor iS = 1:nSamp
                        c(iS) = c(iS) + sum(b.*exp(mult*(...
                            obj.traj(:,1)*obj.vox(iS,1) + ...
                            obj.traj(:,2)*obj.vox(iS,2) + ...
                            obj.traj(:,3)*obj.vox(iS,3))));
%                         p.progress();
                    end
                else
                    nSamp = size(obj.traj,1);
%                     p = Progress(nSamp,1);
                    parfor iS = 1:nSamp
                        c = c + b(iS).*exp(mult*(...
                            obj.traj(iS,1)*obj.vox(:,1) + ...
                            obj.traj(iS,2)*obj.vox(:,2) + ...
                            obj.traj(iS,3)*obj.vox(:,3)));
%                         p.progress();
                    end
                end
%                 p.close();
%                 delete(p);
                
                % Revert back
                obj.isDFT = true;
            end
        end
        
        % b = a'
        function b = ctranspose(obj)
            obj.isDFT = false;
            b = obj;
        end
    end
end