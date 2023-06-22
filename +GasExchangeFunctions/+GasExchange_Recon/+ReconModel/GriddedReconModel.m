classdef GriddedReconModel < GasExchangeFunctions.GasExchange_Recon.ReconModel.ReconModel
    methods
        % Constructor
        function [obj] = GriddedReconModel(systemModel, verbose)
            % Call super constructor to build recon obj
            obj = obj@GasExchangeFunctions.GasExchange_Recon.ReconModel.ReconModel(systemModel, verbose);
        end
        
        function reconVol = reconstruct(obj, data, traj)
            if(obj.verbose)
                disp('Reconstructing...');
            end
            
            % Grid data
            if(obj.verbose)
                disp('Gridding Data...');
            end
            reconVol = obj.grid(data);
            if(obj.verbose)
                disp('Finished gridding Data...');
            end
            
            % Reshape from vector to matrix
            reconVol = reshape(full(reconVol),ceil(obj.system.fullsize));
            
            if(obj.verbose)
                disp('Calculating IFFTN...');
            end
            % Calculate image space
            reconVol = ifftn(reconVol);
            reconVol = fftshift(reconVol);
            if(obj.verbose)
                disp('Finished calculating IFFTN.');
            end
            
            %MMW Start
            reconVol = circshift(circshift(circshift(reconVol,round(obj.PixelShift(1)),3),round(obj.PixelShift(2)),2),round(obj.PixelShift(3)),1);
            %MMW end
            
            if(obj.crop)
                reconVol = obj.system.crop(reconVol);
            end
            
            if(obj.deapodize)
                % Calculate deapodization volume and deapodize
                if(obj.verbose)
                    disp('Calculating k-space deapodization function...');
                end
                deapVol = obj.grid(double(~any(traj,2)));
                %MMW allow deap even for imperfect trajs
                if sum(~any(traj,2)) == 0
                    dist = sqrt(traj(:,1).^2+traj(:,2).^2+traj(:,3).^2);
                    temp = zeros(size(dist));
                    temp(dist<=(min(dist)*1.01)) = 1;
                    deapVol = obj.grid(double(temp));
                end
                
                % Reshape from vector to matrix
                deapVol = reshape(full(deapVol),ceil(obj.system.fullsize));
                if(obj.verbose)
                    disp('Finished calculating k-space deapodization function.');
                end
                
                % Calculate image domain representation
                if(obj.verbose)
                    disp('Calculating Image domain deapodization function...');
                end
                deapVol = ifftn(deapVol);
                deapVol = fftshift(deapVol);
                if(obj.verbose)
                    disp('Calculating Image domain deapodization function...');
                end
                
                %MMW Start
                deapVol = circshift(circshift(circshift(deapVol,round(obj.PixelShift(1)),3),round(obj.PixelShift(2)),2),round(obj.PixelShift(3)),1);
                %MMW end
                
                if(obj.crop)
                    deapVol = obj.system.crop(deapVol);
                end
                
                if(obj.verbose)
                    disp('deapodizing...');
                end
                reconVol = reconVol./deapVol;
                clear deapVol;
                if(obj.verbose)
                    disp('Finished deapodizing.');
                end
            end
            
            if(obj.verbose)
                disp('Finished Reconstructing.');
            end
        end
    end
    methods (Abstract)
        gridVol = grid(obj, data);
    end
end
