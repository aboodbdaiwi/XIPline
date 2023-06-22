%% MATRIXSYSTEMMODEL
%
%   A matrix system model class that stores all interpolation
%   coefficients into a sparse matrix. Note that storage of the
%   interpolation coefficients can take significant memory. If you are
%   memory limited, consider creating a class which calculates
%   interpolation coefficients on the fly, and does not require the memory
%   overhead for the system matrix. The downside to on the fly calculations
%   is that they compute slower in itterative applications, where
%   interpolation coefficients are calculated twice each iteration (once to
%   grid, and once to ungrid)
%
% Author: Scott Haile Robertson
% Website: www.ScottHaileRobertson.com
%
classdef MatrixSystemModel < GasExchangeFunctions.GasExchange_Recon.SysModel.SystemModel
    properties
        A;
        ATrans;
        supersparse;
        issupersparse;
        istranspose = false;
    end
    methods
        % Constructor
        function obj = MatrixSystemModel(traj, overgridFactor, imageSize, proximityObj, ...
                verbosity)
            % Call super constructor to build obj
            obj = obj@GasExchangeFunctions.GasExchange_Recon.SysModel.SystemModel(overgridFactor, imageSize, proximityObj, verbosity);
            
            % Fill in unique string
            obj.unique_string = ['MatMod_' proximityObj.unique_string];
            obj.issupersparse = false;
            
            % Calculate interpolation coefficients of system matrix
            if(obj.verbose)
                disp('Calculating Matrix interpolation coefficients...');
            end
            [sample_idx,voxel_idx,kernel_vals] = obj.proximity.evaluateKernel(traj,...
                obj.overgrid, obj.fullsize);
            if(obj.verbose)
                disp('Finished calculating Matrix interpolation coefficients.');
            end
            
            % Make sparse matrix
            obj.A = sparse(sample_idx,voxel_idx,kernel_vals,size(traj,1),...
                prod(obj.fullsize),length(sample_idx));
            obj.ATrans = obj.A';
        end
        
        function varargout = makeSuperSparse(obj, varargin)
            if(obj.verbose)
                disp(['   Making System Model super-sparse...']);
            end
            % Keep track of non-sparse voxels (non-empty collumns)
            obj.supersparse.fullsize = size(obj.A);
            [nzrows obj.supersparse.nzcols nzvals] = find(obj.A);
            
            % Check if its worth making matrix more sparse
            if(nzvals < 0.75*prod(obj.supersparse.fullsize))
                % Find unique non empty colums for reindexing
                [obj.supersparse.uniqueNzCols,iCol,iUniqueCol] = unique(obj.supersparse.nzcols);
                
                % Make the system matrix more sparse
                obj.A = sparse(nzrows,iUniqueCol,nzvals,...
                    obj.supersparse.fullsize(1), ...
                    length(obj.supersparse.uniqueNzCols),...
                    length(nzvals));
                
                obj.ATrans = obj.A';
                
                % Change sparseness flag
                obj.issupersparse = true;
                
                if(obj.verbose)
                    disp(['   Finished making System Model super-sparse.']);
                end
            else
                obj.issupersparse = false;
                if(obj.verbose)
                    disp(['   Not making System Model sparser, its very sparse already.']);
                end
            end
            varargout{1} = obj;
            
            if(obj.issupersparse)
                % If additional arguments are provided, we need to make those
                % matrices sparse too 
                nVols = length(varargin);
                for iVol = 1:nVols
                    if(obj.verbose)
                        disp(['   Making ' num2str(iVol) '/' num2str(nVols) ' additional volume(s) super-sparse.']);
                    end
                    varargout{iVol+1} = varargin{iVol}(obj.supersparse.uniqueNzCols);
                end
            end
        end
        
        function varargout = revertSparseness(obj, varargin)
            if(obj.issupersparse)
                if(obj.verbose)
                    disp(['   Reverting from super-sparse System Matrix...']);
                end
                % increase matrix dimensions
                tmpSparse = sparse([],[],[],obj.supersparse.fullsize(1),...
                    obj.supersparse.fullsize(1),nnz(obj.A));
                
                % Copy values into new matrix
                tmpSparse(:,obj.supersparse.uniqueNzCols) = obj.A(:,:);
                
                % apply less sparse matrix to object, remove supersparse
                % data and revert sparseness flag
                obj.A = tmpSparse;
                obj.ATrans = obj.A';
                obj.issupersparse = false;
                clear tmpSparse;
                if(obj.verbose)
                    disp(['   Finished reverting from super-sparse System Matrix.']);
                end
                varargout{1} = obj;
                
                % If additional arguments are provided, we need to revert 
                % the sparseness of those matrices too
                nVols = length(varargin);
                for iVol = 1:nVols
                    if(obj.verbose)
                        disp(['   Making ' num2str(iVol) '/' num2str(nVols) ' additional volume(s) super-sparse.']);
                    end
                    varargout{iVol+1} = zeros(obj.supersparse.fullsize(2),1);
                    varargout{iVol+1}(obj.supersparse.uniqueNzCols) = varargin{iVol};
                end
                
                % remove supersparse
                obj.supersparse = [];
            end
        end
        
        % c = a*b
        function c = mtimes(obj,b)         
            if(obj.istranspose)
%                 c = (b'*obj.A)'; % Dont transpose A!
                c = obj.ATrans*b; % Dont transpose A!
            else
                c = obj.A*b;
            end
        end
        
        % b = a'
        function obj = ctranspose(obj)
            obj.istranspose = ~obj.istranspose;
        end
    end
end