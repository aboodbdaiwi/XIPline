%% SYSTEMMODEL
%
%   An abstract class defining a particular system model represenation
%
% Author: Scott Haile Robertson
% Website: www.ScottHaileRobertson.com
%
classdef (Abstract) SystemModel
    properties
        proximity;
        cropsize;
        fullsize;
        overgrid;
        unique_string;
        verbose;
    end
    methods
        % Constructor
        function obj = SystemModel(overgridFactor, cropSize, proximityObj, verbosity)
            obj.verbose = verbosity;
            obj.overgrid = overgridFactor;
            obj.proximity = proximityObj;
            obj.cropsize = cropSize;
            obj.fullsize = ceil(obj.cropsize*obj.overgrid);
        end
        
        function croppedVol = crop(obj,uncroppedVol)
            croppedVol = subvolume(uncroppedVol,...
                [round([0.5*(obj.fullsize-obj.cropsize)+1]); ...
                round([0.5*(obj.fullsize+obj.cropsize)])]);
        end
    end
    
    methods (Abstract)
        c = mtimes(obj,b); %a*b
        b = ctranspose(obj); %a'
    end
end