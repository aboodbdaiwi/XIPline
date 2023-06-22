%VOLUME This object just holds an N-dimmensional volume (a N-D space) and
%       some basic properties of the volume. The benefit to this class is
%       that its an object, so it can be passed by refference to functions
%       (saves memory).
%
%   Copyright: 2012 Scott Haile Robertson.
%   Website: www.ScottHaileRobertson.com
%   $Revision: 1.1 $  $Date: 2012/04/04 $
%       4/5/2013 Added dimmension names
classdef Volume < handle
    properties (SetAccess = private)
        Data;
        Dims;
        NDims=0;
        DimNames;
        ViewPoint;
        Links = {};
        NLinks = 0;
    end
    
    methods
        function this = Volume(varargin)
            if(nargin > 0)
                this = Global.Volume();
                this.NLinks = 0;
                this = this.updateData(varargin{1});
                if(nargin > 1)
                    this = this.setDimmensionNames(varargin{2});
                else
                    dim_names = cell(1,this.NDims);
                    for i=1:this.NDims
                        dim_names{i} = ['Dim' num2str(i)];
                    end
                    this = this.setDimmensionNames(dim_names);
                end
            end
        end
        
        function this = updateData(this, newData)
            this.Data = newData;
            this.Dims = size(this.Data);
            this.NDims = ndims(this.Data);
        end
        
        function this = setDimmensionNames(this, newDimNames)
            if(iscell(newDimNames))
                if(length(newDimNames)==this.NDims)
                    for i=1:this.NDims
                        if(~isa(newDimNames{i},'char'))
                            if(isnumeric(newDimNames{i}))
                                newDimNames{i} = num2str(newDimNames{i});
                            else
                                error('Dimmension names must be strings or numbers.');
                            end
                        end
                    end
                    this.DimNames = newDimNames;
                else
                    error(['Number of dimmension names must equal number of '...
                        'dimmensions in volume']);
                end
            else
                error('Dimmension names must be passed in as a cell.');
            end
        end
        
        function this = setViewPoint(this, viewPoint)
            this.ViewPoint = viewPoint;
        end
        
        function this = linkToVolume(this, linkVol)
            if(isa(linkVol, 'Volume'))
                % See if you're already linked
                if(~this.isLinked(linkVol))
                    % link volume to this volume
                    this.NLinks = this.NLinks + 1;
                    this.Links{this.NLinks} = linkVol;
                end
                
                % link both ways
                if(~(linkVol.isLinked(this)))
                    linkVol = linkVol.linkToVolume(this);
                end
            else
                error('You can only link to another volume.');
            end
        end
        
        function isLinked = isLinked(this, linkVol)
            isLinked = 0; %false
            for i=1:this.NLinks
                if(isequal(this.Links{i},linkVol))
                   isLinked = 1; %true
                   return;
                end
            end
        end
    end
end