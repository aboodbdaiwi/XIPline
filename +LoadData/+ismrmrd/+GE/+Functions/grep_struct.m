function grep_struct(h,str,substruct)
%GREP_STRUCT  Search (compare) strings in structure (fields + values)
%   grep_struct(h,str,substruct)
%         h  Structure for searching
%       str  Search string
% substruct  String containing sub structure (used for recursion)
%  3/2010 Rolf Schulte
%  9/2012 Martin Janich
%  8/2020 Rolf Schulte
if (nargin<2), help(mfilename); return; end
if ~exist('substruct','var'), substruct = inputname(1); end
if isnumeric(str), str = num2str(str); end

fn1 = fieldnames(h);
for l1=1:length(fn1)
    xx1 = getfield(h,fn1{l1});
    if isstruct(xx1)
        grep_struct(xx1,str,[substruct '.' fn1{l1}]);
    else
        if ~isempty(regexpi(fn1{l1},str))
            if ~iscell(xx1)
                fprintf('%s.%s = %s\n',substruct,fn1{l1},mat2str(xx1));
            else
                fprintf('%s.%s{1} = \n',substruct,fn1{l1})
                disp(xx1{1});
            end
        end
        
        % if isrow(xx1)
        nn1 = size(xx1);
        if (length(nn1)<3) && (nn1(1)==1)
            if iscell(xx1)
                for l2=1:nn1(2)
                    if ischar(xx1{l2})
                        if ~isempty(regexpi(xx1{l2},str))
                            fprintf('%s.%s = {''%s''}\n',substruct,fn1{l1},xx1{l2});
                        end
                    else
                        grep_struct(xx1{l2},str,...
                            [substruct '.' fn1{l1} '{' num2str(l2) '}'])
                    end
                end
            else
                if ~isempty(regexpi(num2str(xx1),str))
                    fprintf('%s.%s = %s\n',substruct,fn1{l1},num2str(xx1));
                end
            end
        end
    end
end

end      % main function grep_struct.m
