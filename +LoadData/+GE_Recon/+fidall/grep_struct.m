function grep_struct(h,str,substruct,fname)
%GREP_STRUCT  Search (compare) strings in structure (fields + values)
%   grep_struct(h,str,substruct,fname)
%         h  Structure for searching
%       str  Search string
% substruct  String containing sub structure (used for recursion)
%     fname  Write output to filename not stdout                       ('')
%
%  9/2012 Martin Janich
% 11/2023 Rolf Schulte
if (nargin<2), help(mfilename); return; end


%% input
if isnumeric(str),            str = num2str(str); end
if ~exist('substruct','var'), substruct = []; end
if isempty(substruct),        substruct = inputname(1); end
if ~exist('fname','var'),     fname = ''; end


%% redirect output to file
fid = 1;
if ~isempty(fname)
    if isnumeric(fname)
        fid = fname;
    else
        if ischar(fname)
            fid = fopen(fname,'wt');
        else
            warning('fname');
        end
    end
end


%% go through entries
fn1 = fieldnames(h);
for l1=1:length(fn1)
    xx1 = getfield(h,fn1{l1});
    if isstruct(xx1)
        grep_struct(xx1,str,[substruct '.' fn1{l1}],fid);
    else
        if ~isempty(regexpi(fn1{l1},str, 'once')) || isempty(str)
            if ~iscell(xx1)
                try
                    tmpstr = mat2str(xx1);
                catch 
                    tmpstr = mat2str(xx1(:));
                end
                fprintf(fid,'%s.%s = %s\n',substruct,fn1{l1},tmpstr);
            else
                fprintf(fid,'%s.%s{1} = \n',substruct,fn1{l1});
                disp(xx1{1});
            end
        end
        
        % if isrow(xx1)
        nn1 = size(xx1);
        if (length(nn1)<3) && (nn1(1)==1)
            if iscell(xx1)
                for l2=1:nn1(2)
                    if ischar(xx1{l2})
                        if ~isempty(regexpi(xx1{l2},str, 'once')) || isempty(str)
                            fprintf(fid,'%s.%s = {''%s''}\n',substruct,fn1{l1},xx1{l2});
                        end
                    else
                        grep_struct(xx1{l2},str,...
                            [substruct '.' fn1{l1} '{' num2str(l2) '}'],fid)
                    end
                end
            else
                if ~isempty(regexpi(num2str(xx1),str, 'once'))
                    fprintf(fid,'%s.%s = %s\n',substruct,fn1{l1},num2str(xx1));
                end
            end
        end
    end
end


%% close file
if ~isempty(fname) && ischar(fname)
    fclose(fid);
end


end      % grep_struct.m
