function h = header_anonymise(h,rmflds)
%HEADER_ANONYMISE Remove patient data from p-file header
%     h = header_anonymise(h)
%     h   p-file header
%rmflds   fields to zero 
%         default={'exam.patnameff','exam.dateofbirth','exam.patidff'};
%
%  02/2021 Rolf Schulte
if (nargin<1), help(mfilename); return; end

if ~exist('rmflds','var'), rmflds = []; end
if isempty(rmflds)
    rmflds = {'exam.patnameff','exam.dateofbirth','exam.patidff'};
end

fprintf('Anonymising:\n');
for l=1:length(rmflds)
    fld = rmflds{l};
    fprintf('%g: %s\n',l,fld);
    h = sub_null_field(h,fld);
end

end      % header_anonymise.m



%% call routine recursively to null field
function h = sub_null_field(h,fld)

dots = regexp(fld,'\.');
if ~isempty(dots)
    if isfield(h,(fld(1:(dots(1)-1))))
        h.(fld(1:(dots(1)-1))) = sub_null_field(h.(fld(1:(dots(1)-1))),fld((dots(1)+1):end));
    else
        warning('field %s not found (1)',(fld(1:(dots(1)-1))));
    end
else
    if isfield(h,(fld))
        h.(fld) = ''; 
    else
        warning('field %s not found (2)',fld);
    end
end

end      % sub_null_field
