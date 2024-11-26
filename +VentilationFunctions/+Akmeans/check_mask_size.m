function temp_handles=check_mask_size(temp_handles,check_fields)
%
%
% Copyright W. Zha @2014
%%
if nargin<2
check_fields = {'lobemask','lungsmask',...             
                'lobemaskLLL','lobemaskLUL',...
                'lobemaskRUL','lobemaskRML','lobemaskRLL','defectmask',...
                'lobemask2','lungsmask2',...             
                'lobemask2LLL','lobemask2LUL',...
                'lobemask2RUL','lobemask2RML','lobemask2RLL','defectmask2'};
end 
% Use 3He as the reference
[nRows,nCols,nSlices] = size(temp_handles.he);

for n_field = 1:numel(check_fields)
    this_field = check_fields{n_field};
    if isfield(temp_handles,this_field)
        nfield_slices = size(temp_handles.(this_field),3);
        if nfield_slices>=nSlices
            temp_handles.(this_field)=temp_handles.(this_field)(:,:,1:nSlices);
        else
            
            mask_tmp = zeros(nRows,nCols,nSlices);
            mask_tmp(:,:,1:nfield_slices) = temp_handles.(this_field);
            temp_handles.(this_field)=mask_tmp;
        end
    end
end