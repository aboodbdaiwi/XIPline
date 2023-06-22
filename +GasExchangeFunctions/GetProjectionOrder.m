function [order, ky, kz] = GetProjectionOrder(ListInfo)
%OrderProjections - Parses the .list info from the loadLISTDATA function to
%provide an order of the projections acquired with respect to time

del_indx = [];
%Delete all indices that don't correspond to the dissolved phase data
for i = 1:length(ListInfo)
    if strcmp(ListInfo(i).typ,'NOI')
        del_indx = [del_indx i];
    end
    if ListInfo(i).mix ~= 1
        del_indx = [del_indx i];
    end
    if ListInfo(i).dyn ~=1
        del_indx = [del_indx i];
    end
end
del_indx = unique(del_indx);
ListInfo(del_indx) = [];

%Pull ky and kz from the structures
ky = zeros(1,length(ListInfo));
kz = zeros(1,length(ListInfo));
for i = 1:length(ListInfo)
    ky(i) = ListInfo(i).ky;
    kz(i) = ListInfo(i).kz;
end
%Calculate what the projection number is
order = max(kz) * (ky-1)+kz;

end