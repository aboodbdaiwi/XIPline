function [order, ky, kz] = GetProjectionOrder2(ListInfo)

keep = false(1,length(ListInfo));

for i = 1:length(ListInfo)
    if strcmp(strtrim(ListInfo(i).typ),'STD') && ...
       ListInfo(i).mix == 1 && ...
       ListInfo(i).dyn == 1
        keep(i) = true;
    end
end

ListInfo = ListInfo(keep);

ky = [ListInfo.ky];
kz = [ListInfo.kz];

order = (ky - 1) .* max(kz) + kz;
end