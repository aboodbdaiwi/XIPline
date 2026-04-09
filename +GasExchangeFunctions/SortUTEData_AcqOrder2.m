function sortedData = SortUTEData_AcqOrder2(unsortedData, DataSizes, order)

tempData   = zeros(size(unsortedData), 'like', unsortedData);
sortedData = zeros(size(unsortedData), 'like', unsortedData);

PtsPerKz = DataSizes(2);   % 122
NKz      = DataSizes(3);   % 10
NPro     = PtsPerKz * NKz;
Pts      = 1:PtsPerKz:NPro;
NewPts   = 1:NKz;

for i = 1:PtsPerKz
    tempData(:, NewPts + NKz*(i-1)) = unsortedData(:, Pts + (i-1));
end

for i = 1:NPro
    sortedData(:, i) = tempData(:, order(i));
end

end