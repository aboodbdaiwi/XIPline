function sortedData = SortUTEData_AcqOrder(unsortedData, DataSizes, order)
%SortUTEData_AcqOrder - sorts data using projection order
%  assumes unsortedData is [read out, proj]
%  assumes DataSize is [read out, profs, kz(interleafs/hubs)]

%initialize output
tempData = zeros(size(unsortedData));
sortedData = zeros(size(unsortedData));

%Calculate params
PtsPerKz = DataSizes(2);
NKz = DataSizes(3);
NPro = PtsPerKz * NKz;
Pts = 1:PtsPerKz:NPro;
NewPts = 1:NKz;

%Sort kz
for i = 1:PtsPerKz
    tempData(:,NewPts + NKz*(i-1)) = unsortedData(:,Pts+(i-1));
end

%Final Sort
for i = 1:NPro
    sortedData(:,i) = tempData(:,order(i));
end

end