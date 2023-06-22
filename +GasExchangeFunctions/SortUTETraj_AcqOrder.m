function sortedTraj = SortUTETraj_AcqOrder(unsortedTraj, DataSizes, order)
%SortUTEData_AcqOrder - sorts data using projection order
%  assumes unsortedTraj is [read out, proj]
%  assumes DataSize is [read out, profs, kz(interleafs/hubs)]

%initialize output
tempTraj = zeros(size(unsortedTraj));
sortedTraj = zeros(size(unsortedTraj));

%Calculate params
PtsPerKz = DataSizes(2);
NKz = DataSizes(3);
NPro = PtsPerKz * NKz;
Pts = 1:PtsPerKz:NPro;
NewPts = 1:NKz;

%Sort kz
for i = 1:PtsPerKz
    tempTraj(:,:,NewPts + NKz*(i-1)) = unsortedTraj(:,:,Pts+(i-1));
end

%Final Sort
for i = 1:NPro
    sortedTraj(:,:,i) = tempTraj(:,:,order(i));
end

end