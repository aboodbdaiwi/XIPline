function sortedTraj = SortUTETraj_AcqOrder2(unsortedTraj, DataSizes, order)

tempTraj   = zeros(size(unsortedTraj), 'like', unsortedTraj);
sortedTraj = zeros(size(unsortedTraj), 'like', unsortedTraj);

PtsPerKz = DataSizes(2);
NKz      = DataSizes(3);
NPro     = PtsPerKz * NKz;
Pts      = 1:PtsPerKz:NPro;
NewPts   = 1:NKz;

for i = 1:PtsPerKz
    tempTraj(:,:, NewPts + NKz*(i-1)) = unsortedTraj(:,:, Pts + (i-1));
end

for i = 1:NPro
    sortedTraj(:,:, order(i)) = tempTraj(:,:, i);
end

end