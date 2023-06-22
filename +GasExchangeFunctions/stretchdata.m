function stretchdata = stretchdata(data)

%Function to take oscilatory data of varying amplitude and make it all the
%same amplitude
%[HighPeaks,HighLoc] = findpeaks(data);
%[LowPeaks,LowLoc] = findpeaks(-data);

zerofind = GasExchangeFunctions.zero_cross(data);

stretchdata = data;

for i = 2:length(zerofind)
    low = zerofind(i-1);
    high = zerofind(i);
    dir1 = max(data(low:high));
    dir2 = min(data(low:high));
    if abs(dir2) < abs(dir1)
        Peak = 1;
        stretchdata(low:high) = data(low:high)*Peak/max(data(low:high));
    elseif abs(dir2) > abs(dir1)
        Peak = 1;
        stretchdata(low:high) = data(low:high)*Peak/abs(min(data(low:high)));
    end
end
