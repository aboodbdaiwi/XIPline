function reshapedData = reshapeUTEData(inputData)
%reshapeUTEData - reshapes UTE data from [[read out, profs,
%kz(interleafs/hubs)] to [read out, proj]

%get data sizes
nRO = size(inputData,1);
nProfs = size(inputData,2);
nKz = size(inputData,3);

%create zeros
reshapedData = zeros(nRO,nProfs*nKz);

%fill data
reshapedData(:,:) = reshape(inputData(:,:,:),nRO,[]);

end