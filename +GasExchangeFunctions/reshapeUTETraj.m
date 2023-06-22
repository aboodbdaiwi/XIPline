function reshapedTraj = reshapeUTETraj(inputTraj)
%reshapeUTEData - reshapes UTE trajectories
%from [[3, read out, profs, kz(interleafs/hubs)] to [3, read out, proj]

%get traj sizes
nRO = size(inputTraj,2);
nProfs = size(inputTraj,3);
nKz = size(inputTraj,4);

%create zeros
reshapedTraj = zeros(3,nRO,nProfs*nKz);

%fill data
reshapedTraj(1,:,:) = reshape(inputTraj(1,:,:,:),nRO,[]);
reshapedTraj(2,:,:) = reshape(inputTraj(2,:,:,:),nRO,[]);
reshapedTraj(3,:,:) = reshape(inputTraj(3,:,:,:),nRO,[]);

end




