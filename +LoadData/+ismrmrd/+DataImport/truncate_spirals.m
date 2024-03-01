function [Gas_Fid,Gas_Traj] = truncate_spirals(FID,Traj,scale)

Traj2 = Traj*scale;

rad = squeeze(sqrt(Traj2(1,:,:).^2 + Traj2(2,:,:).^2 + Traj2(3,:,:).^2)); 

toobig = find(rad(:,1) >0.5);

FID(toobig,:) = [];
Traj2(:,toobig,:) = [];

Gas_Traj = Traj2;
Gas_Fid = FID;