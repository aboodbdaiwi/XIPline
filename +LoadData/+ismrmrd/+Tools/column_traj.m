function traj_out = column_traj(traj)

traj_x = reshape(traj(1,:,:),1,[])';
traj_y = reshape(traj(2,:,:),1,[])';
traj_z = reshape(traj(3,:,:),1,[])';

traj_out = [traj_x traj_y traj_z];