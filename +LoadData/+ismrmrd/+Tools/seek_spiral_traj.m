function traj = seek_spiral_traj(twix_obj)

%% Function to create spiral trajectories given a twix object.
%This will look in the parent path of the package for the proper trajectory
%file. If it doesn't find the right measured trajectory file, it will
%calculate theoretical trajectories (which will probably be bad)

%Trajectories for a given protocol should be stored in parent_path/Traj -
%It will work best if each trajectory file name is sensibly named as 
%(name of protocol the traj are for)_traj, or something similar.
%Get Parent path
parent_path = which('Tools.seek_spiral_traj');
idcs = strfind(parent_path,filesep);%determine location of file separators
parent_path = parent_path(1:idcs(end-1)-1);%remove file

%% Identify protocol name:
prot_name = twix_obj.hdr.Config.ProtocolName;
%Replace spaces with underscores:
prot_name = strrep(prot_name,' ','_');
%% Compare protocol name to trajectory names - easiest and best way to find trajectory
all_traj = dir(fullfile(parent_path,'Traj'));

%Make sure these are all .dat files:
del_ind = [];
for i = 1:length(all_traj)
    if ~contains(all_traj(i).name,'.dat')
        del_ind = [del_ind i];
    end
end
all_traj(del_ind) = [];
proper_traj = nan;
if ~isempty(all_traj)
    %First test, compare names of trajectory files to protocol names:
    for i = 1:length(all_traj)
        if contains(all_traj(i).name,prot_name)
            proper_traj = i;
            break;
        end
    end

    %% Compare parameters in twix header - probably not great:
    if isnan(proper_traj)
        for i = 1:length(all_traj)
            traj_twix = mapVBVD(fullfile(parent_path,'Traj',all_traj(i).name));
            if length(traj_twix)>1
                traj_twix = traj_twix{2};
            end
            if Tools.compare_twix(twix_obj,traj_twix)
                proper_traj = i;
                break;
            end
        end
    end
end
%% Read in Trajecories
%If we managed to find the proper trajectories, the variable proper_traj is
%no longer nan, so we can read in trajectories. Otherwise, we'll have to
%calculate theoretical
if ~isnan(proper_traj)
    traj_twix = mapVBVD(fullfile(parent_path,'Traj',all_traj(proper_traj).name));
    if length(traj_twix)>1
        traj_twix = traj_twix{2};
    end
    traj = Tools.spiral_coords_from_dat(traj_twix,twix_obj);
else
    %Hardcode to no delay, although some may be needed - Let's try 
    %Delay = [3 3 3];
    Delay = [12 12 12];
    disp('No Trajectory File Found. Using Theoretical');
    traj = Tools.siemensspiralcoords(twix_obj,Delay);
end




