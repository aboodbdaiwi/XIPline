function [DiceScore] = DiceCoeff(Volume1, Volume2)
%   Author: Abdullah S. Bdaiwi
%   Work email: abdullah.bdaiwi@cchmc.org
%   Personal email: abdaiwi89@gmail.com
%   Website: https://www.cincinnatichildrens.org/research/divisions/c/cpir

% Making variables to calculate dice coeff of exp/insp
    Volume1_Vec = Volume1(:);       % Making matrix into vector
    Volume1_Tot = sum(Volume1_Vec); % Summing total # of voxels in Volume 1
    
    Volume2_Vec = Volume2(:);       % Making matrix into vector
    Volume2_Tot = sum(Volume2_Vec); % Summing total # of voxels in Volume 2
    
    CommonVoxels = sum(Volume1_Vec & Volume2_Vec); % Common voxel locations between volumes
    
% Calculating DSC and writing output
    DiceScore.DSC = (2 * CommonVoxels) / (Volume1_Tot + Volume2_Tot);
    DiceScore.CommonVoxels = CommonVoxels;
    DiceScore.TotalVolume1 = Volume1_Tot;
    DiceScore.TotalVolume2 = Volume2_Tot;
    DiceScore.Volume1 = Volume1;
    DiceScore.Volume2 = Volume2;

%     save(savefile, 'DiceScore');
%     cd(start);
end

