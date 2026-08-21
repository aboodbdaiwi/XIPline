function matrix=grid_interp(data,trajectory,FOV,matrix_dim)
%GRID_INTERP  Interface to C-based gridding routine
% matrix = grid_interp(data, trajectory, FOV, matrix_dim)
% 
%       data: k-space data acquired along the arbitrary trajectory
% trajectory: arbitrary k-space trajectory (normalised to +-0.5)
%             3 x points: 1=kx;2=ky;3=dcf (density compensation function)
%        FOV: Scale k-space to physical field of view [x y] (doubles)
% matrix_dim: Cartesian resolution [x y] (ints)
%
%     matrix: gridded data (not yet FFTed nor deapodised)
%
%Literature: itomi24_799_2005_beatty
%Compilation:
%>> mex grid_interp_mex.c -O
%
% Created   2/2007  Philip Beatty
% Modified  11/2020  Rolf Schulte
if (nargin<1), help(mfilename); return; end

% do error checking here
if isreal(data)
    data = complex(data);
end

% run mex file
matrix = grid_interp_mex(data, trajectory, FOV, matrix_dim);

end      % main function grid_interp.m
