function data = igrid_interp(matrix, trajectory, FOV)
%IGRID_INTERP  Interface to C-based gridding routine
% data = igrid_interp(matrix, trajectory, FOV)
% 
%     matrix: Cartesian k-space data (not yet FFTed nor deapodised)
% trajectory: arbitrary k-space trajectory (normalised to +-0.5)
%             3 x points: 1=kx;2=ky;3=dcf (density compensation function)
%        FOV: [x y] (doubles)
%
%       data: k-space data gridded to the arbitrary trajectory
%
%Literature: itomi24_799_2005_beatty
%
%Compilation:
%>> mex igrid_interp_mex.c -O
%
% Created   2/2007  Philip Beatty
% Modified  2/2007  Rolf Schulte
if (nargin<1), help(mfilename); return; end;

% do error checking here
if isreal(matrix),
    matrix = complex(matrix);
end

% run mex file
data = igrid_interp_mex(matrix, trajectory, FOV);
