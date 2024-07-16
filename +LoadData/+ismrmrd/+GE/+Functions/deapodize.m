function matrix = deapodize(matrix_dim,alpha)
%DEAPODIZE Deapodization function for grid_interp
%
% Created by Philip Beatty
% Modified by Rolf Schulte (changed scaling in last row)

Gx = matrix_dim(1);
Gy = matrix_dim(2);

if ~exist('alpha','var'), alpha = 1.5; end
S = 64;

Nx = 2 * floor(Gx / (2*alpha));
Ny = 2 * floor(Gy / (2*alpha));

gridFOV = [Gx Gy];

matrix = zeros(8);
matrix(5,5) = 1;

trajectory_x = (0:(2*S-1)) ./ (S * Gx);
trajectory_x = [trajectory_x; zeros(size(trajectory_x)); ones(size(trajectory_x))];

trajectory_y = (0:(2*S-1)) ./ (S * Gy);
trajectory_y = [zeros(size(trajectory_y)); trajectory_y; ones(size(trajectory_y))];

half_kernel = igrid_interp(matrix, trajectory_x, gridFOV);
deapod_FT = zeros(S*Gx, 1);
deapod_FT(1:length(half_kernel)) = half_kernel;
deapod_FT((length(deapod_FT)-length(half_kernel)+2):end) = flipud(half_kernel(2:end));
deapod = ifft(deapod_FT);
%deapod_x = 1./ real([flipud(deapod(2:(Nx/2+1))); deapod(1:(Nx/2))]);
deapod_x = (1./ real([flipud(deapod(2:(Nx/2+1))); deapod(1:(Nx/2))]) ...
    .* sub_sinc(((-Nx/2):(Nx/2-1)).' ./ (S * Gx)).^2);

half_kernel = igrid_interp(matrix, trajectory_y, gridFOV);
deapod_FT = zeros(S*Gy, 1);
deapod_FT(1:length(half_kernel)) = half_kernel;
deapod_FT((length(deapod_FT)-length(half_kernel)+2):end) = flipud(half_kernel(2:end));
deapod = ifft(deapod_FT);
%deapod_y = 1./ real([flipud(deapod(2:(Ny/2+1))); deapod(1:(Ny/2))]);
deapod_y = (1./ real([flipud(deapod(2:(Ny/2+1))); deapod(1:(Ny/2))]) ...
    .* sub_sinc(((-Ny/2):(Ny/2-1)).' ./ (S * Gy)).^2);

matrix = zeros(matrix_dim);

x_shift = Gx/2 - Nx/2;
y_shift = Gy/2 - Ny/2;

matrix((1:Nx)+x_shift,(1:Ny)+y_shift)= ...
    ...%repmat(deapod_x,[1 Ny]).*repmat(deapod_y',[Nx 1])./(S*4);
    repmat(deapod_x,[1 Ny]).*repmat(deapod_y',[Nx 1])/(Gx*Gy);

end      % main function deapodize



%% sub-functions
function y = sub_sinc(x)

i=find(x==0);                                                              
x(i)= 1;      % From LS: don't need this is /0 warning is off                           
y = sin(pi*x)./(pi*x);                                                     
y(i) = 1;   

end      % sub-function sub_sinc
