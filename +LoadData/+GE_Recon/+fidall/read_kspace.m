function [k,dcf] = read_kspace(fname,format,mafo)
%READ_KSPACE  Read in kspace data from bin or ascii file with format
% kx1,ky1,dcf1,kx2,ky2,dcf2, ... (like spiral.e output /w/config/kspfile)
%
% [k,dcf] = read_kspace(fname,format)
%   fname   File name for k-space file                          ('kspfile')
%  format   'bin' = float32; 'asc' = ascii                          ('bin')
%    mafo   Machine format                                            ('l')
%           'b' big-endian 
%           'l' little-endian (default, scanner output)
%       k    k-space points, kx+i*ky.
%     dcf    Density compensation factors.
%
% 11/2023 Rolf Schulte
if (nargin<1) && (nargout<1), help(mfilename); return; end


%% iput args
if ~exist('fname','var'),  fname = ''; end
if isempty(fname),         fname = 'kspfile'; end
if ~exist('format','var'), format = ''; end
if isempty(format),        format = 'bin'; end
if ~exist('mafo','var'),   mafo = ''; end
if isempty(mafo),          mafo = 'l'; end


%% read in data
switch lower(format)
    case 'bin'
        fid = fopen(fname,'r',mafo);
        [arr,n] = fread(fid,inf,'float32');
    case 'asc'
        fid = fopen(fname,'r');
        [arr,n] = fscanf(fid,'%f',inf);
end
fclose(fid);
l = floor(n/3);
arr = arr(1:3*l);

if (l < 1)
    warning('No data in file');
else
	fprintf('%d k-space samples in k-space file %s\n',l,fname);
end


%% assign to k and dcf
arr = reshape(arr,3,l);
k = arr(1,:) + 1i*arr(2,:);
dcf = arr(3,:);


end      % read_kspace.m
