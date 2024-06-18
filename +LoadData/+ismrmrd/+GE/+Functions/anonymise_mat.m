function anonymise_mat(filename_mat4anon)
%ANONYMISE_MAT Remove patient data from header structure in matlab mat file
% anonymise_mat(fname)
%  fname   Matlab mat file containing field h 
%
%  02/2021 Rolf Schulte
if (nargin<1), help(mfilename); return; end


%% input
if ~isempty(regexpi(filename_mat4anon,'\.mat$'))
    filename_mat4anon = filename_mat4anon(1:end-4); 
end
if ~exist([filename_mat4anon '.mat'],'file')
    error('fname (=''%s'') not found',[filename_mat4anon '.mat']);
end


%% load mat file
load([filename_mat4anon '.mat']);


%% anonymise h
if ~exist('h','var')
    error('header structure h not found'); 
end
h = header_anonymise(h);


%% save anoymised file
save([filename_mat4anon '_anon.mat']);


end      % anonymise_mat.m 
