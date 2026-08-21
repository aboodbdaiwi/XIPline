function [output] = readh5complex(filename)
%READH5COMPLEX Reads h5 file with \real and \imag datasets and returns a
% complex matrix with single precision
% output = readh5complex(filename)
%
% 12/2024 Kylie Yeung
if (nargin<1), help(mfilename); return; end

if isempty(regexpi(filename,'\.h5$', 'once')), filename = [filename '.h5']; end
if ~exist(filename,'file'), warning('file not found: %s',filename); end

output_real = h5read(filename, '/real');
output_imag = h5read(filename, '/imag');

if ~isa(output_real,'single'), warning('output_real not single'); end
if ~isa(output_imag,'single'), warning('output_imag not single'); end

output = single(output_real + 1i*output_imag);

end      % readh5complex.m
