function saveh5complex(filename,input)
%SAVEH5COMPLEX Saves complex input to a h5 file with \real and \imag 
%  datasets with single precision
% saveh5complex(filename,input)
%
% 12/2024 Kylie Yeung
if (nargin<1), help(mfilename); return; end

if isempty(regexpi(filename,'\.h5$', 'once')), filename = [filename '.h5']; end
if exist(filename,'file'), warning('file existing: %s',filename); end


h5create(filename,'/real',size(input),'Datatype','single')
h5create(filename,'/imag',size(input),'Datatype','single')

h5write(filename,'/real',real(input))
h5write(filename,'/imag',imag(input))

end      % saveh5complex.m
