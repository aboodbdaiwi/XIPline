function rf = read_rho(shape_name)
%READ_RHO  Read in pulse shape from GE waveform file.
%    rf = read_rho(shape_name)
%
%  11/2013  Rolf Schulte
if nargin<1, help(mfilename); return; end

if ~exist(shape_name,'file'), warning('read_rho:file','File not found'); end
rf = rp(shape_name);         % reading pulse shape


function pulse = rp(shape_name)
% fprintf('Reading pulse shape %s\n', shape_name);
fid = fopen(shape_name, 'rb','b');
%header = fread(fid,64,'*char')
% header is 64 bytes long
h1 = fread(fid,4,'*char');
h = fread(fid,3,'int32');
if h(1)~=h(3),
    warning('read_rho:header','HeaderByteCount~=FileByteCount');
    fprintf('\tHeaderByteCount = %g\n\tFileByteCount = %g\n',h(1),h(3));
end

% skip header
fseek(fid,64,'bof');

DataByteCount = h(1) - 64 - 4 - 4;
pulse = fread(fid,DataByteCount/2,'int16');
checksum = fread(fid,inf,'int16');
fclose(fid);

if ~any(length(checksum)==[3 4]), 
    warning('read_rho:checksum','Length of checksum ~= 3 or 4'); 
end

% fprintf('header:\n%s\n',header);
pulse = pulse(:).';
