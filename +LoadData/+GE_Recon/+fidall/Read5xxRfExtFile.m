function [rf,par,hdr_section,amp] = Read5xxRfExtFile(filename)
%READ5XXRFEXTFILE  Read in waveforms from external rf file
% Read RF wave form in 5xx format: header, data, and checksum.
% The 5xx header is 68 bytes long (including checksum). The checksum
% is calculated in little endian integer format. Total file size is
% 64 bytes header + 2*npts bytes data + 4 bytes checksum. If variable
% "par" is omitted then file is written in 5x format.
%
%[rf,par,hdr_section,amp] = Read5xxRfExtFile(filename)
%
% filename                : file name for wave form in 5xx format
% rfamp                   : wave form normalized to interval [-1,1]
% par.pw                  : nominal pulse width [us]
% par.isodelay            : pulse center reference relative to pw
% par.bwcs                : nominal bandwidth in spectral direction
% par.bwsp                : nominal bandwidth in spatial direction
% par.maxB1               : nominal maximum B1 for proton [G]
% par.nlobe
% par.pwlobe
% par.pwlobea
% par.pwlobed
% par.minwidth            : minimum slice thickness for proton [mm]
% par.flip                : nominal flip angle [deg]
% hdr_section.magicnum    : magic number for 5x format ('IPG0')
% hdr_section.hdrlen      : number of bytes in file
% hdr_section.revnum      : revision number (= 1)
% hdr_section.filelen     : equal to hdrlen
% hdr_section.magicnum5xx : magic number for 5xx format
%                          ('HARD', 'SLR0', or 'SPSP')
% hdr_section.isodelay    : par.isodelay*par.pw
% hdr_section.maxB1       : par.maxB1*1e9
% hdr_section.checksum    : checksum
%
% See Get5xChecksum() in psd/wtools/psdrfprof/rfextfile.c and
% wg_bf_calc_chksum() in mgd/agb/pgen/wave/wg_bf_util.c
%
% Martin Janich 5 Oct 2012
% 6/2020 Rolf Schulte
if nargin<1, help(mfilename); return; end

% HDR_LEN = 68;
% HDR_LEN_UNUSED = 48;
HDR_LEN_UNUSED_5xx = 8; % unused bytes in header
% SIZEOF_INT = 4;
% SIZEOF_SHORT = 2;


%% open file
fid = fopen(filename,'rb','b'); % read, big endian
hdr_section.magicnum = fread(fid,4,'*char*1').';
if hdr_section.magicnum~='IPG0'
    warning('hdr_section.magicnum(=%c)~=''IPG0''',hdr_section.magicnum);
end


%% read hdr length (size of hdr section + checksum + data)
hdr_section.hdrlen = fread(fid,1,'int32');
hdr_section.revnum = fread(fid,1,'int32');       % revision number = 1
if hdr_section.revnum~=1
    warning('hdr_section.revnum(=%g)~=1',hdr_section.revnum);
end
hdr_section.filelen = fread(fid,1,'int32');      % filelen = hdrlen
hdr_section.magicnum5xx = fread(fid,4,'*char*1').'; % magic number for 5xx format
switch hdr_section.magicnum5xx
    case 'SLR0', fprintf('Reading 5xx SLR0 file\n');
    case 'HARD', fprintf('Reading 5xx HARD file\n');
    case 'SPSP', fprintf('Reading 5xx SPSP file\n');
    otherwise
        fprintf('Reading standard 5 file\n');
        fprintf('Attention: limited header\n');
        fprintf('Use read_rho.m instead\n');
end


%% read par header
par.pw = fread(fid,1,'int32');
test = fread(fid,1,'int32');  % resolution; if par was not empty this should be 0;
hdr_section.isodelay = fread(fid,1,'int32'); % isodelay
par.bwcs=fread(fid,1,'int32');              % bwcs
par.bwsp=fread(fid,1,'int32');              % bwsp
position_maxB1 = ftell(fid);
hdr_section.maxB1 = fread(fid,1,'int32');
par.maxB1 = hdr_section.maxB1/1e9;          % max_b1
par.nlobe =   fread(fid,1,'int16');         % nlobe
par.pwlobe =  fread(fid,1,'int16');         % pwlobe
par.pwlobea = fread(fid,1,'int16');         % pwlobea
par.pwlobed = fread(fid,1,'int16');         % pwlobed
par.minwidth = fread(fid,1,'int16');        % minwidth
par.flip = fread(fid,1,'int16');            % nom_fa
test2 = fread(fid,HDR_LEN_UNUSED_5xx,'uint8'); % unused Bytes; should be zeros
if test2~=0, warning('test2(=%g)~=0',test2); end


%% read rf amplitude data
temp = fread(fid,'int16');
amp = temp(1:end-2,1);       % discard last 2 points
rf_points = length(amp);
if rf_points/2-floor(rf_points/2) ~= 0
    warning('odd number of waveform points (=%g)',rf_points);
end


%% reading checksum
fseek(fid, -4, 'eof');
hdr_section.checksum = fread(fid,1,'int32');
fclose(fid);                 % closing file


%% checks + scaling
if rf_points~=par.pw/4
    warning('wrong # of rf-points');
end

if par.maxB1>0
    rf = amp/32766*par.maxB1*1d-4;
else
    warning('par.maxB1(=%g)<=0; omitting scaling',par.maxB1)
    rf = amp;
end

end   % main function Read5xxRfExtFile.m
