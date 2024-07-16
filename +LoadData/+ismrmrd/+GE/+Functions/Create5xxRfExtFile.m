function hdr = Create5xxRfExtFile(filename,rfamp,par)
%Create5xxRfExtFile
% Write RF wave form in 5xx format. Write header, data, and checksum.
% The 5xx header is 68 bytes long (including checksum). The checksum
% is calculated in little endian integer format. Total file size is
% 64 bytes header + 2*npts bytes data + 4 bytes checksum. If variable
% "par" is omitted then file is written in 5x format.
%
% hdr = Create5xxRfExtFile(filename, rfamp, par)
%
%       filename   file name for wave form in 5xx format
%          rfamp   wave form normalized to interval          (-1,1)
%     par.rftype   type of rf pulse
%                  'freq' frequency selective, e.g. hard pulse
%                  'spat' spatially selective, e.g. SLR
%                  'spsp' spectral-spatial
%         par.pw   nominal pulse width                       [us]
%   par.isodelay   relative pulse center reference from end  (0,1)
%                  0=perfectly self-refocused
%                  0.5=symmetric pulse
%                  1=max phase
%       par.bwcs   nominal bandwidth in spectral direction   [Hz]
%       par.bwsp   nominal bandwidth in spatial direction    [Hz]
%                  (only required for 'spsp')
%      par.maxB1   nominal maximum B1 for proton             [G]
%                  maxB1 reference in fidall = 0.128128
%      par.nlobe   (only required for 'spsp')
%     par.pwlobe   (only required for 'spsp')
%    par.pwlobea   (only required for 'spsp')
%    par.pwlobed   (only required for 'spsp')
%   par.minwidth   minimum slice thickness for proton        [mm]
%       par.flip   nominal flip angle                        [deg]
%
%   hdr.magicnum   magic number for 5x format ('IPG0')
%     hdr.hdrlen   number of bytes in file
%     hdr.revnum   revision number (= 1)
%    hdr.filelen   equal to hdrlen
%hdr.magicnum5xx   magic number for 5xx format
%                  ('HARD', 'SLR0', or 'SPSP')
%   hdr.isodelay   par.isodelay*par.pw  [us]
%      hdr.maxB1   par.maxB1*1e9
%   hdr.minwidth   par.minwidth*1e3
%   hdr.checksum   checksum
%
% See Get5xChecksum() in psd/wtools/psdrfprof/rfextfile.c and
% wg_bf_calc_chksum() in mgd/agb/pgen/wave/wg_bf_util.c
%
%   2012 Martin Janich
% 8/2020 Rolf Schulte
if (nargin<1), help(mfilename); return; end


%% parameters and checks
HDR_LEN = 68;
HDR_LEN_UNUSED = 48;
HDR_LEN_UNUSED_5xx = 8; % unused bytes in header
% SIZEOF_INT = 4;
SIZEOF_SHORT = 2;

if nargin<2, error('Creat5xxRfExtFile requires >=2 input variables'); end
if nargin<3, par = []; end

npts = length(rfamp);             % number of points
if (npts/2-floor(npts/2))~=0
    error('length(rfamp)(=%g) is odd; use even numbered waveform',npts);
end
if ~isempty(par)
    if abs(par.pw-npts*4)>1d-8
        warning('Pulse width (=%g) must be divisible HW cycle rate (=4) times resolution (=%g)',...
            par.pw,npts);
        input('continue?');
    end
    if par.maxB1>0.128128
        warning('par.maxB1(=%g[G])>0.128128[G] (fidall reference)',...
            par.maxB1);
    end
    if par.maxB1<0.128128/10
        warning('par.maxB1(=%g[G])<0.128128/10[G] (fidall reference)',...
            par.maxB1);
    end
    if par.pw>65535
        warning('par.pw(=%g)>65535',par.pw);
    end
end


%% check rfamp
if max(abs(rfamp(:)))>1
    warning('max(abs(rfamp(:)))(=%g)>1: normalising',max(abs(rfamp(:))));
    rfamp = rfamp/max(abs(rfamp(:)));
end
if max(abs(rfamp(:)))<0.9
    warning('max(abs(rfamp(:)))(=%g)<0.9',max(abs(rfamp(:)))); 
end
if abs(rfamp(1))>eps
    warning('rfamp(1)(=%g)~=0',rfamp(1));
end
if abs(rfamp(end))>eps
    warning('rfamp(end)(=%g)~=0',rfamp(end));
end


%% translate into instruction amplitudes + final bit=1
rfiamp = round(rfamp*32766/2)*2;
if abs(rfiamp(end))>eps 
    warning('rfiamp(end)(=%g)~=0; setting to 1',rfiamp(end)); 
end
if npts>0, rfiamp(end) = 1; end


%% actual file writing
% open file
fid = fopen(filename,'wb','b');   % write big endian

% magic number
hdr.magicnum = 'IPG0';
fwrite(fid,hdr.magicnum,'char*1');

% Set hdr length, to size of hdr section + checksum + data
hdr.hdrlen = HDR_LEN+npts*SIZEOF_SHORT+mod(npts,2)*SIZEOF_SHORT;
fwriteuint32(fid,hdr.hdrlen);

% revision number - number is hardcoded to a 1
hdr.revnum = 1;
fwriteuint32(fid,hdr.revnum);

% filelen = hdrlen
hdr.filelen = hdr.hdrlen;
fwriteuint32(fid,hdr.filelen);

if isempty(par)
    fprintf('Writing 5x file to ''%s''\n',filename);

    % unused section of the header packet
    for i=1:HDR_LEN_UNUSED
        fwrite(fid,0,'uint8');
    end

else
    % check rf pulse type
    if isfield(par,'rftype') && strcmp(par.rftype,'spat')
        hdr.magicnum5xx = 'SLR0'; % magic number for 5xx format
        par.bwsp = 0;
        par.nlobe = 0;
        par.pwlobe = 0;
        par.pwlobea = 0;
        par.pwlobed = 0;
    elseif isfield(par,'rftype') && strcmp(par.rftype,'freq')
        hdr.magicnum5xx = 'HARD'; % magic number for 5xx format
        par.bwsp = 0;
        par.nlobe = 0;
        par.pwlobe = 0;
        par.pwlobea = 0;
        par.pwlobed = 0;
    elseif isfield(par,'rftype') && strcmp(par.rftype,'spsp')
        hdr.magicnum5xx = 'SPSP'; % magic number for 5xx format
        % check mandatory entries in par structure
        if ~isfield(par,'bwsp')
            error('Creat5xxRfExtFile: par.bwsp is not defined.'); 
        end
        if ~isfield(par,'nlobe')
            error('Creat5xxRfExtFile: par.nlobe is not defined.'); 
        end
        if ~isfield(par,'pwlobe')
            error('Creat5xxRfExtFile: par.pwlobe is not defined.');
        end
        if ~isfield(par,'pwlobea')
            error('Creat5xxRfExtFile: par.pwlobea is not defined.');
        end
        if ~isfield(par,'pwlobed')
            error('Creat5xxRfExtFile: par.pwlobed is not defined.');
        end
        if abs(par.pw/length(rfamp)-4)>1d-12
            warning('pw/res=%g: gradient update time ~= 4[us]\n',...
                par.pw/length(rfamp));
        end
    else
        error('Creat5xxRfExtFile: par.rftype is not defined correctly.');
    end
       
    % check mandatory entries in par structure
    fprintf('Writing 5xx %s file %s\n',hdr.magicnum5xx,filename);
    if ~isfield(par,'pw')
        error('Creat5xxRfExtFile: par.pw is not defined.'); 
    end
    if ~isfield(par,'isodelay')
        error('Creat5xxRfExtFile: par.isodelay is not defined.'); 
    end
    if par.isodelay>1, warning('par.isodelay(=%g)>1',par.isodelay); end
    if par.isodelay<0, warning('par.isodelay(=%g)<0',par.isodelay); end
    hdr.isodelay = par.isodelay*par.pw+0.5;
    if ~isfield(par,'bwcs')
        error('Creat5xxRfExtFile: par.bwcs is not defined.'); 
    end
    if ~isfield(par,'maxB1')
        error('Creat5xxRfExtFile: par.maxB1 is not defined.'); 
    end
    hdr.maxB1 = par.maxB1*1e9;
    if ~isfield(par,'minwidth')
        error('Creat5xxRfExtFile: par.minwidth is not defined.'); 
    end
    hdr.minwidth = par.minwidth*1e3;
    if ~isfield(par,'flip')
        error('Creat5xxRfExtFile: par.flip is not defined.'); 
    end
    if par.bwcs<1,  error('par.bwcs(=%g)<1',  par.bwcs); end
    if par.bwcs<10, warning('par.bwcs(=%g)<1',par.bwcs); end
        
    % write to file
    fwrite(fid,hdr.magicnum5xx,'char*1');
    fwriteint32(fid,par.pw);
    fwriteint32(fid,npts);
    fwriteint32(fid,hdr.isodelay);
    fwriteint32(fid,par.bwcs);
    fwriteint32(fid,par.bwsp);
    fwriteint32(fid,hdr.maxB1);
    fwriteuint16(fid,par.nlobe);
    fwriteuint16(fid,par.pwlobe);
    fwriteuint16(fid,par.pwlobea);
    fwriteuint16(fid,par.pwlobed);
    fwriteuint16(fid,hdr.minwidth);
    fwriteuint16(fid,par.flip);
    
    % unused section of the header packet
    for i=1:HDR_LEN_UNUSED_5xx
        fwrite(fid,0,'uint8');
    end
    
end

% Write the amplitude data
fwriteint16(fid,rfiamp); % write big endian short integers

% if uneven length, write extra short
if mod(npts,2)>0
    fwriteint16(fid,0);
end

fclose(fid);

% Write the 5x checksum and close file
hdr.checksum = Get5xChecksum(filename);
fid = fopen(filename,'ab','b'); % append write big endian
fwriteint32(fid, hdr.checksum);
fclose(fid);

fprintf('Completed successfully with checksum =%d\n',hdr.checksum);

end   % main function


%% sub-functions
function fwriteuint32(fid,value)
if value > intmax('uint32') || value < intmin('uint32')
    error('fwriteuint32: value %d outside of range of unsigned 32 bit integers.',value);
end
fwrite(fid,value,'uint32');
end


function fwriteint32(fid,value)
if value > intmax('int32') || value < intmin('int32')
    error('fwriteint32: value %d outside of range of signed 32 bit integers.',value);
end
fwrite(fid,value,'int32');
end


function fwriteuint16(fid,value)
if value > intmax('uint16') || value < intmin('uint16')
    error('fwriteuint16: value %d outside of range of unsigned 16 bit integers.',value);
end
fwrite(fid,value,'uint16');
end


function fwriteint16(fid,value)
if max(value) > intmax('int16') || min(value) < intmin('int16')
    error('fwriteint16: value %d outside of range of signed 16 bit integers.',value);
end
fwrite(fid,value,'int16');
end


function checksum = Get5xChecksum(filename)
% Read big endian long integers and accumulate them. Return the bit
% complement.

% e.g. the checksum of an empty pulse file written with 
% Create5xxRfExtFile('temp.rho',[]) is -1229998010.

fid = fopen(filename,'rb','b'); % open in big endian mode
cs = fread(fid,'int32'); % read the whole file
fclose(fid);
checksum = 0;

for i=1:length(cs)
    checksum = checksum + cs(i);
    % force checksum to range of long integer
    if checksum > intmax('int32')
        % checksum = -bitcmp(checksum,32)-1;
        checksum = -bitcmp(checksum,'uint32')-1;
    end
    if checksum < intmin('int32')
        % checksum = bitcmp(-checksum,32)+1;
        checksum = bitcmp(-checksum,'uint32')+1;
    end
end

checksum = -checksum-1;

end

