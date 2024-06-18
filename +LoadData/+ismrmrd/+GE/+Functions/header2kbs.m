function KBS = header2kbs(h)
%HEADER2KBS  Convert fidall header information into kbs scaling factor for 
% Bloch-Siegert B1+ mapping
% kbs = header2kbs(h)
%    h  header
%       h.rdb_hdr.user30 = Bloch-Siegert pulse frequency   [Hz]
%       h.rdb_hdr.user31 = Bloch-Siegert pulse width       [us]
%  kbs  scaling factor for phase difference                [rad]
%       b1map = sqrt(b1phase/2/KBS);
%
% 1/2016 Rolf Schulte
if (nargin<1), help(mfilename); return; end;


nucleus = h.rdb_hdr.user2;        % nucleus

%% extract Bloch-Siegert parameters from header
off_freq = h.rdb_hdr.user30;      % Bloch-Siegert pulse frequency [Hz]
off_pw = h.rdb_hdr.user31;        % pw blosi pulse [us]
KBS = 4109.6;      % deg/G^2; conversion for Bloch_Siegert.rho (8ms,4kHz)


%% phase-shift constant
% KBS -> deg/G^2 Bloch-Siegert shift conversion from phase to B1 magnitude
KBS = KBS/100^2;                  % convert Gauss into uT
KBS = KBS/(off_freq/4000);        % scale for frequencies other than 4kHz
KBS = KBS*(off_pw/8000);          % scale to pw's other than 8ms
scale = (gyrogamma(nucleus)/gyrogamma('1h'))^2;
KBS = scale*KBS;                  % scale to nucleus
KBS = KBS*pi/180;                 % convert from deg to rad
