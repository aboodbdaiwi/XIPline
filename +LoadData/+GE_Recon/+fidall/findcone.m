%function [g,k,nint] = findcone(RES,FOV,LEN,THETA,PRECISION,TS,OS,SMAX,GMAX,DCF,MINDENS);
%
%   Generates a gradient corresponding to the cone with the given resolution, 
%  field of view, length, and range of polar angles.
%
% g - Cone Gradient
% k - Cone k-space trajectory 
% nint - Number of interleaves (only roughly for a range of polar angles)
% RES - Desired Resolution in mm
% FOV - Desired Field of View in cm
% LEN - Desired Length
% THETA - Range of polar angles (e.g. [0 pi/2] for 0 to 90 degrees)
% PRECISION - Precision of number of interleaves 
% TS -  Sampling period in seconds
% OS -  Oversampling Rate 
% SMAX - Slew rate (in G/cm/s)
% GMAX - Max Gradient Amplitude (in G)
% DCF - Density Compensation Factor
% MINDENS - Minimum Density

function [g,k,nint,lenro] = findcone(RES,FOV,LEN,THETA,PRECISION,TS,OS,SMAX,GMAX,DCF,MINDENS,APOD,REW);
if (nargin<5)
	PRECISION = 0.1;
end
if (nargin<6)
        TS = 0.000004;
end
if (nargin<7)
	OS = 4;
end
if (nargin<8)
        SMAX = 15000;
end
if (nargin<9)
        GMAX = 3.98;
end
if (nargin<10)
	DCF = 0.0;
end
if (nargin<11)
	MINDENS = 0.0;
end
if (nargin<12)
	APOD = [0 1];
end
if (nargin<13)
	REW = 0;
end

NINTlo = min(PRECISION,0.1);
NINThi = 100;
LENlo = 10000000;

curNINT = NINThi;
MAXLEN = LEN+5;

% Try a huge number of interleaves to make sure it's possible to generate
% such a waveform.
[g,k,lenro] = gencone(RES,FOV,1e8,THETA,MAXLEN,TS,SMAX,GMAX,OS,DCF,MINDENS,APOD,REW);
if (length(g)>LEN)
	disp('Sorry, but it is impossible to achieve that resolution in that length of time');
	g = [];
	gr = [];
	nint = 0;
	return
end

[g,k,lenro] = gencone(RES,FOV,curNINT,THETA,MAXLEN,TS,SMAX,GMAX,OS,DCF,MINDENS,APOD,REW);
leng = length(g); 

LENhi = leng;
% Keep doubling the number of interleaves until we get a waveform which
% satisfies the requirements. 
while ((leng>LEN))
	NINTlo = NINThi;
	NINThi = NINThi*2;
	curNINT = NINThi;
	LENhi = leng;
        [g,k,lenro] = gencone(RES,FOV,curNINT,THETA,MAXLEN,TS,SMAX,GMAX,OS,DCF,MINDENS,APOD,REW);
	leng = length(g); 
end

% Do a binary search over the remaining range until we achieve 
% the requested precision.
while ((((NINThi-NINTlo)>PRECISION) || (leng>LEN)) && (LENhi~=LEN))
	curNINT = (NINThi-NINTlo)/2+NINTlo;
        [g,k,lenro] = gencone(RES,FOV,curNINT,THETA,MAXLEN,TS,SMAX,GMAX,OS,DCF,MINDENS,APOD,REW);
	leng = length(g); 
	if (leng>LEN)
             NINTlo = curNINT;
	     LENlo = leng;
	else
             NINThi = curNINT;
	     LENhi = leng;
	end
	%[NINTlo NINThi]
        %[leng lenro]
end

  nint = NINThi;
  [g,k,lenro] = gencone(RES,FOV,nint,THETA,MAXLEN,TS,SMAX,GMAX,OS,DCF,MINDENS,APOD,REW);

