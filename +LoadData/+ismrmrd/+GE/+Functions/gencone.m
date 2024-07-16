function [g,k,rolen] = gencone(RES,FOV,NINT,THETA,MAXLEN,TS,SMAX,GMAX,OVERSAMPLE,DCF,MINDENS,APOD,REW);
%function [g,k] = gencone(RES,FOV,NINT,THETA,MAXLEN,TS,SMAX,GMAX,OVERSAMPLE,DCF,MINDENS);
%
%  This function generates a particular cone, with a given resolution and field of view. 
%
% g gives the desired gradient
% k gives the desired k-space trajectory
% RES - Desired Resolution in mm 
% FOV - Desired Field of View in cm
% NINT - Desired Number of Interleaves (rough for range of polar angles)
% THETA - Range of polar angles (e.g. [0 pi/2] for 0 to 90 degrees)
% MAXLEN - Maximum length output allowed
% TS -  Sampling period in seconds
% SMAX - Slew rate (in G/cm/s)
% GMAX - Max Gradient Amplitude (in G)
% OVERSAMPLE - Amount to oversample 
% DCF - Density Compensation Factor
% MINDENS - Minimum Density

if (nargin<5)
	MAXLEN = 10000;
end
if (nargin<6)
	TS = 0.000004;
end
if (nargin<7)
	SMAX = 15000;
end
if (nargin<8)
	GMAX = 3.98;
end
if (nargin<9)
	OVERSAMPLE = 4;
end
if (nargin<10)
	DCF = 1.0;
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
if (length(FOV)<2)
	FOV(2) = FOV(1);
end
if (length(RES)<2)
	RES(2) = RES(1);
end

str = sprintf('gencone([%g %g],[%g %g],%g,[%g %g],%g,%g,%g,%g,%g,%g,%g,[%g %g]);',...
    RES(1),RES(2),FOV(1),FOV(2),NINT,min(THETA),max(THETA),MAXLEN,TS,...
    SMAX,GMAX,OVERSAMPLE,DCF,MINDENS,APOD(1),APOD(2));
%disp(str);

% Make sure theta is between 1e-5 and PI/2-1e-5 
% (the slight offsets mean we don't need to worry about divide by zero 
% errors everywhere)
THETA = min(THETA,pi/2-1e-5);
THETA = max(THETA,1e-5);

maxtheta = max(THETA);
mintheta = min(THETA);

kmax_xy = 5/RES(1);
kmax_z = 5/RES(2);

% Choose an appropriate angle in between the two thetas
THETArange = atan2(kmax_z*sin(maxtheta),kmax_xy*cos(mintheta));
% Find the kmax at that angle
KMAXrange = sub_intlineellipse(kmax_xy,kmax_z,THETArange);
% Determine the parametric theta for the chosen angle
tparam = atan2((KMAXrange*sin(THETArange)*kmax_xy),...
    (KMAXrange*cos(THETArange)*kmax_z));

% Determine the amount by which the waveforms will be scaled to cover the 
% range of thetas
xyscale = cos(tparam)/cos(mintheta);
zscale = sin(tparam)/sin(maxtheta);

% The circumferential FOV is always just FOV_xy
FOVcirc = FOV(1);

% The radial FOV is the FOV in the radial direction
truemaxtheta = angle(cos(maxtheta)*kmax_xy+1i*sin(maxtheta)*kmax_z);
truemintheta = angle(cos(mintheta)*kmax_xy+1i*sin(mintheta)*kmax_z);

FOVrad_maxtheta = 1/sub_intlineellipse(1/FOV(1),1/FOV(2),truemaxtheta);
FOVrad_mintheta = 1/sub_intlineellipse(1/FOV(1),1/FOV(2),truemintheta);

FOVrad = max(FOVrad_maxtheta,FOVrad_mintheta);

% Figure out the dtheta relationship (intracones) to adjust the density
curFOV = (FOV(1)*FOV(2)/sqrt((FOV(2)^2-FOV(1)^2)*sin(tparam)^2+FOV(1)^2));
curKMAX = (kmax_z*kmax_xy/sqrt((kmax_xy^2-kmax_z^2)*sin(tparam)^2+kmax_z^2));
dthetaCUR = 1/curFOV/curKMAX;
dthetaBASE = 1/FOV(2)/kmax_xy;
densadjust = dthetaCUR/dthetaBASE;

TS = TS/OVERSAMPLE;

str = sprintf('wcc(%g,[%g %g],%g,%g,%g,%g,%g,[%g %g],[%g %g],%g,[%g %g]);\n',...
    pi/2-THETArange,FOVcirc,FOVrad,DCF,KMAXrange,NINT,MAXLEN*OVERSAMPLE,TS,...
    SMAX*TS*xyscale, SMAX*TS*zscale,GMAX*xyscale,GMAX*zscale,...
    MINDENS*densadjust,APOD(1),APOD(2));

%disp(str);

[g,k,len] = wcc(pi/2-THETArange,[FOVcirc FOVrad],DCF,KMAXrange,NINT,...
    MAXLEN*OVERSAMPLE,TS,[SMAX*TS*xyscale SMAX*TS*zscale],...
    [GMAX*xyscale GMAX*zscale],MINDENS*densadjust,APOD);
g(1,1) = 0;
g(1,2) = 0;
g(1,3) = 0;

% Scale the waveforms up so that they can just be multiplied by the simple 
% scaling factor later
if (OVERSAMPLE>1)
	k = [k(1:OVERSAMPLE:len,1:2) k(1:OVERSAMPLE:len,3)];
	g = diff([0 0 0; k]);
else
	g = g(1:len,:);
end

g = [1/xyscale*g(:,1:2) 1/zscale*g(:,3)];
gc = (g(:,1)+1i*g(:,2));
gz = (g(:,3));
rolen = length(g);
if (REW>0)
	if (length(gc)>0) 
              gr = mrewind(sum(gc)*4258*TS,gc(end),GMAX,SMAX,TS);
        else 
	      gr = [];
        end
	if (length(gz)>0)
              grz = prewind(sum(gz)*4258*TS,gz(end),GMAX,SMAX,TS);
        else
	      grz = [];
        end

size(g);
size(real(gr));
size(grz);
if (length(grz)>length(gr))
	gr(length(grz)) = 0;
elseif (length(grz)<length(gr))
	grz(length(gr)) = 0;

end
size(real(gr));
size(grz);
g = [g; real(gr(2:end))' imag(gr(2:end))' grz(2:end)'];
end

k = cumsum(g)*4258*TS;

end      % main function gencone.m


%% sub-functions
% Finds the radius of a line segment starting at the origin
% and at an angle phi, and an ellipse with x-axis a, and y-axis b
%function r = intlineellipse(a,b,phi);
function r = sub_intlineellipse(a,b,phi)

r = b*a/sqrt((a^2-b^2)*sin(phi)^2+b^2);
end


