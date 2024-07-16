function [rf,par]=design_hardpulse4fidall(flip,sv,b1max_ref,dt)
%DESIGN_HARDPULSE4FIDALL  Design block RF excitation pulse
% [rf,par] = design_hardpulse4fidall(flip,sv,b1max_ref,dt)
%     flip   Flip angle                                    [deg] (90)
%b1max_ref   B1max reference of fidall (relative to 1H)    [T]   (0.1281d-4)
%       dt   Sampling dwell time                           [s]   (4d-6)
%
%       rf   RF waveform
%      par   Parameter structure
%
% See also Create5xxRfExtFile
% 7/2019  Rolf Schulte
if nargin<1, help(mfilename); return; end

if ~exist('flip','var'),     flip = []; end
if isempty(flip),            flip = 90; end
if ~exist('sv','var'),       sv = []; end
if isempty(sv),              sv = false; end
if ~exist('b1max_ref','var'),b1max_ref = []; end
if isempty(b1max_ref),       b1max_ref = 0.1281d-4; end
if ~exist('dt','var'),       dt = []; end
if isempty(dt),              dt = 4d-6; end


%% misc pars
nn = flip/180*pi/gyrogamma/b1max_ref/dt;    % #rf-pts for given peak B1
n = ceil((nn+1)/2)*2;                       % add 0 to end and round up
t = n*dt;                                   % pulse duration [s]
rf = [ones(1,n-1),0];                       % rf pulse shape
b1max = b1max_ref*nn/(n-1);                 % b1max of pulse [T]
bwcs = 1200*1d-3/t;                         % spectral bandwidth [Hz]
% pw_rf=1ms->bw_rf1=1200Hz (FWHM of abs(Mxy) in small tip angle) 


%% print pulse info
fprintf('flip=%g[deg]; t=%g[us]; n=%d; dt=%g[us]\n',flip,t*1d6,n,dt*1d6);
fprintf('bwcs=%g[Hz]; b1max=%g[uT]; b1max_ref=%g[uT]\n',...
    bwcs,b1max*1d6,b1max_ref*1d6);


%% exporting waveform
if sv
    fname = sprintf('hardpulse_flip%d_pw%d',round(flip),t*1d6);
    fprintf('Saving RF to ''%s''\n',fname);
    
    par.rftype = 'freq';
    par.pw = n*dt*1d6;
    par.isodelay = 0.5;
    par.bwcs = bwcs;
    par.maxB1 = b1max*1d4;
    par.minwidth = 1;
    par.flip = flip;
    
    Create5xxRfExtFile([fname '.rho'], rf, par);
else
    par = [];
end


end   % main function design_freqselrf.m

