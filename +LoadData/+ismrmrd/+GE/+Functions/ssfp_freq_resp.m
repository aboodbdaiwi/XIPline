function [ww,srf] = ssfp_freq_resp(flip,tr,t1,t2,df0,npc)
%SSFP_FREQ_RESP  Caculate frequency response of balanced SSFP train
% [ww,srf] = ssfp_freq_resp(flip,tr,t1,t2,hz,npc)
%     flip   Excitation flip angle of SSFP train           [deg]
%       tr   Repetition time                               [s]
%       t1   T1 relaxation time                            [s]
%       t2   T2 relaxation time                            [s]
%      df0   Off-resonance frequencies                     [Hz]
%      npc   #phase cycles
%       ww   Complex weights for combining data
%            (reshaped to same arbitrary size as df0)
%      srf   Spectral response function
%
% Literature: jmri38_2_2013_bieri.pdf -> Eq. 4+5
% 12/2020 Rolf Schulte
if (nargin<1), help(mfilename); return; end


%% default parameters + checks
if ~exist('npc','var'), npc = []; end
if isempty(npc),        npc = 1; end
nn = size(df0);
hz = reshape(df0,[prod(nn) 1]);


%% calculate spatial response function
e1 = exp(-tr/t1);       % Eq. 5 wrong, see mrm18_244_1991_mizumoto.pdf
e2 = exp(-tr/t2);
c = e2*(e1-1)*(1+cosd(flip));
d = (1-e1*cosd(flip))-(e1-cosd(flip))*e2.^2;
phi = 2*pi*hz*tr;
if npc>1
    phi = bsxfun(@plus,phi,2*pi*(0:npc-1)/npc);
end
srf = ((1-e1)*(1-e2*exp(-1i*phi))*sind(flip))./(c*cos(phi)+d);


%% normalise to get weights
% ww = conj(srf)./repmat(sum(abs(srf),2),[1 npc]);
ww = bsxfun(@rdivide,conj(srf),sum(abs(srf),2));
% ww = conj(srf);


%% reshape to original matrix size of df0 + #phase cycles
ww = reshape(ww,[nn npc]);


%%
if (nargout==1) && (npc==1)
    warning('function return weight as first output');
end


end      % ssfp_freq_resp.m
