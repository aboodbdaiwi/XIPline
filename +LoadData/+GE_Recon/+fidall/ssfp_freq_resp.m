function [srf,ww] = ssfp_freq_resp(flip,tr,T1,T2,df0,npc,docen)
%SSFP_FREQ_RESP  Caculate frequency response of balanced SSFP train
% [srf,ww] = ssfp_freq_resp(flip,tr,t1,t2,hz,npc)
%     flip   Excitation flip angle of SSFP train           [deg]
%       tr   Repetition time                               [s]
%       T1   T1 relaxation time                            [s]
%       T2   T2 relaxation time                            [s]
%      df0   Off-resonance frequencies                     [Hz]
%      npc   #phase cycles                                              (1)
%    docen   Start with centred band (0-180-0-...)                   (true)
%
%      srf   Spectral response function [1,npc]
%       ww   Complex weights for combining data
%            (reshaped to same arbitrary size as df0)
%
% Literature: jmri38_2_2013_bieri.pdf -> Eq. 4+5
%  4/2024 Rolf Schulte
if (nargin<1), help(mfilename); return; end

%% input
if ~exist('npc','var'), npc = []; end
if isempty(npc),        npc = 1; end
if ~exist('docen','var'), docen = []; end
if isempty(docen),      docen = true; end


%% reformat input + checks
nn = size(df0);
hz = reshape(df0,[prod(nn) 1]);
T1 = permute(T1(:),[2 3 1]);
T2 = permute(T2(:),[2 3 1]);
flip = permute(flip(:),[2 3 1]);
if any(T1<T2), warning('T1<T2'); end
if any(T2<tr), warning('T2<TR'); end
if any(T1<tr), warning('T1<TR'); end


%% off-resonance
phi = 2*pi*hz*tr;
if npc>1
    phi = bsxfun(@plus,phi,2*pi*(0:npc-1)/npc);
end
if docen, phi = phi-pi; end


%% calculate spatial response function
e1 = exp(-tr./T1);       % Eq. 5 wrong, see mrm18_244_1991_mizumoto.pdf
e2 = exp(-tr./T2);
c = e2.*(e1-1).*(1+cosd(flip));
d = (1-e1.*cosd(flip))-(e1-cosd(flip)).*e2.^2;

if length(d)==1
    srf = ((1-e1)*(1-e2*exp(-1i*phi))*sind(flip))./(c*cos(phi)+d);   % Eq.4
else
    srf = bsxfun(@times,(1-e1).*sind(flip),(1-bsxfun(@times,e2,exp(-1i*phi))))./...
        bsxfun(@plus,bsxfun(@times,c,cos(phi)),d);
end


%% weights
if nargout>1
    % normalise to get weights
    ww = bsxfun(@rdivide,conj(srf),sum(abs(srf),2));

    % reshape to original matrix size of df0 + #phase cycles
    % ww = reshape(ww,[nn npc]);
end


end      % ssfp_freq_resp.m
