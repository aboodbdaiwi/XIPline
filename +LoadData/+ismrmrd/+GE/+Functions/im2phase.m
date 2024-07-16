function bp = im2phase(bb,iref)
%IM2PHASE  Coil-combined phase of images
%    bp = im2phase(bb,iref)
%    bb   Complex images           [#x,#y,#z,#time_steps,#coils]
%  iref   Index of reference image (subtract this phase from rest;time_step)
%    bp   Phase of coil combined images [#x,#y,#time_steps]
%
% 1/2020 Rolf Schulte
if (nargin<1), help(mfilename); return; end

if ~exist('iref','var'), iref = []; end
if isempty(iref),        iref = 1; end

nn = size(bb);
if length(nn)~=5
    warning('bb not 5D');
    bp = angle(bb);
    return 
end
if nn(4)==1
    warning('nn(4)==1; returning bp=[]');
    bp = [];
    return
end


%% coil combination + phase calculation
bb = bsxfun(@rdivide,bb,sum(abs(bb),5)); % normalise magnitude for weighting
% bp = angle(sum(bb.*conj(bb(:,:,:,iref,:,:)),5));
bp = angle(sum(bsxfun(@times,bb,conj(bb(:,:,:,iref,:,:))),5));
% subtract phase of reference image
% weighted coil combination
% calculate phase


end   % end main function im2phase.m
