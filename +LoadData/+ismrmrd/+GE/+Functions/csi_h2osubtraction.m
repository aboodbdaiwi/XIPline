function [spec] = csi_h2osubtraction(spec,ind)
%CSI_H2OSUBTRACTION Subtracts water signal using GE product method (CSI data) 
%  [spec] = csi_h2osubtraction(spec,ind,verb)
%        spec   Reconstructed data   size(spec) = [#s,#x,#y,#z,#t,#c]
%               Dim5: 1=spec,2=reference
%         ind   Index list (in spectral domain) for calculating weights
%
% 11/2015  Rolf Schulte
if nargin<1, help(mfilename); return; end

%% misc parameters + warnings
if ~exist('ind','var'),  ind = []; end
si = size(spec);
if length(si)==5, si = [si 1]; end
if length(si)~=6, error('spec not 6D (=%g)',length(si)); end
if si(5)~=2,
    warning('size(spec,5)(=%g) ~= 2',si(5));
    fprintf('Water subtraction requires reference spectrum stored as second echo\n');
end
if isempty(ind), ind = floor(si(1)/2) + (-25:24); end
if sum(ind)<10, warning('sum(ind)(=%g) < 10',sum(ind)); end


%% water subtraction
% [#s,#x,#y,#z,#t,#c]

% scale between spectrum and reference
fac = mean(abs(spec(ind,:,:,:,1,:))./abs(spec(ind,:,:,:,2,:)),1);

% subtract scaled water reference from spectrum
sws = spec(:,:,:,:,1,:) - spec(:,:,:,:,2,:).*repmat(fac,[size(spec,1) 1 1 1 1 1]);

% determine where subtraction made things worse
ind_worse = max(abs(sws),[],1)>max(abs(spec(:,:,:,:,1,:)),[],1);
% if any(ind_worse(:)), sum(ind_worse(:)), end

% use only spectra with improved water suppression
for l2=1:si(2),
    for l3=1:si(3),
        for l4=1:si(4),
            for l6=1:si(6),
                if ~ind_worse(1,l2,l3,l4,1,l6),
                    spec(:,l2,l3,l4,1,l6) = sws(:,l2,l3,l4,1,l6);
                end
            end
        end
    end
end

end
