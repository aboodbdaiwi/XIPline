function data = csi_ecc(data,method)
%CSI_ECC Eddy current correction for CSI signal using phase deconvolution
%  Dividing data (with water suppression) through phase of reference (w/o
%  water suppression to remove phase errors from eddy currents.
%
% data = csi_ecc(data,method)
%    data  Spectrally time-domain, spatially image-domain
%          size(data) = [#s,#x,#y,#z,#t,#c]
%           Dim5: 1=spec,2=reference
%  method  1=take raw phase of reference signal
%          2=fit/smooth raw phase first
%
%  Literature:
%  In Vivo Proton Spectroscopy in Presence of Eddy Currents.
%  U Klose, MRM 14, 26-30 (1990).
%
% 11/2015  Rolf Schulte
if nargin<1, help(mfilename); return; end

%% misc parameters + checks
ind0 = 3:20;
if ~exist('method','var'), method = []; end
if isempty(method),        method = 1; end

si = size(data);
if length(si)==5, si = [si 1]; end
if length(si)~=6, error('spec not 6D (=%g)',length(si)); end
if si(5)~=2,
    warning('size(data,5)(=%g) ~= 2',si(5));
    fprintf('Water subtraction requires reference spectrum stored as second echo\n');
end

%%

switch method,
    case 1,
        % subtraction of the disturbing phase is equal to a point-wise
        %  division of the phase factor
        pha = data(:,:,:,:,2,:)./abs(data(:,:,:,:,2,:));
        data(:,:,:,:,1,:) = data(:,:,:,:,1,:)./pha;
        if ~isempty(ind0),
            pha0 = mean(unwrap(angle(pha(ind0,:,:,:,1,:)),[],1),1);
            data(:,:,:,:,2,:) = data(:,:,:,:,2,:).*repmat(exp(-1i*pha0),[si(1) 1 1 1 1 1]);
        end
    case 2,
        warning('not yet properly implemented');
%         pha = unwrap(angle(data(:,:,:,:,2,:)));
%         x_all = (1:si(2));
%         seg = [10 600];
%         x_fit = (seg(1):seg(2));
%         pol2 = polyfit(x_fit,pha(x_fit),3);
%         pha_fit = polyval(pol2,x_all);
%         fid = fid.*repmat(exp(-1i*pha_fit),[si(1),1]);
%         
%         if plt,
%             clf;
%             subplot(2,1,1);
%             plot(x_all,pha,'b',x_all,pha_fit,'r');
%             legend('phase of reference FID','polynomial fit to phase');
%             subplot(2,1,2);
%             plot(x_all,unwrap(angle(fid)),'b');
%             figure;
%         end
    otherwise, warning('method(=%g) not existing',method);
end


end