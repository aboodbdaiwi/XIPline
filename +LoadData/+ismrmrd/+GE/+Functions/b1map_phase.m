function [b1map,bbabs] = b1map_phase(bb,mask,verb)
%B1MAP_PHASE  Determine phase between Tx channels
%
%[b1map,bbabs] = b1map_phase(bb,mask,verb)
%      bb  Images acquired with different Tx channels (dim=4)
%          image sizes = [n1,n2,n3,nechoes,ncoils]
%    mask  Mask for selecting ROI (logical)
%          if length(mask)==1: threshold multiplicator for generating mask
%          (default=[]; no masking)
%    verb  Verbose mode: 0=off; 0=off; 1=verbose; 
%          2=verb+imagesc_row; 3=verb+imagesc_ind3d (default)
%
%   b1map  B1+ map (complex with mask (0 or 1)); 
%          phases relative to channel1
%          size = [n1,n2,n3,nTX]
%   bbabs  RMS magnitude images [n1,n2,n3,nTX]
%
% 10/2017 Rolf Schulte
if nargin<1, help(mfilename); return; end


%% defaults
[n1,n2,n3,n4,nc] = size(bb);
if ~exist('mask','var'),    mask = []; end
if ~exist('verb','var'),    verb = []; end
if isempty(verb),           verb = 3; end


%% parameter checks
if isreal(bb), warning('bb not complex'); end
if n4<2, warning('n4(=%g)<2',n4); end
if length(size(bb))>5, warning('length(size(bb))(=%g)>5',length(size(bb))); end
if length(mask)>1
    if ~islogical(mask), error('mask not logical'); end
    if any(size(mask,1)~=[n1 n2 n3]), error('size(mask)~=[n1 n2 n3]'); end
end


%% relative B1+ magnitude map
b = abs(sqrt(mean(mean(conj(bb).*bb,5),4)));
b = b/max(b(:));             % normalise


%% generate mask, if length(mask)==1
if length(mask)==1
    if verb>0, fprintf('Generating mask\n'); end
    b_vec = sort(b(:)); nnn = length(b_vec);
    threshold = mask*mean(b(1:floor(0.2*nnn)));
    if verb>0, fprintf('  threshold = %g\n',threshold); end
    mask = b>threshold;
    mask = imfill(mask,'holes');
    % mask = ~imfill(~mask,'holes');
end
if isempty(mask), mask = true(n1,n2,n3); end
if verb>0
    suma = sum(mask(:));
    fprintf('Mask contains %g voxel (%.4g%%)\n',suma,suma/(n1*n2*n3)*100);
end
b1map = repmat(single(mask),[1 1 1 n4]);


%% relative B1+ phase map
% b1p_phase = angle(bsxfun(@times,conj(bb(:,:,:,1,:)),bb));
% if nc>1                      % weighted coil combination
%     if verb>0, fprintf('Weighted coil combination\n'); end
%     weight = abs(sqrt(mean(conj(bb).*bb,4)));
%     weight = repmat(weight./repmat(sum(weight,5),[1 1 1 1 nc]),[1 1 1 n4 1]);
%     b1p_phase = sum(weight.*b1p_phase,5);
% end
if nc>1       % weighted phase average
    bb = bsxfun(@rdivide,bb,sum(abs(bb),5));     % normalise coils
    b1p_phase = angle(sum(bsxfun(@times,conj(bb(:,:,:,1,:)),bb),5));
else
    b1p_phase = angle(bsxfun(@times,conj(bb(:,:,:,1,:)),bb));
end


%% print info
if verb>0
    fprintf('Relative B1+ map:\n');
    for l4=2:n4
        tmp = b1p_phase(:,:,:,l4);
        fprintf('\tphase (channel %g) = %g +- %g\n',...
            l4,mean(tmp(mask)),std(tmp(mask)));
    end
end


%% add phase to b1map
b1map = b1map.*exp(1i*b1p_phase);


%% plotting
if verb>1,
    figure(20); clf;
    if verb==2, imagesc_row(b,[0 1],'',true);
    else, imagesc_ind3d(b,'',[0 1]); end
    title('B');
    colorbar;

    figure(21); clf;
    if verb==2, imagesc_row(mask,[0 1],'',true);
    else, imagesc_ind3d(mask,'',[0 1]); end
    title('B1+ Mask');
    colorbar;

    for l4=2:n4
        figure(20+l4); clf;
        if verb==2, imagesc_row(angle(b1map(:,:,:,l4)),[-pi pi],'',true);
        else, imagesc_ind3d(angle(b1map(:,:,:,l4)),'',[-pi pi]); end
        title(['B1+ Phase channel ' num2str(l4)]);
        colorbar;
    end
    drawnow;
end

%% RMS image
if nargout>1
    bbabs = sqrt(mean(b.^2,4));
end


end  % main function b1map_phase.m

