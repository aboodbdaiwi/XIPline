function b1map = relative_b1map(bb,mask,do_fiex,verb)
%RELATIVE_B1MAP  Determine relative B1+ maps from two images acquired with 
%   transmitting only on a single Tx channel
%    b1map = relative_b1map(bb,mask,do_fiex,verb)
%        bb  Images acquired with Tx on I&Q (dim=4) and Rx on I&Q (dim=5)
%            image sizes = [n1,n2,n3,nechoes,ncoils]
%      mask  Mask for selecting ROI (logical)
%            if length(mask)==1: threshold multiplicator for generating mask
%            (default=[]; no masking)
%   do_fiex  Voxel-wise polynomial smoothing (fiex3d.m)   (default=false)
%      verb  Verbose mode (0=off; 1=verbose; 2=verbose+plotting (default))
%     b1map  B1+ map (complex)
%            size = [n1,n2,n3,2]
%
% 10/2017 Rolf Schulte
if nargin<1, help(mfilename); return; end


%% defaults
[n1,n2,n3,n4,nc] = size(bb);
if ~exist('mask','var'),    mask = []; end
if ~exist('do_fiex','var'), do_fiex = []; end
if isempty(do_fiex),        do_fiex = false; end
if ~exist('verb','var'),    verb = []; end
if isempty(verb),           verb = 2; end


%% parameter checks
if isreal(bb), warning('bb not complex'); end
if n4~=2, warning('n4(=%g)~=2',n4); end
if nc~=2, warning('nc(=%g)~=2',nc); end
if length(mask)>1
    if ~islogical(mask), error('mask not logical'); end
    if any(size(mask,1)~=[n1 n2 n3]), error('size(mask)~=[n1 n2 n3]'); end
end


%% relative B1+ magnitude map
b1map(:,:,:,1) = sqrt(abs(bb(:,:,:,1,1)));
b1map(:,:,:,2) = sqrt(abs(bb(:,:,:,2,2)));
b1map = b1map/max(b1map(:));           % normalise


%% generate mask, if length(mask)==1
if length(mask)==1
    if verb>0, fprintf('Generating mask\n'); end
    b1map_vec = sort(b1map(:)); nnn = length(b1map_vec);
    threshold = mask*mean(b1map(1:floor(0.2*nnn)));
    if verb>0, fprintf('  threshold = %g\n',threshold); end
    mask = ((b1map(:,:,:,1)>threshold) & (b1map(:,:,:,2)>threshold));
    mask = imfill(mask,'holes');
    mask = ~imfill(~mask,'holes');
else
    mask = true(n1,n2,n3);
end
if verb>0,
    suma = sum(mask(:));
    fprintf('Mask contains %g voxel (%.4g%%)\n',suma,suma/(n1*n2*n3)*100);
end
b1map = bsxfun(@times,b1map,mask);


%% relative B1+ phase map
% b1p_phase = angle(conj(bb(:,:,:,1,:)).*bb(:,:,:,2,:));
% if nc>1,  % weighted coil combination
%     weight = sqrt(mean(conj(bb).*bb,4));
%     weight = weight./repmat(sum(weight,5),[1 1 1 1 nc]);
%     b1p_phase = sum(weight.*b1p_phase,5);
% end
if nc>1       % weighted phase average
    bb = bsxfun(@rdivide,bb,sum(abs(bb),5));     % normalise coils
    b1p_phase = angle(sum(bsxfun(@times,conj(bb(:,:,:,1,:)),bb),5));
else
    b1p_phase = angle(bsxfun(@times,conj(bb(:,:,:,1,:)),bb));
end


%% apply voxel-wise polynomial smoothing
if do_fiex,
    for l4=1:2,
        b1map(:,:,:,l4) = fiex3d(b1map(:,:,:,l4),[],0);
    end
end


%% print info
if verb>0,
    fprintf('Relative B1+ map:\n');
    fprintf('\tmagnitude = %g +- %g\n',mean(b1map(mask)),std(b1map(mask)));
    fprintf('\tphase     = %g +- %g\n',mean(b1p_phase(mask)),std(b1p_phase(mask)));
end

%% add phase to Q (channel 2)
b1map(:,:,:,2) = b1map(:,:,:,2).*exp(1i*b1p_phase);


%% plotting
if verb>1,
    figure(11); clf;
    imagesc_row(abs(b1map(:,:,:,1)),[0 1],'',true);
    title('B1+ Magnitude Map (I)');
    colorbar;

    figure(12); clf;
    imagesc_row(abs(b1map(:,:,:,2)),[0 1],'',true);
    title('B1+ Magnitude Map (Q)');
    colorbar;

    figure(13); clf;
    imagesc_row(angle(b1map(:,:,:,2)),[-pi pi],'',true);
    title('B1+ Phase Map (Q)');
    colorbar;
    
    figure(21); clf;
    imagesc_ind3d(abs(b1map(:,:,:,1)),'',[0 1]);
    title('B1+ Magnitude Map (I)');
    colorbar;

    figure(22); clf;
    imagesc_ind3d(abs(b1map(:,:,:,2)),'',[0 1]);
    title('B1+ Magnitude Map (Q)');
    colorbar;

    figure(33); clf;
    imagesc_ind3d(angle(b1map(:,:,:,2)),'',[-pi pi]);
    title('B1+ Phase Map (Q)');
    colorbar;
    drawnow;
end
