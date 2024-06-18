function [b1map,bbabs] = blosi_b1map(b,h,threshold,do_circ_mask,verb,cocome)
%BLOSI_B1MAP  Extract B1+ map from images acquired with Bloch-Siegert pulse
%
% [b1map,bbabs] = blosi_b1map(b,h,threshold,do_circ_mask,verb,cocome)
%        b  Images acquired with +/- Bloch-Siegert frequency
%           size of images (mtx,mtx,#exc,#metabolites,#slice,#coils)
%           +BS in b(:,:,1:2:end,:,:,:)
%           -BS in b(:,:,2:2:end,:,:,:)
%           if acquired with fidall cv rfs_dudri=3
%           Ch1 in b(:,:,[1,2],:,:,:)
%           Ch2 in b(:,:,[3,4],:,:,:)
%        h  header
%           or KBS (Bloch-Siegert scaling factor from header2kbs.m)  [rad]
%threshold  if given, mask results
%  do_circ  apply circular mask (FOV of spirals in circle)
%     verb  Verbose (1+2) & Plotting (2) (default=2)
%   cocome  Coil combination method: 
%           0=off; 1=blosi_phase; 2=angle(sum(normalised coils))
%
%    b1map  Transmit B1 map                                   [uT]
%           [n1,n2,n3,nTX]
%    bbabs  RMS magnitude images [n1,n2,n3]
%
% See also HEADER2KBS.
% 10/2017 Rolf Schulte
if (nargin<1), help(mfilename); return; end;

if ~exist('threshold','var'), threshold = []; end
if ~exist('do_circ_mask','var'),do_circ_mask = []; end
if isempty(do_circ_mask),     do_circ_mask = true; end
if ~exist('verb','var'),      verb = []; end
if isempty(verb),             verb = 2; end
if ~exist('cocome','var'),    cocome = []; end
if isempty(cocome),           cocome = 2; end


%% misc parameters and checks
if verb>0, fprintf('assuming +/-BS in exc=1/2; reshaping data\n'); end
bp = b(:,:,1:2:end,:,:,:);
bn = b(:,:,2:2:end,:,:,:);

if any(size(bn)~=size(bp)), error('size(bn)~=size(bp)'); end
if isreal(b), error('isreal(b)'); end
[mtx1,mtx2,nexc,nmet,nslice,ncoils] = size(bn);
if ncoils==1, cocome = 0; end
if isstruct(h),
    if verb>0, fprintf('extracting KBS from header\n'); end
    KBS = header2kbs(h);
else
    if verb>0, fprintf('KBS = h\n'); end
    KBS = h;
end
if verb>0,
    fprintf('mtx1=%g, mtx2=%g, nexc=%g, nmet=%g, nslice=%g, ncoils=%g\n',...
        mtx1,mtx2,nexc,nmet,nslice,ncoils);
    fprintf('KBS=%g\n',KBS);
    fprintf('rfs_dudri=%g (0=normal;1=I,2=Q,3=alternate I-Q,4=I+Q)\n',...
        h.rdb_hdr.user35);
end


%% calculation of Bloch-Siegert phase difference
if (cocome==2)
    bn1 = bsxfun(@rdivide,bn,sum(abs(bn),6)/2);  % normalise magnitude for weighting
    bp1 = bsxfun(@rdivide,bp,sum(abs(bp),6)/2);  % normalise magnitude for weighting
    blosi_phase = angle(sum(conj(bn1).*bp1,6));
else
    blosi_phase = angle(conj(bn).*bp);
end

%% extend possible phases from 0..pi to 0..2*pi 
% -> add 2*pi to negative values; problematic with low SNR
% blosi_phase = blosi_phase + 2*pi.*(blosi_phase<0);


%% B1+ phase map
% if multi-Tx-channels (nexc>=4)
% phase difference between Tx channels
% all relative to Tx channel 1; ie b1phase of Tx channel 1 = 0
b1phase = [];
if nexc<2,
    fprintf('nexc<4; setting b1phase = []\n');
else
    if cocome==2
        b1phase = zeros(mtx1,mtx2,nexc,nmet,nslice);
    else
        b1phase = zeros(size(b));
    end
    for l=2:nexc,
        % phase difference between Tx channels for both BS+ and BS-
        pp = (l-1)*2+(1:2);
        if (cocome==2)
            btmp = bsxfun(@rdivide,b,sum(abs(b),6));
            tmp = angle(sum(conj(btmp(:,:,pp,:,:,:)).*btmp(:,:,1:2,:,:,:),6));
        else
            tmp = angle(conj(b(:,:,pp,:,:,:)).*b(:,:,1:2,:,:,:));
        end
        b1phase(:,:,pp,:,:,:) = tmp;
    end
end


%% averaging
if (cocome>0)
    bb = abs(sqrt(mean((conj(bn).*bn),6)+mean((conj(bp).*bp),6)));
else
    bb = abs(sqrt((conj(bn).*bn)+(conj(bp).*bp)));
end
if nmet>1
    warning('nmet(=%g)>1; averaging',nmet);
    blosi_phase = mean(blosi_phase,4);
    bb = sqrt(mean(bb.^2,4));
end
if (cocome==1)
    fprintf('Coil combination of blosi_phase (weighted by coil sensitivities)\n');
    ww = sqrt((conj(bn).*bn)+(conj(bp).*bp));
    ww = ww./repmat(sum(ww,6),[1 1 1 1 1 ncoils]);
    blosi_phase = sum(ww.*blosi_phase,6);
    if ~isempty(b1phase), 
        b1phase = sum(repmat(ww,[1 1 2 1 1 1]).*b1phase,6);
    end
end


%% magnitude B1+ map
b1map = real(sqrt(blosi_phase/2/KBS));


%% combine phase map + add to b1map
if ~isempty(b1phase),
    b1phase = (b1phase(:,:,1:2:end,:,:,:)+b1phase(:,:,2:2:end,:,:,:))/2;
    b1map = b1map.*exp(1i*b1phase);
end


%% masking
if ~isempty(threshold),
    threshold = threshold*(mean(bb(:))/2);
    mask = logical((blosi_phase>0).*(bb>threshold));
else
    mask = true(size(blosi_phase));
end
if do_circ_mask
    mask = logical(bsxfun(@times,mask,circular_mask([mtx1,mtx2])));
end
b1map = b1map.*mask;

%% reshaping
% size(b1map):  [mtx1,mtx2,nexc,nmet,nslice,ncoils] -> [mtx1,mtx2,mtx3,ntx]
tmp = zeros(mtx1,mtx2,nslice,nexc);
for l3=1:nexc
    for l5=1:nslice
        tmp(:,:,l5,l3) = b1map(:,:,l3,1,l5,1);
    end
end
b1map = tmp;

bbabs = sqrt(mean(bb.^2,3));
bbabs = reshape(bbabs,[mtx1,mtx2,nslice]);


%% plotting
if verb>1,
    figure(10); clf
    imagesc_row(abs(b1map),[0 13]);
    title('B1+ [uT]','FontSize',14);
    xlabel('slices'); ylabel('Tx coils')
    colorbar
    
    figure(11); clf
    % imagesc_row(reshape(bb,[mtx1,mtx2,nexc,nslice]));
    imagesc_row(bbabs);
    colormap gray
    title('Magnitude image','FontSize',14);
    xlabel('slices'); % ylabel('Tx coils')
    % colorbar
    
    if ~isempty(b1phase),
        figure(12); clf
        % imagesc_row(squeeze(b1phase),[-pi pi]);
        imagesc_row(angle(b1map(:,:,:,2:end)),[-pi pi]);
        title('B1 phase [rad]','FontSize',14);
        xlabel('slices'); ylabel('Tx coils')
        colorbar
    end
end

%% summarise results
if verb>0,
    fprintf('BloSi Phase = %.3g +- %.3g [deg]\n',...
        mean(blosi_phase(mask))*180/pi,std(blosi_phase(mask))*180/pi);
    fprintf('B1+ = %.3g +- %.3g [uT]\n',mean(b1map(mask)),std(b1map(mask)));
end

