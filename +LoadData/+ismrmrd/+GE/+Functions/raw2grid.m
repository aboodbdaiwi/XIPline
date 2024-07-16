function [out_data] = raw2grid(data,chp,k,shft,bw_filt,dec,ind,delay,...
    bw_acq,plt,ind_pc)
%out_data = raw2cg(data,chp,k,shft,bw_filt,dec,ind,delay,bw_acq,plt,ind_pc)
% Reformat data from P-files into format used for recons
%
%     data  dim1=#excitations
%           dim2=#samples (acq in one view)
%           dim3=1 (#phases)
%           dim4=#echoes
%           dim5=#slices
%           dim6=#coils
%      chp  chop (phase cycling)
%           if true -> multiply each other echo with -1 (default=true)
%        k  k-space positions (req for shft) (<0.5)
%     shft  timings for linear phase shift [y x] or [imag(k) real(k)]
%  bw_filt  Filter sampling BW (fft -> lowpass -> ifft): 0 < bw_filt < 1
%      dec  downsample (take every dec'th point) dec=2,3,4... 
%           filter with bw<=1/dec first
%      ind  index list (to select data to use) logical
%           dimension (1,kspace)
%    delay  Gradient delay time [s] (somewhere around 60d-6)        [s]
%   bw_acq  (Full) acq bandwidth (=2*(+/-BW)) (req for delay)       [Hz]
%      plt  Plotting (default=true)
%   ind_pc  Index for phase correction applied along dim2
%           ind_pc{l}{1}->extract phase; ind_pc{l}{2}->apply phase to
%
% out_data  dim1=#coils
%           dim2=#indexed kspace data
%           dim3=#excitations (if size(ind,1)==1)
%                #time steps  (otherwise)
%           dim4=#slices
%
% 7/2019 Rolf Schulte
%
% See also AK_DOWNSAMPLE,ISCHOP.
if (nargin<1), help(mfilename); return; end


%% input defaults + checks
[nexc,ns,n3,n4,nslices,ncoils] = size(data);
if n3~=1,           warning('size(data,3)(=%g) ~= 1',n3); end
if ~islogical(chp), error('~islogical(chp)'); end
if ~exist('k','var'),        k = []; end
if ~exist('shft','var'),     shft = []; end
if ~exist('bw_filt','var'),  bw_filt = []; end
if ~exist('dec','var'),      dec = []; end
if isempty(dec),             dec = 1; end
if dec<1, error('dec<1'); end
if ~exist('ind','var'),      ind = []; end
if ~isempty(ind), ns1 = sum(ind(:)); end
if ~exist('delay','var'),    delay = []; end
if ~exist('plt','var'),      plt = []; end
if isempty(plt),             plt = true; end
if ~exist('ind_pc','var'),   ind_pc = []; end
mak = max(abs(k(:)));
if mak>0.6, warning('max(abs(k(:)))(=%g)>0.6', mak); end
if mak>5,    error('max(abs(k(:)))(=%g)>5',mak); end


%% shifting data
if ~isempty(shft)
    if isempty(k) 
        warning('k empty/not existing'); 
    end
    if size(k,1)~=1, error('size(k,1)~=1'); end
    shft = shft(:);
    if ~(length(shft)==2 || length(shft)==3)
        error('length(shft) not 2 or 3'); 
    end
    if ~isreal(k)       % convert complex k-space to [1,npts,2(x+y)]
        if size(k,3)~=1
            error('k complex and size(k,3)~=1');
        end
        tmp = zeros(1,size(k,2),2);
        tmp(1,:,1) = real(k); tmp(1,:,2) = imag(k);
        k = tmp;
    end
    if length(size(k))~=3, error('length(size(k))~=3'); end
    if size(k,3)~=length(shft)
        error('size(k,3)~=length(shft)'); 
    end
    if size(k,2)~=ns1
        warning('size(k,2)(=%g)~=ns1(=%g)',size(k,2),ns1);
    end
end
if dec~=1, warning('dec not anymore implemented; ignoring'); end
if ~isempty(bw_filt)
    warning('bw_filt not anymore implemented; ignoring'); 
end


%% plotting
if plt
    figure; 
    subplot(4,1,1);
    [tmp,l1] = max(abs(data(:,3,1)));
    fprintf('Chosing l1=%g for displaying data(l1,:,...)\n',l1);
    plot_complex(linspace(-1,1,ns),fftshift(ifft(data(l1,:,1),[],2),2));
    set(gca,'xtickLabel','');title('');ylabel('original');
    pos=get(gca,'Position');set(gca,'Position',pos.*[1 1 1 1.2]);
end


%% unchop data first
if chp && nexc>1
%     for l1=2:2:nexc,
%         data(l1,:) = -data(l1,:);
%     end
    data(2:2:nexc,:) = -data(2:2:nexc,:);
end


%% adjust gradient delay time: shift data by linear phase in f-domain
if ~isempty(delay)
    if delay~=0
        fprintf('Compensating for gradient delay = %g [us]\n',...
            delay*1d6);
        phafu = fftshift(exp(2*pi*1i*(-(ns/2):(ns/2-1))/ns*bw_acq*delay));
        % data = fft(ifft(data,[],2).*repmat(phafu,[nexc 1 1 n4 nslices ncoils]),[],2);
        data = fft(bsxfun(@times,ifft(data,[],2),phafu),[],2);
    end
end


%% Extract phase from data and phase correct data
if ~isempty(ind_pc)
    fprintf('Performing a phase correction\n');
    for l1=1:nexc
        for l3=1:n3
            for l4=1:n4
                for l5=1:nslices
                    for l6=1:ncoils
                        for l2=1:length(ind_pc{1})
                            pc = mean(unwrap(angle(data(l1,ind_pc{1}{l2},l3,l4,l5,l6))));
                            % mc = mean(abs(data(l1,ind_pc{1}{l2},l3,l4,l5,l6)));
                            data(l1,ind_pc{2}{l2},l3,l4,l5,l6) = ...
                                data(l1,ind_pc{2}{l2},l3,l4,l5,l6)*exp(-1i*pc);
                        end
                    end
                end
            end
        end
    end
end


%% selecting indexed data
if ~isempty(ind)
    if size(data,1)<size(ind,1) 
        error('size(data,1)(%g)<size(ind,1)(%g): increase cv4',...
            size(data,1),size(ind,1)); 
    end
    if size(data,2)~=size(ind,2)
        fprintf('size(data,2)=%g\nsize(ind,2)=%g\n',...
            size(data,2),size(ind,2));
        if size(data,2)>size(ind,2)
            warning('size(data,2)>size(ind,2); adapting/increasing ind');
            tmp = false(size(ind,1),size(data,2)); 
            tmp(:,1:size(ind,2)) = ind;
            ind = tmp;
        else
            fprintf('size(data,2)(%g)<size(ind,2)(%g): truncating ind\n',...
                size(data,2),size(ind,2));
            tmp = ind(:,(size(data,2)+1):end);
            fprintf('sum(ind)=%g; discarding %g actual data pts\n',...
                sum(ind(:)), sum(tmp(:)));
            if sum(tmp(:))>0
                warning('truncation contains actual data (%g pts)\n',...
                    sum(tmp(:)));
                fprintf('Increase sampling (opuser1)!!!\n');
                fprintf('Zero-filling data instead\n');
                dtmp = complex(zeros(nexc,size(ind,2),n3,n4,nslices,ncoils,class(data)));
                dtmp(:,1:ns,:,:,:,:) = data;
                data = dtmp;
                clear dtmp
            else
                ind = ind(:,1:size(data,2));
            end
        end
    end
    fprintf('Selecting only indexed data (%g / %g)\n',...
        sum(ind(:)),size(data,2)*size(ind,1));
    
    nintl = size(ind,1);     % #of interleaves for encoding one image
    if nintl==1              % every repetition same ind (eg IDEAL reco)
        data = data(:,ind,:,:,:,:);
    else
        nrep = floor(nexc/nintl);      % #time steps
        dtmp = complex(zeros(nrep,ns1,1,n4,nslices,ncoils,class(data)));

        for lr=1:nrep
            for l4=1:n4
                for lsli=1:nslices
                    for lc=1:ncoils
                        tmp = data(((lr-1)*nintl+(1:nintl)),:,1,l4,lsli,lc).';
                        dtmp(lr,:,1,l4,lsli,lc) = tmp(ind.');
                    end
                end
            end
        end
        data = dtmp;
        nexc = nrep;
    end
    if any(isnan(data(:))) || isempty(data) || any(isinf(data(:)))
        error('data contains nan,inf or empty');
    end
end


%% plotting
if plt
    subplot(4,1,2);
    plot_complex(linspace(-1,1,size(data,2)),fftshift(ifft(data(1,:,1),[],2),2));
    set(gca,'xtickLabel','');title('');ylabel('delayed + indexed');
    pos=get(gca,'Position');set(gca,'Position',pos.*[1 1 1 1.2]);
end


%% data modulation for spatial shifting
if ~isempty(shft)
    switch length(shft)
        case 2
            fprintf('Off-isocenter data modulation shft=(%g,%g)\n',...
                shft(1),shft(2));
            tmp = exp(2*pi*1i*k(1,:,1)*shft(1) + 2*pi*1i*k(1,:,2)*shft(2));
        case 3
            fprintf('Off-isocenter (3D) data modulation shft=(%g,%g,%g)\n',...
                shft(1),shft(2),shft(3));
             tmp = exp(2*pi*1i*k(1,:,1)*shft(1) + ...
                 2*pi*1i*k(1,:,2)*shft(2) + ...
                 2*pi*1i*k(1,:,3)*shft(3));
    end
    % tmp = (exp(2*pi*1i*real(k)*shft(1) + 2*pi*1i*imag(k)*shft(2)));
    % data = data.*repmat(tmp,[nexc 1 1 n4 nslices ncoils]);
    data = bsxfun(@times,data,tmp);
end


%% plotting
if plt
    subplot(4,1,3); 
    plot_complex(linspace(-1,1,ns1),fftshift(ifft(data(1,:,1),[],2),2));
    set(gca,'xtickLabel','');title('');ylabel('modulated');
    pos=get(gca,'Position');set(gca,'Position',pos.*[1 1 1 1.2]);
end


%% lowpass filter data; reduce (noise) effective sampling BW
% skipping, as it has little effect on SNR (even with bw_filt=0.33)
% oversampling in gridding.m or gridding3d.m seems to do the job as well
% would need to be adapted for size(ind,1)>1 -> fft of interleaves only
% if ~isempty(bw_filt)
%     ns2 = size(data,2);
%     xx = floor(bw_filt*ns2/2+0.5);
%     lbs = zeros(1,ns2); lbs(floor(ns2/2)+(-xx:xx-1)+1)=1;
%     lbs = ifftshift(real(ifft(ifftshift(gausswin(ns2,ns2/70).'.*fftshift(fft(lbs))))));
%     for li=1:nintl
%         data = fft(ifft(data,[],2).*repmat(lbs,[nexc 1 1 1 nslices ncoils]),[],2);
%     end
% end
% if ((dec>1) || ~isempty(bw_filt))
%     fprintf('Apply (Cartesian) lowpass filter: bw_filt = %g\n',bw_filt);
%     data = ak_downsample(data,dec,bw_filt,false);
%     ns2 = size(data,2);
% else
    ns2 = ns1;
% end


%% plotting
if plt
    subplot(4,1,4);
    plot_complex(linspace(-1,1,ns2),fftshift(ifft(data(1,:,1),[],2),2));
    set(gca,'xtickLabel','');title('');ylabel('filtered -> final');
    pos=get(gca,'Position');set(gca,'Position',pos.*[1 1 1 1.2]);
end


%% reformat data
fprintf('Reformat data\n');
out_data = complex(zeros(ncoils,ns2,nexc*n4,nslices,class(data)));
for lc=1:ncoils
    for le=1:nexc
        for lt=1:n4    % number of time steps (rfs_nreps)
            for ls=1:nslices
                % out_data(lc,:,le,ls) = data(le,:,1,1,ls,lc);
                out_data(lc,:,(le+(lt-1)*nexc),ls) = data(le,:,1,lt,ls,lc);
            end
        end
    end
end


%% checks
if any(isnan(out_data(:))), error('out_data contains NaN'); end
if any(isinf(out_data(:))), error('out_data contains inf'); end


end    % main function raw2grid.m
