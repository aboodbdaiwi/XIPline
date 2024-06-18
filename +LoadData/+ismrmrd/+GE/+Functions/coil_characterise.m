function coil_characterise(d,h,out_fname,pix_fname,fig_off)
%COIL_CHARACTERISE  Determine SNR for fidcsi acquisition
% if multi-channel -> individual SNR and after SNR-optimal coil combination
% if noise and multi-channel 
%    -> calculate noise correlation + covariance matrix
%
% coil_characterise(d,h,out_fname,pix_fname)
%         d  data (or P-file name)
%         h  header (or empty if name is given)
% out_fname  write output to file
% pix_fname  save picture as png
%
% Literature (noise covariance): ISMRM2011,#4548.
% 10/2019 Rolf Schulte
if (nargin<1), help(mfilename); return; end


%% default parameters
plt = true;
if ~exist('fig_off','var'),   fig_off = []; end
if isempty(fig_off),          fig_off = 0; end
if ~exist('out_fname','var'), out_fname = []; end
if ~exist('pix_fname','var'), pix_fname = []; end
if ~isempty(pix_fname)
    if strcmpi(pix_fname,'false'), plt = false; end
    if ~isempty(regexpi(pix_fname,'\.png'))
        pix_fname = pix_fname(1:end-4);
    end
    if ~isempty(regexpi(pix_fname,'\.7'))
        pix_fname = pix_fname(1:end-2);
    end
end
if isempty(out_fname)
    outid = 1;                         % direct to stdout
else
    outid = fopen(out_fname,'w+');     % direct to file
end

ind_signal = (3:8);          % points of FID used for signal (area)
pts_noise = 200;             % #points at end of FID for noise
pts_svds = 32;               % #largest points for svds coil phase calc


%% reading in data, if pfile name is given
if ~isnumeric(d)
    if exist(d,'file')
        pfname = d;
        [d,h] = read_p(d);
    else
        fprintf(outid,'Warning: strange input d\n');       
    end
end


%% misc
[nexc,npts,n3,n4,nslices,ncoils] = size(d);
if n3>1
    warning('n3(=%d)>1; ignoring/truncating dim3',n3);
    d = d(:,:,1,:,:,:);
end
if n4>1
    warning('n4(=%d)>1; ignoring/truncating dim4',n4);
    d = d(:,:,:,1,:,:);
end
if pts_noise>npts
    warning('pts_noise(=%d)>npts(=%d)',pts_noise,npts);
    pts_noise = npts;
end
if pts_svds>npts
    warning('pts_svds(=%d)>npts(=%d)',pts_svds,npts);
    pts_svds = npts;
end
    
ind_noise = ((npts-pts_noise+1):npts);
dd = d;

chp = ischop(h);
bw = h.rdb_hdr.user0;                  % (full) BW [Hz]
hz = (-npts/2:npts/2-1)/npts*bw;       % freq axis
t = (1:size(d,2))/bw*1d3;              % time-axis [ms]
if chp, d(2:2:end,:,:,:,:,:) = -d(2:2:end,:,:,:,:,:); end
% d = mean(d,1);


%% SNR calculation
if ncoils>1   % individual for getting the weights
    fprintf(outid,'slice coil signal noise SNR(1exc)\n');
    snr = zeros(nslices,ncoils);
    
    for ls=1:nslices
        for lc=1:ncoils
            % determine SNR
            signal = mean(abs(mean(d(:,ind_signal,1,1,ls,lc),1)),2);
            tmp = d(:,ind_noise,1,1,ls,lc); tmp = tmp(:).';
            noise = std(tmp);
            snr(ls,lc) = signal/noise;
            fprintf(outid,'%g\t%g\t%g\t%g\t%g\n',ls,lc,signal,noise,snr(ls,lc));
        end
        if plt
            figure(ls+fig_off); clf
            subplot(3,1,1);
            plot(t,squeeze(abs(mean(d(:,:,1,1,ls,:),1))));
            ylabel('coil FIDs (abs)'); axis tight;
            title(sprintf('SNR=%s', num2str(snr(ls,:),'%.2g ')));
        end
        fprintf(outid,'%g sqrt(sum((snr(ls,:).^2),2))=%g\n',ls,...
            sqrt(sum((snr(ls,:).^2),2)));
    end
end


%% SNR-optimal coil combination
if ncoils>1
    ww = snr./(sum(snr,2)*ones(1,ncoils));  % coil weighting
    
    % pc = mean(unwrap(angle(mean(d(1,ind_pts,1,1,:,:),1))),2);
    % pc = reshape(pc,[nslices ncoils]);      % phase correction
    % cocobi = permute(ww.*exp(-1i*pc),[3 4 5 6 1 2]);
    
    pc = zeros(nslices,ncoils);
    tmp = mean(d(:,:,:,:,ls,:),1); tmp = mean(abs(tmp.*conj(tmp)),6);
    [~,ii] = sort(tmp,'descend');
    ii = ii(1,1:pts_svds);
    for ls=1:nslices
        A = reshape(mean(d(:,ii,1,1,ls,:),1),[pts_svds ncoils]);
        [~,~,V] = svds(A,1);
        pc(ls,:) = V./abs(V);
        % pc(ls,:) = V;
    end
    % ww = ones(size(ww));
    cocobi = permute(ww.*pc,[3 4 5 6 1 2]);
    
    d = bsxfun(@times,d,cocobi);
    d = sum(d,6);
end
d = reshape(d,[nexc npts nslices]);

snrCMB = NaN(1,nslices);
if ncoils>1, fprintf(outid,'\nSNR with SNR-optimal coil combination\n'); end


%% loop through slices
for ls=1:nslices
    signal = mean(abs(mean(d(:,ind_signal,ls),1)),2);
    tmp = d(:,ind_noise,ls); tmp = tmp(:).';
    % tmp = mean(d(:,ind_noise,ls),1); tmp = tmp(:).';
    noise = std(tmp);
    snrCMB(1,ls) = signal/noise;
    fprintf(outid,'slice=%g signal=%g noise=%g ',ls,signal,noise);
    fprintf(outid,'SNR(1exc)=%g\n',snrCMB(1,ls));
    if plt 
        figure(ls+fig_off);
        if ncoils>1
            subplot(3,1,2);
        else
            subplot(2,1,1); 
        end
        plot(t,abs(mean(d(:,:,ls),1)).','b',...
            t(ind_signal),abs(mean(d(:,ind_signal,ls),1)).','gx',...
            t(ind_noise), abs(mean(d(:,ind_noise,ls),1)).','mx');
        ylabel('FID (abs)'); axis tight; xlabel('time [ms]');
        title(sprintf('SNR=%.3g (1 exc); pulse=%g; thick=%g [mm]', ...
            snrCMB(1,ls),h.image.user14,h.image.slthick));
        if ncoils>1
            subplot(3,1,3);
        else
            subplot(2,1,2); 
        end
        plot(hz,abs(fftshift(fft(conj(mean(d(:,:,ls),1)),[],2),2)),'b'); 
        xlabel('freq [Hz]'); ylabel('Spectrum'); axis tight;
        if exist('pfname','var'), set(gcf,'name',pfname(end-7:end)); end

        % print figures
        if ~isempty(pix_fname)
            print('-dpng',[pix_fname '_' num2str(ls) '.png']);
        end
    end
end


%% if multi-channel noise scan -> noise correlation matrix
if (ncoils>1) && all(snrCMB<1)
    fprintf(outid,'\nAttention: multi-channel coil and all noise:\n');
    fprintf(outid,'   calculating noise correlation + covariance\n');
    
    tmp = reshape(dd,[nexc*npts*nslices ncoils]);
    noisecor = abs(corrcoef(tmp));     % noise correlation
    noisecov = abs(tmp'*tmp);          % noise covariance
    noisecov = noisecov/mean(diag(noisecov));
    
    % print info
    fprintf(outid,'noisecor = \n');
    tmp = num2str(noisecor,'%.3g  ');
    for lc=1:ncoils, fprintf(outid,'%s\n',tmp(lc,:)); end
    fprintf(outid,'noisecov = \n');
    tmp = num2str(noisecov,'%.3g  ');
    for lc=1:ncoils, fprintf(outid,'%s\n',tmp(lc,:)); end
    
    fprintf(outid,'\nnoise correlation: max=%.2g%%; min=%.2g%%; mean=%.2g%%\n\n',...
        max(max(noisecor-eye(ncoils)))*100,min(min(noisecor))*100,...
        (sum(sum(noisecor))-ncoils)/(ncoils*(ncoils-1))*100);
    
    fprintf(outid,'noise covariance:\n');
    fprintf(outid,'\ton-diagonal:  max=%.3g; min=%.3g; mean=%.3g; std=%.3g\n',...
        max(diag(noisecov)),min(diag(noisecov)),...
        mean(diag(noisecov)),std(diag(noisecov)));
    tmp = noisecov + diag(inf(ncoils,1));
    tmp = tmp(~isinf(tmp));
    fprintf(outid,'\toff-diagonal: max=%.3g; min=%.3g; mean=%.3g; std=%.3g\n',...
        max(tmp),min(tmp),mean(tmp),std(tmp));
    
    % plotting
    if plt
        figure(ls+10); clf
        imagesc(noisecor,[0 1]); colorbar;
        title('Noise Correlation')
        % print figures
        if ~isempty(pix_fname)
            print('-dpng',[pix_fname '_' num2str(ls) '_noisecor.png']);
        end
        
        figure(ls+20); clf
        mm = (2+max(noisecov(:)))/2;
        imagesc(noisecov,[0 mm]); colorbar;
        title('Noise Covariance')
        % print figures
        if ~isempty(pix_fname)
            print('-dpng',[pix_fname '_' num2str(ls) '_noisecov.png']);
        end
    end
    
end


%% close file output
if ~isempty(out_fname), fclose(outid); end


end      % main function coil_characterise
