function snrCMB=calc_snr(d,h)
%CALC_SNR  Calculate SNR from raw FID data
% snr = calc_snr(d,h)
%   d   Raw data (size=[nexc,npts,tmp1,tmp2,nslices,ncoils])
%   h   GE p-file header structure
% snr   SNR
%
%  3/2023 Rolf Schulte
if nargin<1, help(mfilename); return; end


%% default parameters
outid = 1;
ind_pts = (3:35);
ind_signal = (3:15);    % points of FID used for signal (area)
pts_noise = 100;        % #points at end of FID for noise
[nexc,npts,n3,n4,nslices,ncoils] = size(d);
if n3~=1, warning('n3(=%g)~=1',n3); end
if n4~=1, warning('n4(=%g)~=1',n4); end
ind_noise = ((npts-pts_noise):npts);

chp = ischop(h);
if chp, d(2:2:end,:,:,:,:,:) = -d(2:2:end,:,:,:,:,:); end
d = mean(d,1);
fprintf('TG=%g; R1=%g; R2=%g; cv14=%g\n',...
    h.rdb_hdr.ps_mps_tg,h.rdb_hdr.ps_mps_r1,h.rdb_hdr.ps_mps_r2,...
    h.image.user14);


%% SNR calculation
if ncoils>1         % individual for getting the weights
    fprintf(outid,'slice coil signal noise SNR corSNR\n');
    for ls=1:nslices
        for lc=1:ncoils
            % determine SNR
            signal = mean(abs(d(:,ind_signal,1,1,ls,lc)),2);
            noise = std(d(:,ind_noise,1,1,ls,lc));
            snr(ls,lc) = signal/noise;
            fprintf(outid,'%g\t%g\t%g\t%g\t%g\t',ls,lc,signal,noise,snr(ls,lc));
            % SNR corrected: 1 exc
            snr_corr(ls,lc) = snr(ls,lc)/sqrt(nexc);
            fprintf(outid,'%g\n',snr_corr(ls,lc));
        end
    end
end


%% SNR optimal coil combination
if ncoils>1
    % ww = sum(snr_corr,1); ww = ww/sum(ww)
    ww = snr_corr./(sum(snr_corr,2)*ones(1,ncoils));  % coil weighting
    pc = mean(unwrap(angle(d(1,ind_pts,1,1,:,:))),2);
    pc = reshape(pc,[nslices ncoils]);   % phase correction
    for ls=1:nslices
        for lc=1:ncoils
            d(1,:,1,1,ls,lc) = ww(ls,lc)*d(1,:,1,1,ls,lc)*exp(-1i*pc(ls,lc));
        end
    end
    d = sum(d,6);
end
d = reshape(d,[1 npts nslices]);

snrCMB = NaN(1,nslices);
% snr_corrCMB = NaN(1,nslices);

if ncoils>1, fprintf('\nSNR with SNR-optimal coil combination\n'); end


%% loop through slices
for ls=1:nslices
    signal = mean(abs(d(:,ind_signal,ls)),2);
    noise = std(d(:,ind_noise,ls));
    snrCMB(ls) = signal/noise;
    fprintf(outid,'signal=%g noise=%g SNR=%g\n',signal,noise,snrCMB(1,ls));
end
