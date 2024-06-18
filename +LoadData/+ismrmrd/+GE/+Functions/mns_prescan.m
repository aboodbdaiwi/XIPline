function [tg,f0,snr,lb] = mns_prescan(d,h,fn_log,plt,fname,do_sort)
%MNS_PRESCAN  Determine TG (Bloch-Siegert), SNR, f0, lb; acquired w/ fidall
% [tg,f0,snr,lw] = mns_prescan(d,h,fn_log,plt,fname,do_sort)
%
%         d  Raw (p-file) data  (or pfile fname)
%         h  Header from p-file (or empty)
%    fn_log  Save summary into this log-file                           ([])
%       plt  Plotting                                                (true)
%     fname  Print pix + write function output to file ([]->stdout)    ([])
%   do_sort  Sort phase difference according to |FID|                (true)
%
%        tg  Transmit gain
%        f0  Centre frequency [Hz]
%       snr  Normalised SNR (1 excitation & 10mm slice)
%        lw  Linewidth [Hz]
% 
% Literature: NMRiB DOI: 10.1002/nbm.1657
% 2/2021 Rolf Schulte
if (nargin<1), help(mfilename); return; end


%% input defaults
if ~exist('fn_log','var'), fn_log = []; end
if ~exist('fname','var'),  fname = []; end
if ~exist('plt','var'),    plt = []; end
if isempty(plt),           plt = true; end
if ~exist('do_sort','var'),do_sort = []; end
if isempty(do_sort),       do_sort = true; end
if isempty(fname)
    outid = 1;                                   % direct to stdout
    msgbx = false;
else
    if ischar(fname)
        if ~isempty(regexpi(fname,'\.7$')),   fname = fname(1:end-2); end
        if ~isempty(regexpi(fname,'\.h5$')),  fname = fname(1:end-3); end
        if ~isempty(regexpi(fname,'\.mat$')), fname = fname(1:end-4); end
        if ~isempty(regexpi(fname,'\.png$')), fname = fname(1:end-4); end
        if ~isempty(regexpi(fname,'\.txt$')), fname = fname(1:end-4); end
        outid = fopen([fname '.txt'],'w+');      % direct to file
        msgbx = true;                            % display message box
    else
        fname = [];
        outid = 1;
        msgbx = false;
    end
end
export_dcm = true;                     % export prescan window to dicom database


%% reading in data, if pfile name is given
if ~isnumeric(d)
    if exist(d,'file')
        [d,h] = read_p(d);
    else
        msg = sprintf('Warning: strange input d');
        fprintf(outid,'%s\n',msg);
        if msgbx, warndlg(msg,'mns_prescan'); end
    end
end
if size(d,3)~=1, warning('size(d,3)(=%d)~=1',size(d,3)); end
if size(d,4)~=1, warning('size(d,4)(=%d)~=1',size(d,4)); end


%% default parameters
n_pts = 30;             % #FID-points used for TG calculation
ind_signal = (3:10);    % points of FID used for signal (area)
pts_noise = 200;        % #points at end of FID for noise
f0_current = h.rdb_hdr.ps_mps_freq/10;
slthick = h.image.slthick;
iwidth = 10;             % width of centre-of-mass for freq determination

if ~isempty(fn_log)
    if ~exist(fn_log,'file')
        fileid = fopen(fn_log,'a');
        fprintf(fileid,'date time slice nucleus phase std ');
        fprintf(fileid,'TG(curr change set) ');
        fprintf(fileid,'f0(curr change set) lb signal nexc SNR pulse thick\n');
    else
        fileid = fopen(fn_log,'a');
    end
end


%% check data
if any(isinf(d(:)))
    warning('isinf(d(:)); setting to 0'); 
    d(isinf(d)) = 0; 
end
if any(isnan(d(:)))
    warning('isnan(d(:)); setting to 0'); 
    d(isnan(d)) = 0;
end
if isempty(regexpi(h.image.psd_iname,'fidall'))
    msg = sprintf('Warning: psd_iname = %s; wrong psd!? (use fidall)', ...
        h.image.psd_iname);
    fprintf(outid,'%s\n',msg);
    if msgbx, warndlg(msg,'mns_prescan'); end
end
[nexc,npts,~,~,nslices,ncoils] = size(d);
n_fft = max([8192 2*npts]);            % interpolate to n_fft
nexc_noise = h.rdb_hdr.user32;
if nexc_noise>0
    fprintf('Noise scans = %d\n',nexc_noise);
    if nexc_noise>nexc
        error('nexc_noise(=%d)>nexc(=%d)',nexc_noise,nexc); 
    end
    nexc = nexc-nexc_noise;
    if nexc_noise==1
        % last data frame
        ind1_noise = nexc+nexc_noise;
    else
        % skip first noise scan to avoid remaining signal
        ind1_noise = nexc+(2:nexc_noise);
    end
    ind2_noise = [];
else
    if npts<pts_noise
        warning('npts(=%g) < pts_noise(=%g); setting pts_noise = %g',...
            npts,pts_noise,npts-10);
        pts_noise = npts-10;
    end
    ind2_noise = ((npts-pts_noise):npts);
end
if (abs((nexc/2)-floor(nexc/2))>1d-10)
    msg = sprintf('Error: nexc(=%g) not multiple of 2',nexc);
    fprintf(outid,'%s\n',msg);
    if msgbx, hid = errordlg(msg,'mns_prescan'); uiwait(hid); end
    if ~isempty(fname),  fclose(outid); end
    if ~isempty(fn_log), fclose(fileid); end
    error(msg);
end


%% misc variables
B1max = h.rdb_hdr.user34*100; % B1max [uT] of exc pulse (scaled to pw & flip)
nucleus = h.rdb_hdr.user2;             % nucleus
if B1max<0.01
    msg = sprintf('Warning: B1max (=%g)<0.01 [uT]; rhuser34=%g', ...
        B1max,h.rdb_hdr.user34);
    fprintf(outid,'%s\n',msg);
    if msgbx, warndlg(msg,'mns_prescan'); end
    B1max = 12.8079*gyrogamma(1)/abs(gyrogamma(nucleus)); 
    % cv14=1, flip_rf1=90°(pw_rf1=1800)
end

chp = ischop(h);                       % rf chopping
bw = h.rdb_hdr.user0;                  % (full) BW [Hz]
t = (1:size(d,2))/bw*1d3;              % time-axis [ms]

off_freq = h.rdb_hdr.user30;           % Bloch-Siegert pulse frequency [Hz]
off_pw = h.rdb_hdr.user31;             % pw blosi pulse [us]
nbp = h.rdb_hdr.user19;                % # of blosi off-freq pulses/sign
B1_desired = h.rdb_hdr.user33*100;     % B1 of Bloch-Siegert pulse [uT]
KBS = 4109.6;      % deg/G^2; conversion for Bloch_Siegert.rho (8ms,4kHz)

if nbp<1
    msg = sprintf('Error: blosi_npulses(=%g) < 1',nbp);
    fprintf(outid,'%s\n',msg);
    if msgbx, hid = errordlg(msg,'mns_prescan'); uiwait(hid); end
    if ~isempty(fname),  fclose(outid); end
    if ~isempty(fn_log), fclose(fileid); end
    error(msg);
end

if 2*nbp>nexc
    msg = sprintf('Error: 2*blosi_npulses(=%g) > nexc(=%g)',nbp,nexc);
    fprintf(outid,'%s\n',msg);
    if msgbx, hid = errordlg(msg,'mns_prescan'); uiwait(hid); end
    if ~isempty(fname),  fclose(outid); end
    if ~isempty(fn_log), fclose(fileid); end
    error(msg);
end


%% phase-shift constant
% KBS -> deg/G^2 Bloch-Siegert shift conversion from phase to B1 magnitude
KBS = KBS/100^2;             % convert Gauss into uT
KBS = KBS/(off_freq/4000);   % scale for frequencies other than 4kHz
KBS = KBS*(off_pw/8000);     % scale to pw's other than 8ms
scale = (gyrogamma(nucleus)/gyrogamma('1h'))^2;
KBS = scale*KBS;             % scale to nucleus


%% remove RF chopping
if chp, d(2:2:end,:,:,:,:,:) = -d(2:2:end,:,:,:,:,:); end


%% sort acquisitions with positive/negative frequency
% to opposite order of pos,neg vs. neg,pos 4 kHz, if B1 comes out negative
% or imaginary, just switch 'in' and 'ip' below
dip = zeros(1,npts,1,1,nslices,ncoils); 
din = zeros(1,npts,1,1,nslices,ncoils);
for l=1:nbp
    dip = dip + mean(d(l:2*nbp:nexc,:,1,1,:,:),1);        % +off, averaged FID
    din = din + mean(d((l+nbp):2*nbp:nexc,:,1,1,:,:),1);  % -off, averaged FID
end
dip = dip/nbp; din = din/nbp;


%% generate index list for phase-difference
if do_sort
    fid = sqrt(mean(mean(abs(conj(dip).*dip+conj(din).*din),5),6));
    [~,ind_pts] = sort(fid,'descend');
    ind_pts = sort(ind_pts(1,1:n_pts));
else
    ind_pts = (1:n_pts)+2;
end


%% multi-channel SNR calculation
if ncoils>1  % individual for getting the weights
    % calculate noise
    noise_coils = zeros(1,ncoils);
    nnn = [];
    for lc=1:ncoils
        if nexc_noise>0
            tmp = d(ind1_noise,:,:,:,:,lc);
        else
            tmp = d(:,ind2_noise,:,:,:,lc);
        end
        noise_coils(1,lc) = std(tmp(:));
        nnn = [nnn tmp(:)];
    end
    % calculate signal
    dtmp = (abs(dip)+abs(din))/2;
    signal_coils = reshape(mean(dtmp(:,ind_signal,1,1,:,:),2),[nslices ncoils]);
    snr_coils = bsxfun(@rdivide,signal_coils,noise_coils);
    
    % calculate noise correlation and covariance
    sub_noise_correlation(nnn,outid,fname,plt);
end


%% multi-channel svds coil combination
if ncoils>1
    [din,dip,d] = sub_coil_combine(din,dip,d,n_pts);
end
din = reshape(din,[1 npts nslices]);
dip = reshape(dip,[1 npts nslices]);


%% calculation of phase difference
B1phase = angle(conj(din).*dip)*180/pi;
tg = NaN(1,nslices); f0_change = tg; f0 = tg; lb = tg; snr = tg;


%% loop through slices
for ls=1:nslices
    % plotting
    if nslices>1
        fprintf(outid,'\nSlice = %g\n',ls);
        if plt, figure(ls); clf; end
    else
        figure;
    end
    if plt
        subplot(3,1,1);
        plot(t,abs(din(:,:,ls)).','r',t,abs(dip(:,:,ls)).','b',...
            t(ind_signal),abs(din(:,ind_signal,ls)).','mx',...
            t(ind2_noise), abs(din(:,ind2_noise,ls)).','mx',...
            t(ind_signal),abs(dip(:,ind_signal,ls)).','cx',...
            t(ind2_noise), abs(dip(:,ind2_noise,ls)).','cx');
        ylabel('FID (abs)'); axis tight;
        grid on
        
        subplot(3,1,2);
        plot(t,B1phase(1,:,ls),'m',(ind_pts)/bw*1d3,B1phase(1,ind_pts,ls),'gx');
        xlabel('time [ms]'); ylabel('Phase diff [deg]'); 
        axis([1 max(t) -220 220]);
        grid on
    end
    
    % calculte noise
    if nexc_noise>0
        tmp = d(ind1_noise,:,:,:,:);
    else
        tmp = d(:,ind2_noise,:,:,:);
    end
    noise = std(tmp(:));
    
    % calculate signal
    dtmp = (abs(din)+abs(dip))/2;
    signal = mean(dtmp(:,ind_signal,ls),2);
    
    % determine SNR
    snr(1,ls) = signal/noise;
    fprintf(outid,'SNR = %g/%g = %g\n',signal,noise,snr(1,ls));
    if ncoils>1
        fprintf(outid,'SNR coils = %s\n',num2str(snr_coils(ls,:),'%.4g '));
    end
    
    % calculate TG values
    mb1p = mean(B1phase(1,ind_pts,ls));
    sb1p = std(B1phase(1,ind_pts,ls));
    if sb1p>abs(mb1p)
        msg = sprintf('Warning: std(phase)(=%g) > abs(mean(phase))(=%g)',sb1p,abs(mb1p));
        fprintf(outid,'%s\n',msg);
        if msgbx, warndlg(msg,'mns_prescan'); end
    end
    if mb1p>150
        msg = sprintf('Warning: mean(phase)(=%g) > 150',mb1p);
        fprintf(outid,'%s\n',msg);
        if msgbx, warndlg(msg,'mns_prescan'); end
    end
    if mb1p<10
        msg = sprintf('Warning: mean(phase)(=%g) < 10',mb1p);
        fprintf(outid,'%s\n',msg);
        if msgbx, warndlg(msg,'mns_prescan'); end
    end
    b1 = sqrt(mb1p/2/KBS);
    b1err1 = sqrt((mb1p+sb1p)/2/KBS);
    b1err2 = sqrt((mb1p-sb1p)/2/KBS);
    tg_change = -10*20*log10(b1/B1_desired);
    tg_change_err1 = -10*20*log10(b1err1/B1_desired);
    tg_change_err2 = -10*20*log10(b1err2/B1_desired);
    tg_current = h.rdb_hdr.ps_mps_tg;
    tg(ls) = tg_current + tg_change;
    fprintf(outid,'Blosi phase = %.4g +/- %.4g [deg]\n',mb1p,sb1p);
    fprintf(outid,'Nucleus = %g; pulse=%g; thick=%g\n', ...
        nucleus,h.image.user14,slthick);
    fprintf(outid,'B1: desired = %.4g [uT]; actual = %.4g (%.4g...%.4g)[uT]\n',...
        B1_desired,b1,b1err2,b1err1);
    fprintf(outid,'TG: current = %g; change = %g (%g...%g); setting = %g\n',...
        tg_current,round(tg_change),round(tg_change_err1),...
        round(tg_change_err2),round(tg(ls)));
    if plt
        sign_str = '+'; if tg_change<0, sign_str = ''; end
        title_str = sprintf('TG=%g (%g%c%g); phi=%.4g+-%.3g', ...
            round(tg(ls)),tg_current,sign_str,round(tg_change),...
            mb1p,sb1p);
        if round(tg(ls))>200
            title_str = sprintf('%s; FA*=%.3g (for TG=200)',...
                title_str,10^((tg(ls)-200)/200));
        end
        title(title_str);
    end

    % spectra; Fourier interpolated for higher precision
    hz = -(-n_fft/2:n_fft/2-1)/n_fft*bw;    % freq axis
    ssn = fftshift(fft(conj(din(1,:,ls)),n_fft,2),2);
    ssp = fftshift(fft(conj(dip(1,:,ls)),n_fft,2),2);
    
    % centre of mass f0 determination
    [~,imax] = max(abs(ssn));
    ind = (imax-iwidth):(imax+iwidth);
    if any(ind<1)
        warning('ind<1: truncating');
        ind = ind(ind>0);
    end
    if any(ind>n_fft)
        warning('ind>n_fft(=%d): truncating',n_fft);
        ind = ind(ind<=n_fft);
    end
    f0n = sum(hz(ind).*abs(ssn(1,ind)))/sum(abs(ssn(1,ind)));
    [~,imax] = max(abs(ssp));
    ind = (imax-iwidth):(imax+iwidth);
    if any(ind<1)
        warning('ind<1: truncating');
        ind = ind(ind>0);
    end
    if any(ind>n_fft)
        warning('ind>n_fft(=%d): truncating',n_fft);
        ind = ind(ind<=n_fft);
    end
    f0p = sum(hz(ind).*abs(ssp(1,ind)))/sum(abs(ssp(1,ind)));
    
    % fprintf('f0n=%g; f0p=%g\n',f0n,f0p);
    f0_change(ls) = -(f0n+f0p)/2;    
    f0(ls) = round(f0_current + f0_change(ls));
    fprintf(outid,'f0: current = %d; change = %d; setting = %d [Hz]\n',...
        f0_current,round(f0_change(ls)),f0(ls));
    
    % determine linewidths (full width at half maximum)
    [lwn,outn] = fwhm(ssn,bw,[],false,false);
    [lwp,outp] = fwhm(ssp,bw,[],false,false);
    % fprintf('lwn=%g; lwp=%g\n',lwp,lwn);
    lb(ls) = (lwn+lwp)/2;

    fprintf(outid,'fwhm = %.4g [Hz]\n',lb(ls));
    
    % plot spectra
    if plt
        subplot(3,1,3);
        plot(hz,abs(ssn),'r', hz,abs(ssp),'b', ...
            f0_change(ls)*[1 1],[-0.02 1.2]*max(abs([ssn ssp])),'g',...
            outn.freqs(1,:),abs(ssn(1,outn.ind_max(1,:))),'m*',...
            [outn.fp(1,:) outn.fm(1,:)],...
            repmat(abs(ssn(1,outn.ind_max(1,:)))/2,[2 1]),'m*-',...
            outp.freqs(1,:),abs(ssp(1,outp.ind_max(1,:))),'c*',...
            [outp.fp(1,:) outp.fm(1,:)],...
            repmat(abs(ssp(1,outp.ind_max(1,:)))/2,[2 1]),'c*-');
        xlabel('freq [Hz]'); ylabel('Spectrum [a.u.]'); axis tight;
        legend('-\omega','+\omega');
        grid on
        
        df = f0(ls)-f0_current;
        if df<0
            sign_str = '-';
        else
            sign_str = '+';
        end
        title(sprintf('f0=%d (%c%g)[Hz]; lb=%.4g [Hz]', ...
            f0(ls),sign_str,abs(df),lb(ls)));
    end
    
    if plt
        subplot(3,1,1); 
        title(sprintf('SNR=%.3g (1 exc); pulse=%g; thick=%g [mm]', ...
            snr(1,ls),h.image.user14,slthick));
        
        figstr = sprintf('P%05d Exam%d Series%d',...
            h.image.rawrunnum,h.exam.ex_no,h.series.se_no);
        if nslices>1, figstr = sprintf('%s Slice%d',figstr,ls); end
        set(gcf,'name',figstr);

        % print figures
        if ~isempty(fname)
            if isdeployed
                renderer = '-painters';
            else
                renderer = '-opengl';
            end
            print('-dpng',[fname '_' num2str(ls) '.png'],'-r600',renderer);
            if export_dcm
                inp.InstanceNumber = ls;
                write_scdicom([fname '_' num2str(ls) '.dcm'],gcf,h,inp);
            end
        end
    end
    
    % store results into file
    if ~isempty(fn_log) && (snr(1,ls)>0.1)
        fprintf(fileid,'%s %s %d %g %.4g %.4g ',...
            h.rdb_hdr.scan_date,h.rdb_hdr.scan_time,ls,nucleus,mb1p,sb1p);
        fprintf(fileid,'%g %g %g ',...
            tg_current,round(tg_change),round(tg(ls)));
        fprintf(fileid,'%d %d %d %.4g ',...
            f0_current,round(f0_change(ls)),round(f0(ls)),lb(ls));
        fprintf(fileid,'%.4g %d %.4g %g %.2g\n',...
            signal,nexc,snr(1,ls),h.image.user14,slthick);
    end
end


%% close output files
if ~isempty(fname),  fclose(outid); end
if ~isempty(fn_log), fclose(fileid); end


end      % mns_prescan.m


%% sub-functions
% svds coil combination
function [din,dip,dd] = sub_coil_combine(din,dip,d,pts4svds)

%% miscelaneous parameters
[~,ns,~,~,n5,nc] = size(d);

%% coil combination
dtmp = sqrt(sum(sum((din.*conj(din)+dip.*conj(dip)),5),6));
[~,ind] = sort(dtmp);
if pts4svds<ns
    ind_svds = sort(ind((ns-pts4svds+1):ns));
else
    ind_svds = (1:ns);
end
ww = zeros(1,1,1,1,n5,nc);

for l5=1:n5
    A = squeeze(din(1,ind_svds,1,1,l5,:));
    [~,~,V,flag] = svds(A,1);
    ww(1,1,1,1,l5,:) = V;
    if flag~=0, fprintf('l5=%g; flag=%g\n',l5,flag); end
end
ww = bsxfun(@rdivide,ww,sum(abs(ww),6));
dd = sum(bsxfun(@times,d,ww),6);
din = sum(bsxfun(@times,din,ww),6);
dip = sum(bsxfun(@times,dip,ww),6);


end      % sub_coil_combine


%% calculate noise correlation and covariance
function sub_noise_correlation(nnn,outid,fname,plt)

ncoils = size(nnn,2);
if isdeployed
    renderer = '-painters';
else
    renderer = '-opengl';
end


%% calculate noise correlation and covariance
noisecor = abs(corrcoef(nnn));               % noise correlation
noisecov = abs(nnn'*nnn);                    % noise covariance
noisecov = noisecov/mean(diag(noisecov));


%% plotting
if plt
    if isempty(fname)
        figure(101); clf
    else
        fid = figure('Visible','off');
    end
    imagesc(noisecor,[0 1]); colorbar;
    title('Noise Correlation')
    % print figures
    if ~isempty(fname)
        print('-dpng',[fname '_noisecor.png'],renderer);
        close(fid);
    end
    
    if isempty(fname)
        figure(102); clf
    else
        fid = figure('Visible','off');
    end
    mm = (2+max(noisecov(:)))/2;
    imagesc(noisecov,[0 mm]); colorbar;
    title('Noise Covariance')
    % print figures
    if ~isempty(fname)
        print('-dpng',[fname '_noisecov.png'],renderer);
        close(fid);
    end
end


%% print info
fprintf(outid,'noisecor = \n');
tmp = num2str(noisecor,'%.3g  ');
for lc=1:ncoils, fprintf(outid,'%s\n',tmp(lc,:)); end
fprintf(outid,'noisecov = \n');
tmp = num2str(noisecov,'%.3g  ');
for lc=1:ncoils, fprintf(outid,'%s\n',tmp(lc,:)); end

fprintf(outid,'\nnoise correlation: max=%.2g%%; min=%.2g%%; mean=%.2g%%\n',...
    max(max(noisecor-eye(ncoils)))*100,min(min(noisecor))*100,...
    (sum(sum(noisecor))-ncoils)/(ncoils*(ncoils-1))*100);

fprintf(outid,'noise covariance:\n');
fprintf(outid,'\ton-diagonal:  max=%.3g; min=%.3g; mean=%.3g; std=%.3g\n',...
    max(diag(noisecov)),min(diag(noisecov)),...
    mean(diag(noisecov)),std(diag(noisecov)));
tmp = noisecov + diag(inf(ncoils,1));
tmp = tmp(~isinf(tmp));
fprintf(outid,'\toff-diagonal: max=%.3g; min=%.3g; mean=%.3g; std=%.3g\n\n',...
    max(tmp),min(tmp),mean(tmp),std(tmp));


end      % sub_noise_correlation
