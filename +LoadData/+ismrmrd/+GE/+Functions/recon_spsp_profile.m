function [b1,f_ax,z_ax] = recon_spsp_profile(d,h,wf_name,spsp_df,pix_fname,msgbx)
%RECON_SPSP_PROFILE  Reconstruct 2D pulse profile
%   Use Cartesian readout (eg gz_fov200_mtx512_1h_gmax23_smax75 or 
%   gz_fov200_mtx512_13c_gmax23_smax75) with an NMR test tube aligned 
%   along z, and step through excitation frequencies
%   Also reconstructs 2D spatial, if csi phase encodes on
%
% [b1,f_ax,z_ax] = recon_spsp_profile(d,h,wf_name,spsp_df,pix_fname,msgbx)
%
%         d  Raw (p-file) data  (or pfile fname)
%         h  Header from p-file (or empty)
%   wf_name  Location of waveform .mat file
%   spsp_df  df of frequency axis  (opt; default=h.rdb_hdr.user48) [Hz]
% pix_fname  Print pix into this file as png (opt,[]); 
%     msgbx  Display output in message box (opt,default=false)
%
%        b1  Measured profile
%      f_ax  Frequency axis
%      z_ax  Spatial axis
%
% 7/2020  Rolf Schulte
if (nargin<1), help(mfilename); return; end


%% input variables
if ~exist('spsp_df','var')
    spsp_df = []; 
else
    if ischar(spsp_df), spsp_df = str2double(spsp_df); end
end

if ~exist('pix_fname','var'), pix_fname = []; end
if ~isempty(pix_fname)
    if ~isempty(regexpi(pix_fname,'\.png'))
        pix_fname = pix_fname(1:end-4);
    end
    if ~isempty(regexpi(pix_fname,'\.7'))
        pix_fname = pix_fname(1:end-2);
    end
end
if ~exist('msgbx','var'), msgbx = false; end
if ischar(msgbx)
    msgbx = str2num(msgbx); 
    if isempty(msgbx), msgbx = false; end
end


%% reading in data, if pfile name is given
if ~isnumeric(d)
    if exist(d,'file')
        [d,h] = read_p(d);
    else
        if msgbx
            hid = warndlg('strange input d','recon_spsp_profile');
            uiwait(hid);
        end
        warning('recon_spsp_profile:read_p','strange input d');
    end
end


%% parameters + checks
if isempty(spsp_df) || isnan(spsp_df)
    spsp_df = h.rdb_hdr.user48;
end

wf = load(wf_name);          % load waveform
n1 = size(d,1);
n2 = size(d,2);

if n1<2
    msg = sprintf('size(d,1)=%g\n',n1);
    if msgbx, hid = warndlg(msg,'recon_spsp_profile'); uiwait(hid); end
    warning('recon_spsp_profile:n1',msg); 
end
if n2<2
    msg = sprintf('Warning: size(d,2)=%g\n',n2);
    if msgbx, hid = warndlg(msg,'recon_spsp_profile'); uiwait(hid); end
    warning('recon_spsp_profile:n2',msg);
end
d = squeeze(d);

chp = ischop(h);                                 % rf chopping
if chp, d(2:2:end,:,:) = -d(2:2:end,:,:); end    % remove RF chopping

spsp_nfreqs = h.rdb_hdr.user22;
if spsp_nfreqs>0
    if isodd(spsp_nfreqs)
        msg = sprintf('spsp_nfreqs (=%g) is odd',spsp_nfreqs);
        if msgbx, hid = errordlg(msg,'recon_spsp_profile'); uiwait(hid); end
        error(msg);
    end
    if spsp_df == 0
        warning('recon_spsp_profile:spsp_df','spsp_df == 0 -> setting to 4'); 
        spsp_df = 4;
    end
    f_ax = (-spsp_nfreqs/2:spsp_nfreqs/2-1)*spsp_df; % freq axis
    xlb = 'freq [Hz]';
else
    f_ax = 1:size(d,1);
    spsp_nfreqs = size(d,1);
    xlb = 'acqs';
end
mtx = wf.in.mtx_res;                        % matrix size
z_ax = (-mtx/2:mtx/2-1)/mtx*h.rdb_hdr.fov;  % z axis [mm]


%% reconstruction
d1 = d(:,wf.ind,:);                         % apply index list
if n1<spsp_nfreqs, f_ax = f_ax(1:n1); end   % truncate f-axis
if n1>spsp_nfreqs                           % average signal
    fac = n1/spsp_nfreqs;
    if abs(fac-floor(fac))>1d-15
        msg = sprintf('size(d,1) (=%g) not multiple of spsp_nfreqs (=%g)\n',...
            n1,spsp_nfreqs);
        if msgbx, hid = errordlg(msg,'recon_spsp_profile'); uiwait(hid); end
        error(msg);
    end
    tmp = zeros(spsp_nfreqs,mtx,size(d1,3));
    for l=1:fac, tmp = tmp+d1(((1:spsp_nfreqs)+(l-1)*spsp_nfreqs),:,:); end
    d1 = tmp;
end
b1 = fftshift(ifft(ifftshift(d1,2),[],2),2);     % actual reconstruction
title_str = sprintf('cv14=%g',h.image.user14);   % title string


%% recon additional phase encoding dimension
csires = max([h.image.user5 h.image.user6 h.image.user7]);
if csires>1
    fprintf('CSI detected: cv5(rl)=%g cv6(ap)=%g cv7(si)=%g\n',...
        h.image.user5,h.image.user6,h.image.user7);
    fprintf('FFT also along dim1\n');
    
    b1 = fftshift(ifft(ifftshift(b1,1),[],1),1);
    f_ax = (-csires/2:csires/2-1)/csires*h.rdb_hdr.fov;  % x axis [mm]
    xlb = 'x [mm]';
end


%% plotting
n3 = size(d,3);
clf
set(gcf,'DefaultAxesFontSize',12)
cmax = 0.9*max(abs(b1(:)));
figstr = sprintf('P%05d Exam%d Series%d',...
    h.image.rawrunnum,h.exam.ex_no,h.series.se_no);
set(gcf,'name',figstr);

for l3=1:n3
    if n3>1, subplot(n3,1,l3); end
    imagesc(f_ax,z_ax,abs(b1(:,:,l3)).',[0 cmax]);
    pos = get(gca,'Position');
    set(gca,'Position',[pos(1)-0.02 pos(2)+0.01 1.1*pos(3) 1.2*pos(4)]);
    if n3==l3
        xlabel(xlb);
        ylabel('z [mm]');
    else
        set(gca,'XTickLabel','');
    end
    if l3==1, title(title_str); end
end


%% saving pix + data
if ~isempty(pix_fname)
    print([pix_fname '.png'],'-dpng','-r600','-painters');
    save([pix_fname '.mat'],'-mat','b1','f_ax','z_ax','wf_name',...
        'spsp_df','h','figstr');
end

end   % main function recon_spsp_profile.m
