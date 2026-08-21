function eval_pol(d_h,h_h,d_t,h_t,t_diss_scan,ind_t1)
%EVAL_POL  Evaluate polarisation level for GE FIDCSI sequence
%  FIDCSI (eg hardpulse; sw=5000; npt=2048; dda=0; nex=1; TG const.)
%  Hyperpol FID: reduced R1/R2, flip_rf1=1.4; pw_rf1=8; optr=5s; cv4=32
%  Thermal  FID: R1=13;R2=30; flip_rf1=90; pw_rf1=512; optr=2s; cv4=512
%                doping to reduce T1: ~3:100 DOTAREM (or alike)
%                (e.g., 2ml PA + ~60mg DOTAREM)
%  Signal level obtained by averaging first points in FID
%   -> Requires one predominant resonance:
%   -> Pyruvate solution should be about pH neutral
%
%    eval_pol(d_h,h_h,d_t,h_t,t_diss_scan,ind_t1)
%         d_h   Raw data or file name hyperpolarised signal
%         h_h   Header hyperpolarised signal
%         d_t   Raw data or file name thermally polarised signal
%         h_t   Header thermal signal
% t_diss_scan   Time from dissolution to acquisition [s]
%      ind_t1   Index of last (or first+last) acq for fitting T1
%
% 11/2022  Rolf Schulte
%See also READ_MR_RAWDATA.
if (nargin<3), help(mfilename); return; end


%% reading in data, if pfile name is given
if ~isnumeric(d_h)
    if ~exist(d_h,'file') && ~isempty(d_h)
        warning('strange input for d (=''%s''); file not existing?',d); 
    end
    if ~isempty(h_h), warning('~isempty(h_h)'); end
    if strcmpi(d_h(end-3:end),'.mat')
        tmp = load(d_h); d_h = tmp.d; h_h = tmp.h; clear tmp;
    else
        [d_h,h_h] = read_p(d_h);
    end
end
if ~exist('h_t','var'), h_t = []; end
if ~isnumeric(d_t)
    if ~exist(d_t,'file') && ~isempty(d_t)
        warning('strange input for d (=''%s''); file not existing?',d); 
    end
    if ~isempty(h_t), warning('~isempty(h_t)'); end
    if strcmpi(d_t(end-3:end),'.mat')
        tmp = load(d_t); d_t = tmp.d; h_t = tmp.h; clear tmp;
    else
        [d_t,h_t] = read_p(d_t);
    end
end
if length(size(d_h))~=2, error('d_h matrix must be 2D'); end
if length(size(d_t))~=2, error('d_t matrix must be 2D'); end


%% get parameters from header or input
if ~exist('t_diss_scan','var'), t_diss_scan = []; end
if ~exist('ind_t1','var'),      ind_t1 = []; end
switch length(ind_t1)
    case 0, ind_t1 = [1 size(d_h,1)];
    case 1, ind_t1 = [1 ind_t1];
    case 2
    otherwise
        error('length(ind_t1)(=%d)~=1',length(ind_t1));
end
if ind_t1(1)<1, error('ind_t1(1)<1'); end
if ind_t1(2)>size(d_h,1), error('ind_t1(2)>size(d_h,1)'); end
if any(abs(ind_t1-floor(ind_t1))>1d-15), error('ind_t1 not integer'); end

indave = (3:15);                  % index for averaging first point in FID
flip_rf1_h = h_h.image.mr_flip;   % flip angles [deg]
flip_rf1_t = h_t.image.mr_flip;
r1_h = h_h.rdb_hdr.ps_mps_r1;     % analog receiver gain
r1_t = h_t.rdb_hdr.ps_mps_r1;
r2_h = h_h.rdb_hdr.ps_mps_r2;     % digital receiver gain
r2_t = h_t.rdb_hdr.ps_mps_r2;
tg_h = h_h.rdb_hdr.ps_mps_tg;     % transmit gain [10TG=1dB]
tg_t = h_t.rdb_hdr.ps_mps_tg;
f0 = h_t.rdb_hdr.ps_mps_freq/10;  % Larmor frequency [Hz]
tr_h = h_h.image.tr*1d-6;         % repetition time
tr_t = h_t.image.tr*1d-6;
t_h = (0:h_h.rdb_hdr.user1-1)/h_h.rdb_hdr.user0;  % acq time
t_t = (0:h_t.rdb_hdr.user1-1)/h_t.rdb_hdr.user0;  % acq time


%% print info
fprintf('\n*************************\nParameters\n');
fprintf('   hyperpol: R1=%g, R2=%g, TG=%g, flip_rf1=%g [deg], TR=%g [s]\n',...
    r1_h,r2_h,tg_h,flip_rf1_h,tr_h);
fprintf('\t\tnucleus = %g; pulse=%g; thick=%g\n', ...
    h_h.rdb_hdr.user2,h_h.image.user14,h_h.image.slthick);
fprintf('    thermal: R1=%g, R2=%g, TG=%g, flip_rf1=%g [deg], TR=%g [s]\n',...
    r1_t,r2_t,tg_t,flip_rf1_t,tr_t);
fprintf('\t\tnucleus = %g; pulse=%g; thick=%g\n', ...
    h_t.rdb_hdr.user2,h_t.image.user14,h_t.image.slthick);

if (tg_h ~= tg_t), warning('different TGs'); end


%% reconstructions
[s_h,~,hz_h] = sub_echo2spec(d_h,h_h);
sig_h = max(mean(abs(d_h(:,indave)),2));

[s_t,dave_t,hz_t] = sub_echo2spec(d_t,h_t);
sig_t = mean(abs(dave_t(indave)));


%% scalings
scale_rg = 2^((r1_t - r1_h)/2+(r2_t - r2_h));     % scale for receiver gain
scale_flip = sin(flip_rf1_t*pi/180)/sin(flip_rf1_h*pi/180); % scale for FA


%% pysical parameters
hq = 1.05457d-34;                 % [Js]
T = 313;                          % [K]
kB = 1.380658d-23;                % [J/K]


%% calculate polarisation levels
pol_t = 2*pi*hq*f0/(2*kB*T);                    % thermal polarisation
pol_h = sig_h/sig_t*pol_t*scale_flip*scale_rg;  % hyper-polarisation


%% fit to extract T1 from hyperpolarised data
tt = ((ind_t1(1):ind_t1(2))-1)*tr_h;
t1decay = mean(abs(d_h((ind_t1(1):ind_t1(2)),indave)),2).';

p = polyfit(tt,log(t1decay),1);
t1fit = polyval(p,tt);
t1 = -1/p(1);


%% plotting 
subplot(3,2,1); 
plot(t_h*1d3,abs(d_h(1,:)),'b',t_h(indave)*1d3,abs(d_h(1,indave)),'g-x'); 
axis tight; ylabel('hyper pol'); xlabel('FID [ms]'); grid on;
subplot(3,2,2); 
plot(hz_h,abs(s_h).'); axis tight; xlabel('spectrum'); grid on;
set(gca,'XDir','reverse');   % spectroscopy convention
subplot(3,2,3); 
plot(t_t*1d3,abs(dave_t),'b',t_t(indave)*1d3,abs(dave_t(indave)),'g-x'); 
axis tight; ylabel('thermal pol'); xlabel('FID [ms]'); grid on;
subplot(3,2,4); 
plot(hz_t,abs(mean(s_t,1))); axis tight; xlabel('spectrum'); grid on;
set(gca,'XDir','reverse');   % spectroscopy convention
subplot(3,2,5); 
plot(tt,t1decay,'b-x',tt,exp(t1fit),'r-x'); 
axis tight; xlabel('T1 decay [s]'); grid on;
subplot(3,2,6); 
plot(tt,log(t1decay),'b-x',tt,t1fit,'r-x'); 
axis tight; xlabel('T1 fit (log) [s]'); grid on;


%% print output
fprintf('\n***********\n');
fprintf('T1 = %.3g [s]\n',t1);
fprintf('Thermal polarisation = %.3g [ppm]\n',pol_t*1d6);
fprintf('Hyper-polarisation\n');
fprintf('\tuncorrected:      \t%.3g [%%]\t\n',pol_h*1d2);

if ~isempty(t_diss_scan)
    fprintf('\tcorr exp(-%.3g/%.3g):\t%.3g [%%]\t\n',...
        t_diss_scan,t1,pol_h*1d2/exp(-t_diss_scan/t1));
end

end      % eval_pol.m


%% sub-function echo2spec
function [spec,fid_ave,hz] = sub_echo2spec(fid,header)
% reconstruct spectra
bw = header.rdb_hdr.user0;           % acquisition bandwidth
si = size(fid);                      % size of matrix
if length(si)~=2, error('fid not 2D'); end

% unchop data (2-step phase cycling -> data stored with opposing signs)
if sub_isodd(header.rdb_hdr.data_collect_type), chp = false;
else, chp = true; end
if chp, fid = fid.*((-1).^(0:si(1)-1).'*ones(1,si(2))); end
fid = conj(fid);                     % direction used for most programs
spec = fftshift(fft(fid,[],2),2);    % actual reconstruction
fid_ave = mean(fid,1);               % average FID
hz = (-si(2)/2:si(2)/2-1)/si(2)*bw;  % frequency axis [Hz]

end      % sub_echo2spec.m


%% sub-function isodd
function bool = sub_isodd(number)
% ISODD  Returns true for odd and false for even input numbers
%
% 7/2004 Rolf Schulte
if ((number-floor(number)) ~= 0), warning('Input number not integer'); end
if ((number/2-floor(number/2)) ~= 0)
    bool = true;
else
    bool = false;
end

end      % sub_isodd