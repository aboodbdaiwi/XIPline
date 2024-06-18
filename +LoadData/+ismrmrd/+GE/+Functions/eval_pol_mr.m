function eval_pol_mr(fname_h,fname_t,t_diss_scan,ind_t1,flip_rf1_h,flip_rf1_t)
%EVAL_POL_MR  Evaluate polarisation level for GE FIDCSI sequence
%  FIDCSI (eg hardpulse; sw=5000; npt=2048; dda=0; nex=1; TG const.)
%  Hyperpol FID: reduced R1/R2, flip_rf1=1.4; pw_rf1=8; optr=5s; cv4=32
%  Thermal  FID: R1=13;R2=30; flip_rf1=90; pw_rf1=512; optr=2s; cv4=512
%                doping to reduce T1: ~3:100 DOTAREM (or alike)
%                (e.g., 2ml PA + ~60mg DOTAREM)
%  Signal level obtained by averaging first points in FID
%   -> Requires one predominant resonance:
%   -> Pyruvate solution should be about pH neutral
%
%   eval_pol_mr(fname_h,fname_t,t_diss_scan,ind_t1,flip_rf1_h,flip_rf1_t)
%      fname_h  P-file name hyperpolarised signal (or *.mat)
%      fname_t   (")   thermally polarised signal (or *.mat)
%  t_diss_scan  Time from dissolution to acquisition  (opt)  [s]
%       ind_t1  Index of last acq for fitting T1      (opt)
%   flip_rf1_h  Flip angle of hyperpolarised FID      (opt)  [deg]
%   flip_rf1_t  Flip angle of thermally polarised FID (opt)  [deg]
%                by default taken from header
% 6/2008  Rolf Schulte
% 7/2013  Rolf Schulte
%See also READ_MR_RAWDATA.
if (nargin<1), help(mfilename); return; end;

%% reading data
if strcmpi(fname_h(end-3:end),'.mat')
    tmp = load(fname_h); d_h = tmp.d; h_h = tmp.h; clear tmp;
else
    [d_h,h_h] = read_MR_rawdata(fname_h); 
end
if strcmpi(fname_t(end-3:end),'.mat')
    tmp = load(fname_t); d_t = tmp.d; h_t = tmp.h; clear tmp;
else
    [d_t,h_t] = read_MR_rawdata(fname_t);
end
if length(size(d_h))~=2, error('d_h matrix must be 2D'); end
if length(size(d_t))~=2, error('d_t matrix must be 2D'); end


%% get parameters from header or input
if ~exist('t_diss_scan','var'), t_diss_scan = []; end
if ~exist('ind_t1','var'),      ind_t1 = []; end
if ind_t1>size(d_h,1), error('ind_t1>size(d_h,1)'); end
if abs(ind_t1-floor(ind_t1))>1d-15, error('ind_t1 not integer'); end
if ~exist('flip_rf1_h','var'),  flip_rf1_h = []; end
if ~exist('flip_rf1_t','var'),  flip_rf1_t = []; end

indave = (3:15);                  % index for averaging first point in FID
if isempty(ind_t1), ind_t1 = size(d_h,1); end  % last index for fitting T1
% flip angles [deg]
if isempty(flip_rf1_h), flip_rf1_h = h_h.image.mr_flip; end
if isempty(flip_rf1_t), flip_rf1_t = h_t.image.mr_flip; end
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

fprintf('\n*************************\nParameters\n');
fprintf('   hyperpol: R1=%g, R2=%g, TG=%g, flip_rf1=%g [deg], TR=%g [s]\n',...
    r1_h,r2_h,tg_h,flip_rf1_h,tr_h);
fprintf('    thermal: R1=%g, R2=%g, TG=%g, flip_rf1=%g [deg], TR=%g [s]\n',...
    r1_t,r2_t,tg_t,flip_rf1_t,tr_t);

if (tg_h ~= tg_t), warning('eval_pol_3t:tg','different TGs'); end


%% reconstructions
[s_h,dave_h,hz_h] = echo2spec_sub(d_h,h_h);
sig_h = max(mean(abs(d_h(:,indave)),2));

[s_t,dave_t,hz_t] = echo2spec_sub(d_t,h_t);
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
tt = (0:ind_t1-1)*tr_h;
t1decay = mean(abs(d_h((1:ind_t1),indave)),2).';

p = polyfit(tt,log(t1decay),1);
t1fit = polyval(p,tt);
t1 = -1/p(1);

%% plotting 
subplot(3,2,1); 
plot(t_h*1d3,abs(d_h(1,:)),'b',t_h(indave)*1d3,abs(d_h(1,indave)),'g-x'); 
axis tight; ylabel('hyper pol'); xlabel('FID [ms]');
subplot(3,2,2); 
plot(hz_h,abs(s_h).'); axis tight; xlabel('spectrum');
set(gca,'XDir','reverse');   % spectroscopy convention
subplot(3,2,3); 
plot(t_t*1d3,abs(dave_t),'b',t_t(indave)*1d3,abs(dave_t(indave)),'g-x'); 
axis tight; ylabel('thermal pol'); xlabel('FID [ms]');
subplot(3,2,4); 
plot(hz_t,abs(mean(s_t,1))); axis tight; xlabel('spectrum');
set(gca,'XDir','reverse');   % spectroscopy convention
subplot(3,2,5); 
plot(tt,t1decay,'b-x',tt,exp(t1fit),'r-x'); 
axis tight; xlabel('T1 decay [s]');
subplot(3,2,6); 
plot(tt,log(t1decay),'b-x',tt,t1fit,'r-x'); 
axis tight; xlabel('T1 fit (log) [s]');

%% print output
fprintf('\n***********\n');
fprintf('T1 = %.3g [s]\n',t1);
fprintf('Thermal polarisation = %.3g [ppm]\n',pol_t*1d6);
fprintf('Hyper-polarisation\n');
fprintf('\tuncorrected:      \t%.3g [%%]\t\n',pol_h*1d2);

if ~isempty(t_diss_scan),
    fprintf('\tcorr exp(-%.3g/%.3g):\t%.3g [%%]\t\n',...
        t_diss_scan,t1,pol_h*1d2/exp(-t_diss_scan/t1));
end


%% subfunctions
function [spec,fid_ave,hz] = echo2spec_sub(fid,header)
% reconstruct spectra
bw = header.rdb_hdr.user0;           % acquisition bandwidth
si = size(fid);                      % size of matrix
if length(si)~=2, error('fid not 2D'); end

% unchop data (2-step phase cycling -> data stored with opposing signs)
if isodd_sub(header.rdb_hdr.data_collect_type), chp = false;
else chp = true; end
if chp, fid = fid.*((-1).^(0:si(1)-1).'*ones(1,si(2))); end
fid = conj(fid);                     % direction used for most programs
spec = fftshift(fft(fid,[],2),2);    % actual reconstruction
fid_ave = mean(fid,1);               % average FID
hz = (-si(2)/2:si(2)/2-1)/si(2)*bw;  % frequency axis [Hz]

function bool = isodd_sub(number)
% ISODD  Returns true for odd and false for even input numbers
%
% 7/2004 Rolf Schulte
if number-floor(number) ~= 0
  warning('isodd_sub:int','Input number not integer');
end
if number/2-floor(number/2) ~= 0
  bool = true;
else
  bool = false;
end
