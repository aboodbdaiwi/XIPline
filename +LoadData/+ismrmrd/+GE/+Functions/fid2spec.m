function [spec,par,fid,hz,ppm,t] = fid2spec(fid,header,zf,lb,bw,chp,verb,coil_comb)
%FID2SPEC  Reconstruct single voxel echo spectrum
%[spec,par,fid,hz,ppm,t] = fid2spec(data,header,zf,lb,bw,chp,verb,coil_comb)
%
%     data  Free induction decay (sampled along 2nd dim)
%   header  GE header file structure
% *** following parameters are optional ***
%       zf  zero fill to this value
%       lb  linebroadening (exp gauss)                [Hz]
%       bw  resample data to this BW                  [Hz]
%      chp  chopping              (default by header)
%     verb  Verbose mode               (default=true)
%coil_comb  Combine coils (via svds)   (default=true) 
%
%     spec  spectrum (not averaged; same size as data)
%      par  parameter structure as used for plot_spec
%      fid  data (t-domain; chopping + conj)
%       hz  x-axis [Hz]
%      ppm  x-axis [ppm]
%        t  time-axis (for FID)
%
% 9/2015 Rolf Schulte
%See also PLOT_SPEC, PAR2HEADER.
if nargin<1, help(mfilename); return; end

if ~exist('zf','var'), zf = []; end
if ~exist('lb','var'), lb = []; end
if ~exist('bw','var'), bw = []; end
if ~exist('chp','var'), chp = []; end
if ~exist('verb','var'), verb = []; end
if isempty(verb), verb = true; end
if ~exist('coil_comb','var'), coil_comb = []; end
if isempty(coil_comb), coil_comb = true; end

if isempty(chp), chp = ischop(header,verb); end
if chp && (size(fid,1)>1), fid(2:2:end,:) = -fid(2:2:end,:); end

if coil_comb && size(fid,6)>1
    fid = mrs_coil_combine(fid,[],[],[],verb);
end

par = header2par(header);         % convert header to par structure
fid = conj(fid);                  % changing dir of axis
% "conj" is used instead of "ifft" to keep consistency with spectroscopy 
% (incl GAMMA), where the most common convention is spec=fft(fid).
% Note: in MRI on GE systems (in P-files), the convention is spec=ifft(fid). 


%% parameter checks
[n1,n2,n3] = size(fid);
si = size(fid);
if n2~=par.samples, warning('samples inconsistent'); end
if isfield(header,'rdb_hdr')
    if header.rdb_hdr.nframes~=n1
        warning('header.rdb_hdr.nframes~=n1');
        fprintf('Setting par.rows to size(d,1) (= %g)\n',n1);
        par.rows = n1;
    end
end

if isempty(zf), zf = n2; end      % no zero-filling
if zf~=si(2), si(2)=zf; end
spec = NaN(si);
if ~isempty(bw)    % prepare some stuff for down-sampling
    if bw>par.bw, error('not implemented: can only sample down'); end
    nn2orig = zf*bw/par.bw;       % new #samples (float)
    nn2 = 2^ceil(log2(nn2orig));  % round up to power of 2
    si_bw = si; si_bw(2) = nn2;
    t1 = (0:zf-1)*1/par.bw;       % orig time
    t2 = linspace(0,1-1/nn2orig,floor(nn2orig))*nn2orig/bw; % new time
    fid_bw = zeros(si_bw);
    spec_bw = NaN(si_bw);
end

for l3=1:n3
    % line broadening
    if ~isempty(lb)
        if verb, fprintf('Apodising\n'); end
        apo = lb_fun(n2,par.sample_frequency,lb);
        fid(:,:,l3) = fid(:,:,l3).*repmat(apo,[n1 1]);
    end

    spec(:,:,l3) = fftshift(fft(fid(:,:,l3),zf,2),2);  % actual recon

    % resampling: FID scaling the same; spec scaled by bw/par.bw
    if ~isempty(bw)
        spec(:,[1:ceil((zf-nn2orig)/2),floor((zf+nn2orig)/2):zf],l3) = 0;
        tmp = ifft(ifftshift(spec(:,:,l3),2),[],2);
        
        for l1=1:n1
            fid_bw(l1,1:floor(nn2orig),l3) = interp1(t1,tmp(l1,:),t2);
            % Here: same fid scaling (area under curve constant, when
            % downsampling)
            % Add "*par.bw/bw;", if same spectral scaling desired
        end
        spec_bw(:,:,l3) = fftshift(fft(fid_bw(:,:,l3),[],2),2); 
    end
end


%% write new info into par structure
if ~isempty(bw)
    fid = fid_bw;
    spec = spec_bw;
    par.samples = nn2;
    par.sample_frequency = bw;
    par.bw = bw;
else
    par.samples = zf;
end
if any(isnan(spec)), warning('spec contains NaN'); end


%% x-axes
n = par.samples;
hz = (-n/2:n/2-1)/n*par.bw;
ppm = hz/par.f0*1d6;
t = (0:header.rdb_hdr.user1-1)/header.rdb_hdr.user0;


end      % fid2spec.m
