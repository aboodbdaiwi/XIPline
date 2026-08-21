function freq = plot_spec_noise(d,h,sv,absfreq,thresh,nmax,frqtxt,max_spec)
%PLOT_SPEC_NOISE Plot spectra in different ways to visualise noise
%                Use >=64 repetitions, BW>=125kHz; little/no RF
%   freq = plot_spec_noise(d,h,sv,absfreq,thresh,nmax)
%      d   raw P-file data
%          optionally filename to P-file or image archive
%      h   header structure
%     sv   Save results; file name or logical                       (false)
%absfreq   Plot x-axis in absolute frequencies                       (true)
% thresh   Factor for determining maxima: thresh*mean(spec(:))          (2)
%   nmax   Restrict #peaks to nmax                                    (100)
% frqtxt   Plot text with max frequencies                           (false)
%max_spec  Maximum of spectrum (for plotting)                          ([])
%   freq   Absolute frequencies [Hz]
%
%  6/2024 Rolf Schulte
if nargin<1, help(mfilename); return; end

if ~exist('sv','var'),    sv = []; end
if isempty(sv),           sv = false; end
if ~exist('absfreq','var'),absfreq = []; end
if isempty(absfreq),      absfreq = true; end
if ~exist('thresh','var'),thresh = []; end
if isempty(thresh),       thresh = 2; end
if ~exist('nmax','var'),  nmax = []; end
if isempty(nmax),         nmax = 100; end
if ~exist('frqtxt','var'), frqtxt = []; end
if isempty(frqtxt),       frqtxt = false; end
if ~exist('max_spec','var'), max_spec = []; end

docog = false;     % refine peak freq via centre-of-gravity


%% reading in data, if pfile name is given
if ~isnumeric(d) || isempty(d)
    if ~exist(d,'file') && ~isempty(d)
        warning('strange input for d (=''%s''); file not existing?',d); 
    end
    [d,h] = read_p(d);
end


%% filenames
if ~exist('fname','var'), fname = sprintf('P%.5d.7',h.image.rawrunnum); end
tmp = regexpi(fname,'/');         % if directory, take name only
if ~isempty(tmp), fname = fname((tmp(end)+1):end); end
tmp = regexpi(fname,'\');         % same for win dir separators
if ~isempty(tmp), fname = fname((tmp(end)+1):end); end

if ischar(sv)
    % if filename given in sv field
    plt_fname = sv;
    if ~isempty(regexpi(plt_fname,'\.7$', 'once')),  plt_fname = plt_fname(1:end-2); end
    if ~isempty(regexpi(plt_fname,'\.h5$', 'once')), plt_fname = plt_fname(1:end-3); end
    sv = true;
else
    % remove suffix
    if sv
        tmp = regexpi(fname,'\.');
        if ~isempty(tmp)
            plt_fname = fname(1:(tmp(end)-1));
        else
            plt_fname = fname;
        end
    end
end


%% spectral reconstruction
[n1,n2,n3,n4,n5,n6] = size(d);
if length(size(d))>6
    warning('length(size(d))(=%d)>6',length(size(d)));
end
if n3>1, warning('size(d,3)(=%d)>1',n3); end
if n4>1, warning('size(d,4)(=%d)>1',n4); end
if n5>1, warning('size(d,5)(=%d)>1',n5); end
nex = h.image.nex;
if nex~=1, warning('h.image.nex(=%g)~=1',nex); end
if any([n3,n4,n5,nex]>1), pause(1); end

r1 = h.rdb_hdr.ps_mps_r1;
r2 = h.rdb_hdr.ps_mps_r2;
scale = 2^((r1-12)/2+(r2-30))*1d6; % scale signal to R1=12,R2=30

% bug in version<21.7.2023: fft -> freq axis inverted
spec = ifftshift(ifft(d,[],2),2)/scale;
spec = sum(sum(abs(conj(spec).*spec),1),6)/(n1*n6)*n2;
if n5>1, spec = sum(spec,5)/n5; end


%% absolute frequency
bw = h.rdb_hdr.spectral_width;
f0 = h.rdb_hdr.ps_mps_freq/10;
hz = (-n2/2:n2/2-1)/n2*bw;
if absfreq
    hz_plt = hz + f0;
else
    hz_plt = hz;
end


%% determine spur frequencies
mspec = mean(spec(:));
ind = []; freq = []; amp = [];
if thresh>1d-6
    [mm,ii] = sort(spec,'descend');
    lmax = 0;
    for ll=1:n2
        if mm(ll)<thresh*mspec, break; end
        if lmax>=nmax, break; end
        ismax = true;
        if ii(ll)<n2
            if mm(ll)<spec(ii(ll)+1), ismax = false; end
        end
        if ii(ll)>1
            if mm(ll)<spec(ii(ll)-1), ismax = false; end
        end
        if ismax
            lmax = lmax+1;
            ind(1,lmax) = ii(ll);
            if ~docog, freq(1,lmax) = hz(ii(ll)); end
            amp(1,lmax) = mm(ll);
        end
    end

    % refine via centre-of-gravity
    if docog
        for ll=1:length(ind)
            if ind(1,ll)>1 && ind(1,ll)<n2
                ii = ind(1,ll) + (-1:1);
                freq(1,ll) = sum(hz(ii).*spec(1,ii))/sum(spec(1,ii));
            else
                waning('ind(1,ll)(=%d) ==1 or %d',ind(1,ll),n2);
                freq(1,ll) = hz(ind(1,ll));
            end
        end
    end
end


%% plotting
figid = figure; clf
figstr = sprintf('nuc=%d f0=%.6g[MHz] BW=%g[kHz] se%d P%d',...
    h.image.specnuc,h.rdb_hdr.ps_mps_freq*1d-7,h.rdb_hdr.user0*1d-3,...
    h.series.se_no,h.image.rawrunnum);
if nex>1, figstr = sprintf('%s nex=%d',figstr,nex); end
if h.rdb_hdr.dacq_data_type~=2
    figstr = sprintf('%s EDR=off',figstr);
end
set(gcf,'name',figstr);

% spec = log(spec);
plot(hz_plt,spec,'',...
    [min(hz_plt) max(hz_plt)],thresh*mspec*[1 1],'r:',...
    freq+absfreq*f0,amp,'xm');
if isempty(max_spec), max_spec = max(spec); end
axis([min(hz_plt) max(hz_plt) 0 max_spec]);
set(gca,'xdir','reverse'); grid on;
xlabel('freq [Hz]'); 
% title('power spectrum');
title(figstr);
if frqtxt && ~isempty(freq)
    for ll=1:length(freq)
        text(freq(ll)+absfreq*f0,amp(ll),num2str(round(freq(ll)+f0)),...
            'FontSize',8);
    end
end


%% info
istr = sprintf('Power spectrum:\nspec = ifftshift(ifft(d,[],2),2)/scale;\n');
istr = sprintf('%sspec = sum(sum(abs(conj(spec).*spec),1),6)/(n1*n6)*n2;\n',istr);
istr = sprintf('%sscale = 2^((r1-12)/2+(r2-30))*1d6 = %g\n',istr,scale);
istr = sprintf('%s------------------\n',istr);
istr = sprintf('%sBW(cv0)   = %g [Hz]\npts(cv1)  = %g\n#exc(cv4) = %g\n',...
    istr,h.rdb_hdr.user0,h.rdb_hdr.user1,h.rdb_hdr.user4);
istr = sprintf('%sdf  = %g [Hz]\n',istr,h.rdb_hdr.user0/h.rdb_hdr.user1);

istr = sprintf('%sR1  = %d\nR2  = %d\nf0  = %d [Hz]\nnex = %d\n',...
    istr,r1,r2,round(f0),nex);
istr = sprintf('%s------------------\n',istr);
istr = sprintf('%smean(d(:))    = %g +1i*%g\n',istr,real(mean(d(:))),imag(mean(d(:))));
istr = sprintf('%sstd(d(:))     = %g\n',istr,std(d(:)));
istr = sprintf('%smean(spec(:)) = %g\n',istr,mean(spec(:)));
istr = sprintf('%sstd(spec(:))  = %g\n',istr,std(spec(:)));
if thresh>1d-6
    % [mm,ii] = sort(spec,'descend');
    mspec = mean(spec(:));
    istr = sprintf('%sthreshold     = %g\n',istr,thresh*mspec);
    istr = sprintf('%snmax          = %g\n',istr,nmax);
    istr = sprintf('%s------------------\n',istr);
    istr = sprintf('%sind\tabs freq\trel freq\tvalue\n',istr);
    for ll=1:length(ind)
        istr = sprintf('%s%d\t%d\t%d\t\t%g\n',istr,ll,round(freq(ll)+f0),...
            round(freq(ll)),amp(ll));
    end
end
fprintf(istr);


%% saving to file
if sv>0
    if isdeployed
        print([plt_fname '_psn.png'],'-dpng','-r300','-painters');
    else
        print([plt_fname '_psn.png'],'-dpng','-r300');
    end
    ff = fopen([plt_fname '_info.txt'],'wt');
    fprintf(ff,'%s',istr);
    fclose(ff);
    saveas(figid,[plt_fname '_psn.fig'],'fig');
    if isdeployed
        write_scdicom([plt_fname '.dcm'],gcf,h);
    end
end


%% output
if exist('freq','var'), freq = freq+f0; end


end      % plot_spec_noise.m
