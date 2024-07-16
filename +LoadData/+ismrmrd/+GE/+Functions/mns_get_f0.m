function mns_get_f0(d,h,fname)
%MNS_GET_F0  Extract f0 from 1H MRS data and
% mns_get_f0(d,h,fname)
%     d  Raw data (or filename for loading)
%     h  p-file header structure
% fname  Save picture and text output
%
% 12/2020 Rolf Schulte
if (nargin<1), help(mfilename); return; end


%% input + default parameters
if ~exist('fname','var'), fname = []; end
plt = true;                  % plot to matlab figure
do_dcm = true;               % export prescan window to dicom database
iwidth = 10;                 % width of centre-of-mass for freq determination


%% reading in data, if pfile name is given
if ~isnumeric(d)
    if exist(d,'file')
        [d,h] = read_p(d);
    else
        warning('strange input d/file not existing');
    end
end


%% strip possible suffix from fname
if ~isempty(fname)
   if ~isempty(regexpi(fname,'\.png$')), fname = fname(1:end-4); end 
   if ~isempty(regexpi(fname,'\.txt$')), fname = fname(1:end-4); end 
   if ~isempty(regexpi(fname,'\.h5$')),  fname = fname(1:end-3); end 
   if ~isempty(regexpi(fname,'\.7$')),   fname = fname(1:end-2); end 
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
[nexc,npts,n3,n4,n5,ncoils] = size(d);
if n3~=1, warning('n3(=%d)~=1; rms combination',n3); end
if n4~=1, warning('n4(=%d)~=1; rms combination',n4); end
if n5~=1, warning('n5(=%d)~=1; rms combination',n5); end
if ncoils~=1, fprintf('ncoils(=%d)>1; rms combination\n',ncoils); end
zf = max([8192 2*npts]);               % interpolate to zf


%% misc variables + checks
if h.image.specnuc~=1
    error('h.image.specnuc(=%g)~=1',h.image.specnuc);
end
bw = h.rdb_hdr.user0;                  % (full) BW [Hz]
f0cur = h.rdb_hdr.ps_mps_freq/10;


%% reconstruct data
if ischop(h) && (nexc>1), d(2:2:nexc,:) = -d(2:2:nexc,:); end
if ~isempty(regexpi(h.image.psd_iname,'presscsi'))
    % extract unsuppressed water signal
    nref = h.rdb_hdr.user19;           % #reference frames (w/o WS)
    nspec = h.image.user4/h.image.nex; % #spectral frames  (w/ WS)
    if nexc~=(nref+nspec)
        warning('nexc(=%g)~=(nref(=%g)+nspec(=%g))',nexc,nref,nspec);
    end
    d = d(1:nref,:,:,:,:,:,:);
    nexc = nref;
end
spec = ifftshift(fft(d,zf,2),2);
spec = mean(spec,1);                   % mean along repetitions
spec = sqrt(mean(mean(mean(mean(abs(spec.*conj(spec)),6),5),4),3)); % rms along rest
spec = spec/max(spec);                 % normalise
hz = -(-zf/2:zf/2-1)/zf*bw;            % freq axis [Hz]


%% centre of mass f0 determination
[~,imax] = max(spec);
ind = (imax-iwidth):(imax+iwidth);
if any(ind<1)
    warning('ind<1: truncating');
    ind = ind(ind>0);
end
if any(ind>zf)
    warning('ind>zf(=%d): truncating',zf);
    ind = ind(ind<=zf);
end
df0 = sum(hz(ind).*spec(1,ind))/sum(spec(1,ind));
f0 = f0cur+df0;
fprintf('1H: f0 = %g+%g = %d [Hz]\n',f0cur,df0,round(f0));


%% determine linewidths (full width at half maximum)
[lw,out] = fwhm(spec,bw,[],false,false);
fprintf('fwhm = %.4g [Hz]\n',lw);


%%
[~,str] = mns_freq_converter(f0);
f13Clac = mns_freq_converter(f0,'13C','lactate',false);
f129Xe = mns_freq_converter(f0,'129Xe','lung',false);


%% plot spectra
if plt
    plot(hz,spec,'b',...
         df0*[1 1],[-0.02 1.2],'g:',...
         out.freqs(1,:),spec(1,out.ind_max(1,:)),'m.',...
         [out.fp(1,:) out.fm(1,:)],...
         repmat(spec(1,out.ind_max(1,:))/2,[2 1]),'m.-');
    xlabel('freq [Hz]'); ylabel('Spectrum [a.u.]'); axis tight;
    grid on
    set(gca,'xdir','reverse');
    text(-bw/40,1.15,sprintf('1H H2O f0=%d[Hz]',round(f0)),'FontSize',12);
    text(-bw/40,1.08,sprintf('13C Lac f0=%d[Hz]',round(f13Clac)),'FontSize',12);
    text(-bw/40,1.01,sprintf('129Xe lung f0=%d[Hz]',round(f129Xe)),'FontSize',12);
    
    
    if ~isempty(fname)
        if isdeployed
            print([fname '.png'],'-dpng','-r600','-painters');
        else
            print([fname '.png'],'-dpng','-r600');
        end
        if do_dcm
            write_scdicom([fname '.dcm'],gcf,h);
        end
    end
end



%%
if ~isempty(fname)
    outid = fopen([fname '.txt'],'w+');     % direct to file
    fprintf(outid,str);
    fclose(outid);
end



end      % mns_get_f0.m
