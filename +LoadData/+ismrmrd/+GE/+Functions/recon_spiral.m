function [bb,bbabs] = recon_spiral(d,h,wfn,mtx_reco,delay,lbt,spks,...
    fname,plt,ave,freqs,metstr,f0bymax,T2s,export_dcm)
%RECON_SPIRAL Reconstruct single-shot spiral trajectory
%[bb,bbabs] = recon_spiral(d,h,wfn,mtx,delay,lbt,spks,fname,plt,ave,...
%                 freqs,metstr,f0bymax,T2s,export_dcm)
%                                                               (default)
%         d  Raw (p-file) data  (or pfile fname)
%         h  Header from p-file (or empty)
%       wfn  Location of waveform .mat file
%       mtx  Reconstruction matrix size (Fourier interpolation) ([])
%     delay  Gradient-acquisition delay [us]                    (0) 
%            1d or 2d=[WHOLE,ZOOM]
%       lbt  Gaussian linebroadening time [Hz]                  (0)
%      spks  Remove spike noise                                 (false)
%     fname  Print <fname>.png and save reco as <fname>.mat     ([]) 
%            also export dicom if template.dcm is found
%       plt  Plotting: 0=off(in background), 1=on               (true)
%       ave  Average data                                       (false)
%     freqs  Freqs for IDEAL [Hz]          ([392 267 177 0 -324];LHAPB)
%            highest peak=0Hz
%            if empty -> look for values in exrf<cv14>.cvs/ideal<cv21>.cvs
%    metstr  Metabolite names (cells with strings)              ([])
%   f0bymax  Determine f0 of maximum peak and centre spectra to it 
%            (overwritten by ideal<cv21>.cvs)                   (true)
%       T2s  Include T2* decay in CSM (IDEAL inversion) [ms]    ([])
%export_dcm  Save images in dicom format                        (1)
%
%        bb  Reconstructed data (mtx,mtx,#exc,#metabolites,#slice,#coils)
%     bbabs  RMS-coil combined images  (mtx,mtx,#exc,#metabolites,#slice)
%
% 12/2020  Rolf Schulte
if (nargin<1), help(mfilename); return; end


%% input variables
if ~exist('mtx_reco','var'), mtx_reco = []; end
if isempty(mtx_reco),        mtx_reco = 0; end
if length(mtx_reco)~=1
    warning('length(mtx)(=%g)~=1; taking first value only',length(mtx_reco));
    mtx_reco = mtx_reco(1);
end
if ~exist('delay','var'),    delay = []; end
if isempty(delay),           delay = 0; end
if ~exist('lbt','var'),      lbt = []; end
if isempty(lbt),             lbt = 0; end
if ~exist('spks','var'),     spks = false; end
if ischar(spks)
    spks = str2num(spks);
    if isempty(spks),        spks = false; end
end
if ~exist('fname','var'),    fname = []; end
if ~isempty(fname)
    if ~isempty(regexpi(fname,'\.7')), fname = fname(1:end-2); end
end
if ~exist('plt','var'),      plt = []; end
if isempty(plt),             plt = true; end
if ~exist('ave','var'),      ave = false; end
if ~exist('freqs','var'),    freqs = []; end
if ~exist('metstr','var'),   metstr = []; end
if ~exist('f0bymax','var'),  f0bymax = []; end
if isempty(f0bymax),         f0bymax = true; end
if ~exist('T2s','var'),      T2s = []; end
if ~exist('export_dcm','var'),export_dcm = []; end
if isempty(export_dcm),      export_dcm = true; end
timerVal = tic;              % record execution time
if isdeployed
    renderer = '-painters';
else
    renderer = '-opengl';
end


%% reading in data, if pfile name is given
if ~isnumeric(d)
    if exist(d,'file')
        [d,h] = read_p(d);
    else
        warning('strange input d/file not existing');
    end
end


%% reading waveform file
if isempty(wfn), error('wfn empty'); end
if isstruct(wfn)
    wf = wfn;
else
    if ~isempty(regexpi(wfn,'\.wav$')), wfn = wfn(1:end-4); end
    if isempty(regexpi(wfn,'\.mat$')),  wfn = [wfn '.mat']; end
    if ~exist(wfn,'file'), error('file not found: wfn=%s',wfn); end
    wf = load(wfn);          % load waveform
end


%% different WHOLE and ZOOM delay times (delay->array)
if length(delay)>1
    switch h.rdb_hdr.grad_mode
        case 1, delay = delay(1);
        case 2, delay = delay(2);
        otherwise
            error('length(delay)(=%g)>1 & h.rdb_hdr.grad_mode(=%g) ~= 1 or 2',...
                length(delay),h.rdb_hdr.grad_mode);
    end
end


%% special reconstructions
reco = 'normal';
if h.rdb_hdr.user22>0
    reco = 'spsp';
end
if h.rdb_hdr.user21>0
    reco = 'ideal';
end


%% load recon frequencies + metabolite names
if isempty(freqs)
    if h.rdb_hdr.user40>0
        % frequencies from header.rdb_hdr.user4X
        for l=1:h.rdb_hdr.user40
            freqs(l) = eval(['h.rdb_hdr.user4' num2str(l)]);
        end
        % metabolite names from file with cv pars (3rd column)
        switch reco
            case 'spsp',  cvfname = sprintf('exrf%g.cvs',h.image.user14);
            case 'ideal', cvfname = sprintf('ideal%g.cvs',h.image.user21);
        end
        if ~exist(cvfname,'file')
            cvfname = [getenv('homedir') '/waveforms/' cvfname];
        end
        if exist(cvfname,'file')
            fprintf('File with pars found\n\t%s\n',cvfname);
            fp = fopen(cvfname,'r');
            tmp = textscan(fp,'%s %f %s');
            fclose(fp);
            for l=1:length(freqs) 
                for ll=1:length(tmp{1})
                    if strfind(tmp{1}{ll},['_freq' num2str(l)])
                        ms = tmp{3}{ll};
                        if isempty(ms), ms = num2str(tmp{2}(ll)); end
                        metstr{l} = ms;
                    end
                    if strfind(tmp{1}{ll},'f0bymax')
                        f0bymax = logical(tmp{2}(ll));
                    end
                end
            end
        else
            fprintf('NO file with pars found\n\t%s\n',cvfname);
         end
    else
        if isempty(metstr)
            [freqs,metstr] = freqs13c;
        else
            freqs = freqs13c;
        end
        if f0bymax, freqs = freqs - freqs(4); end
    end
end

if isempty(metstr)
     for l=1:length(freqs), metstr{l} = num2str(freqs(l)); end
end


%% check for fields required for reconstruction
if ~isfield(wf,'k'),   error('wf.k not existing '); end
if ~isfield(wf,'dcf'), error('wf.dcf not existing'); end
if isempty(wf.dcf),    error('wf.dcf is empty'); end
if any(isnan(wf.dcf)), error('wf.dcf contains NaN'); end
if any(isinf(wf.dcf)), error('wf.dcf contains inf'); end
if ~isfield(wf,'ind'), error('wf.ind not existing'); end
if ~isfield(wf,'t'),   error('wf.t not existing'); end
if ~isfield(wf,'npix'),error('wf.npix not existing'); end
% if ~isfield(wf,'fov'), error('wf.fov not existing'); end

mtx_acq = wf.npix;           % nominal matrix resolution of acquisition
if mtx_reco <= 0, mtx_reco = mtx_acq; end
k = double(wf.k*mtx_acq/mtx_reco); % scale k-space -> Fourier interpolation
fov = h.rdb_hdr.fov*1d-3;    % actual field-of-view [m]


%% misc variables to checkings
bw = h.rdb_hdr.user0;
% n_freqs = h.rdb_hdr.user22;
% n_exc = h.rdb_hdr.user4;
% n_ts = n_exc/n_freqs;
% n_slice = size(d,5);
[nexc,npts,n3,n4,nslices,ncoils] = size(d);
xloc = h.rdb_hdr.user26;
yloc = h.rdb_hdr.user27;
zloc = h.rdb_hdr.user28;
rot90fac = h.data_acq_tab.rotate;
trnsps = h.data_acq_tab.transpose;

if ~((trnsps(1)==0) || (trnsps(1)==3))
    warning('h.data_acq_tab.transpose (=%g) not 0 or 3: unknown image orientation',trnsps(1));
end
if n3~=1, warning('size(d,3)(=%g)~=1',n3); end
if nexc<size(wf.ind,1)
    warning('nexc(=%g)<size(wf.ind,1)(=%g)',nexc,size(wf.ind,1));
end
tmp = nexc/size(wf.ind,1);
if abs(tmp-floor(tmp))>1d-10 
    warning('nexc(=%g) not multiple of size(ind,1)(=%g); truncating data',...
        nexc,size(wf.ind,1));
    d = d(1:floor(tmp)*size(wf.ind,1),:,:,:,:,:);
    nexc = size(d,1);
end


%% pre-processing of data
if spks, d = remove_spikes(d,[],false); end
if ((abs(xloc)>0.01) || (abs(yloc)>0.01))
    shft = mtx_reco*[xloc yloc]*1d-3/fov;
else
    shft = [];
end
% tmp=sort(wf.dcf);
% cart_down_fac = 1.1*mean(tmp(floor(length(tmp)*0.75):floor(length(tmp)*0.95)));
% if cart_down_fac>1, cart_down_fac = []; end
% if cart_down_fac<0.1
%     warning('cart_down_fac<0.1; setting to 0.1');
%     cart_down_fac = 0.1; 
% end
cart_down_fac = [];
ind_pc = [];       % if index is stored in wf -> perform phase correction
if isfield(wf,'ind_pc'), ind_pc = wf.ind_pc; end
dd = raw2grid(d,ischop(h),k,shft,cart_down_fac,[],wf.ind,delay*1d-6,bw,...
    false,ind_pc);

% dim1=#coils; dim2=#indexed kspace data; dim3=#excitations; dim4=#slices
if lbt>0
    % dd = ak_apodise(dd,[],[],wf.t(:).',lbt,false); 
    t = double(wf.t.');
    t = t(:).';
    dd = ak_apodise(dd,[],[],t,lbt,false); 
end
if ave, dd = sum(dd,3); nexc = 1; reco = 'normal'; end


%% optional IDEAL spiral CSI spectral reconstruction
if strcmpi(reco,'ideal')
    fidfirst = h.rdb_hdr.user39;       % first frame FID (gradients nulled)
    if fidfirst==0, f0bymax = false; end
    nte = h.rdb_hdr.user21;            % # of shifted echoes
    dte = h.rdb_hdr.user20*1d-6;       % echo time shift
    % reconstruct spectrum
    dspec = double(d(1:(fidfirst+nte):end,:,:,:,:,:));
    nzf = max(npts,8192);
    spec = fftshift(fft(dspec,nzf,2),2);
    fax = -bw*(-nzf/2:nzf/2-1)/nzf;
    
    % look for max signal -> Pyr -> all other shifts relative to this
    if f0bymax
        [tmp1,i1]=max(abs(spec),[],1);
        [tmp2,i2]=max(abs(tmp1),[],2);
        [tmp4,i4]=max(abs(tmp2),[],4);
        [tmp5,i5]=max(abs(tmp4),[],5);
        [tmp6,i6]=max(abs(tmp5),[],6);
        imax = [i1(i2(i5(i6))) , i2(i5(i6)) , 1 , 1 , i5(i6) , i6];
        spvec = abs(spec(imax(1),:,imax(3),imax(4),imax(5),imax(6)));
        % figure;plot(fax,spvec);set(gca,'XDir','reverse');
        iwidth = ceil(max([2 , 20/(bw/nzf) ])); % width at least +-20Hz or +-2pts
        cog_ind = ((imax(2)-iwidth):(imax(2)+iwidth));
        cog_freq = sum(fax(cog_ind).*spvec(1,cog_ind))/sum(spvec(1,cog_ind));
        fprintf('Pyruvate frequency = %g [Hz]\n',cog_freq);
    else
        cog_freq = 0;
    end
    freq = freqs + cog_freq;
    freq = freq(:);
    
    % misc
    nmeta = length(freq);
    nk = length(k);
    t = (0:nte-1).'*dte;
    A = exp(-1i*2*pi*t*freq.');          % CS kernel
    if ~isempty(T2s),A = A.*repmat(exp(-t/T2s),[1,nmeta]); end
    pinv_A = pinv(A);
    % d1 = zeros(1,nk,nmeta);
    % tt = repmat(t,[1 nk])+repmat((0:nk-1)/bw,[nte 1]);
    tt = (0:nk-1)/bw;
    
    % plotting
    if fidfirst
        ss = sqrt(sum(spec.*conj(spec),6));    % coil combination
        ss = ss/(max(ss(:)))+0.1;
        if plt || ~isempty(fname)
            for l5=1:size(spec,5)
                if plt
                    figure;
                else
                    fid = figure('Visible','off');
                end
                set(gcf,'DefaultAxesFontSize',14);
                figstr = sprintf('P%05d Exam%d Series%d Spectra Slice%g',...
                    h.image.rawrunnum,h.exam.ex_no,h.series.se_no,l5);
                tmp = log(ss(:,:,1,1,l5).')+0.3*repmat((1:size(ss,1)),[size(ss,2) 1]);
                plot(fax-cog_freq,tmp,'');
                hold on
                plot(freqs.'*[1 1],[-10 10],'r:','LineWidth',2);
                plot(-cog_freq*[1 1],[-10 10],'b:','LineWidth',2);
                hold off
                axis([round([min(freqs)-300 max(freqs)+300]/10)*10 min(tmp(:)) max(tmp(:))]);
                xlabel('Frequency [Hz]');
                ylabel('Spectra (log) / time course');
                
                f0 = h.rdb_hdr.ps_mps_freq/10;
                tstr=['freqs=' num2str(freqs) ...
                    '; df=' num2str(-floor(cog_freq+0.5)) ...
                    '; f0=' num2str(f0) '[Hz]'];
                %    '\newline' num2str(floor(freqs+f0-cog_freq*0.5))];
                tstr=regexprep(tstr,'\s+',' ');
                title(tstr,'FontSize',12);
                
                set(gca,'XDir','reverse');
                drawnow;
                set(gcf,'name',figstr);
                if ~isempty(fname)
                    print([fname '_spec_ls' num2str(l5)],'-dpng','-r600',renderer);
                    if export_dcm>0
                        inpscd.InstanceNumber = l5;
                        write_scdicom([fname '_spec_ls' num2str(l5) '.dcm'],...
                            gcf,h,inpscd);
                    end
                end
                if ~plt, close(fid); end
            end
        end
    end
end
clear d dspec spec


%%
switch reco
    case 'normal'
        % ntimesteps = nexc;
        ntimesteps = size(dd,3);
        nmeta = 1;
    case 'spsp'
        nmeta = h.rdb_hdr.user22;
        ntimesteps = floor(nexc/nmeta);
        % ntimesteps = floor(size(dd,3)/nmeta);
        if size(wf.t,1)>1         % spsp combined with multi-arm spiral
            arms = size(wf.t,1);
            nmeta = nmeta/arms;
            if abs(nmeta-floor(nmeta))>1d-10
                warning('nmeta(=%g) not multiple of arms(=%g)',nmeta,arms);
            end
            tmp = metstr; clear metstr;
            for l=1:nmeta, metstr{l} = tmp{((l-1)*arms+1)}; end
        end
    case 'ideal'
        % ntimesteps = floor(nexc/(fidfirst+nte));
        ntimesteps = floor(nexc/(fidfirst+nte))*n4;
end
bb = zeros(mtx_reco,mtx_reco,ntimesteps,nmeta,nslices,ncoils);
% b = zeros(mtx_reco,mtx_reco,ncoils);
% nk = length(k);tt = (0:nk-1)/bw;freq=920;


%% actual reconstruction
for lt=1:ntimesteps
    for ls=1:nslices
        for lc=1:ncoils
            switch reco
                case 'normal'
                    b1 = gridding(dd(lc,:,lt,ls),k,wf.dcf,mtx_reco);
                    % b1 = grid_matlab(dd(lc,:,lt,ls),k,mtx_reco);
                    if trnsps(1)==3
                        bb(:,:,lt,1,ls,lc) = rot90(fliplr(b1).',rot90fac(1));
                    else
                        bb(:,:,lt,1,ls,lc) = rot90(fliplr(b1),rot90fac(1));
                    end
                case 'spsp'
                    for lm=1:nmeta
                        ll = lm + (lt-1)*nmeta;
                        b1 = gridding(dd(lc,:,ll,ls),k,wf.dcf,mtx_reco);
                        % b1 = gridding(dd(lc,:,ll,ls).*exp(1i*2*pi*tt*freq),...
                        %    k,wf.dcf,mtx_reco);
                        if trnsps(1)==3
                            bb(:,:,lt,lm,ls,lc) = rot90(fliplr(b1).',rot90fac(1));
                        else
                            bb(:,:,lt,lm,ls,lc) = rot90(fliplr(b1),rot90fac(1));
                        end
                    end
                    if isempty(metstr), metstr = {'Lac','Pyr','BC','Ala'}; end
                case 'ideal'
                    ind_t = (lt-1)*(nte+fidfirst) + fidfirst + (1:nte);
                    d1 = reshape(((pinv_A*squeeze(dd(lc,:,ind_t,ls)).').* ...
                        exp(1i*2*pi*freq*tt)).',[1,nk,nmeta]);
                    for lm=1:nmeta
                        b1 = gridding(d1(1,:,lm),k,wf.dcf,mtx_reco);
                        if trnsps(1)==3
                            bb(:,:,lt,lm,ls,lc) = rot90(fliplr(b1).',rot90fac(1));
                        else
                            bb(:,:,lt,lm,ls,lc) = rot90(fliplr(b1),rot90fac(1));
                        end
                    end
            end
        end
    end
end
clear dd


%% coil combine
if ncoils>1
    bbabs = sqrt(mean(bb.*conj(bb),6));
else
    bbabs = abs(bb);
end


%% saving result as mat file
if ~isempty(fname)
    if isempty(metstr), metstr = ''; end
    bb = single(bb);   % convert from double (64bit) to single (32bit) precision
    bbabs = abs(single(bbabs));
    if ~exist('ss','var'), ss = []; end
    if ~exist('fax','var'), fax = []; end
     
    save([fname '.mat'],'-mat','-v7.3',...
        'bb','h','wfn','mtx_reco','delay','lbt','spks','fname',...
        'cart_down_fac','metstr','bbabs','ss','fax','bbabs');
end


%% plotting
if plt || ~isempty(fname)
    if nmeta==1
        % plotting for pure imaging
        if plt
            figure;
        else
            fid = figure('Visible','off');
        end
        
        rshp = false;
        if (ntimesteps==1) && (nslices>8), rshp = true; end
        set(gcf,'DefaultAxesFontSize',14);
        imagesc_row(reshape(bbabs,[mtx_reco,mtx_reco,ntimesteps,nslices]),...
            [],'ind',rshp,true);
        colormap(gray);
        figstr = sprintf('P%05d Exam%d Series%d',...
            h.image.rawrunnum,h.exam.ex_no,h.series.se_no);
        set(gcf,'name',figstr);
        if ~isempty(fname)
            print(fname,'-dpng','-r600',renderer);
        end
        if ~plt, close(fid); end
    else
        % plotting for metabolic images
        for l5=1:nslices
            if plt
                figure;
            else
                fid = figure('Visible','off');
            end
            set(gcf,'DefaultAxesFontSize',14);
            scale = 1.5*ones(ntimesteps,1)* ...
                squeeze(max(max(max(bbabs(:,:,:,:,l5),[],1),[],2),[],3)).';
            imagesc_row(bbabs(:,:,:,:,l5),[],1./scale);
            figstr = sprintf('P%05d Exam%d Series%d Slice%d',...
                h.image.rawrunnum,h.exam.ex_no,h.series.se_no,l5);
            ylstr = '';
            if strcmpi(reco,'ideal')
                for lm=nmeta:-1:1, ylstr = sprintf('%s  %s',ylstr,metstr{lm}); end
            end
            if strcmpi(reco,'spsp')
                for lm=nmeta:-1:1, ylstr = sprintf('%s  %s',ylstr,metstr{lm}); end
            end
            xlabel('time steps');
            ylabel(ylstr);
            pos = get(gca,'Position');
            set(gca,'Position',pos.*[0.6 1 1.15 1]);
            set(gcf,'name',figstr);
            if ~isempty(fname)
                print([fname '_ls' num2str(l5)],'-dpng','-r600',renderer);
            end
            if ~plt, close(fid); end
        end
    end
end


%% exporting images as dicom into scanner database
if ~isempty(fname) && (export_dcm>0)
    tmp = sprintf('MNSRP;TS%g,MET%g,SL%g',ntimesteps,nmeta,nslices);
    if nmeta>1
        tmp = [tmp ';'];
        for l=1:length(metstr), tmp = [tmp metstr{l}(1)]; end
    end
    tmp = [h.series.se_desc ';' tmp ';' reco];
    inp.SeriesDescription = tmp;
    if export_dcm>99
        inp.SeriesNumber = h.exam.ex_numseries + export_dcm;
    end
    % matlab-based dicomwrite.m without template
    % dicom header generated from p-file header
    write_dicom(fname,bbabs,h,inp);    
end


%%
fprintf('recon_spiral.m finished: ');
fprintf('total reconstruction time = %.2f [s]\n',toc(timerVal));


end      % main function recon_spiral.m
