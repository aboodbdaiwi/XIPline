function [spec,hz,s1d,ppm] = recon_csi(d,h,zf,lb,do_ecc,ws,fname,epsi,verb,...
    do_coco,abs_fid,freq,lcm_basis,d_ref,cor_intlv,do_ori,shiftdir)
%RECON_CSI  Reconstruct CSI + EPSI data
%[spec,hz,s1d,ppm]=recon_csi(d,h,zf,lb,do_ecc,ws,fname,epsi,verb,...
%    do_coco,abs_fid,freq,lcm_basis,d_ref,cor_intlv,do_ori,shiftdir)
%                                                               (default)
%        d  Raw p-file data [(nx*ny*nz),ns,1,nt,1,nc] (ns=spectral pts;nc=coils)
%           Interleaved reference scan: nt=1 (1=mrs,2=h2o)
%        h  Header structure
%       zf  Zero fill to:   [#s,#x,#y,#z]
%       lb  Line broadening [s,x,y,z] (arbitrary units)
%   do_ecc  Phase deconvolution eddy current correction         (false)
%           requires H2O reference scan: 
%           if ((nt==2)&&(h.rdb_hdr.user23==1)) -> default=true
%           if ~isempty(d_ref) -> default=true; 
%       ws  Water suppression: 0=off; 1=HSVD; 2=subtraction     (0)
%    fname  Print <fname>.png and save reco as <fname>.mat      ([]) 
%           also export max-image and plots as dicom
%     epsi  EPSI reconstruction (file or structure or from header)
%     verb  Verbose mode + plotting 
%           (0=off;1=verb;2=verb+plot;3=ext verb+plot)          (2)
%  do_coco  Do coil combination                                 (true)
%  abs_fid  Spectral Fourier transform of abs(fid)              (false)
%     freq  Reconstruct spectral domain with matrix inversion   ([])
%lcm_basis  Filename of LCModel basis -> save files             ([])
%    d_ref  Reference scan (w/o water suppression)
%cor_intlv  Phase correction for epsi interleaves               (true)
%   do_ori  Image orientation: x-y-z to AP-RL-SI                (true)
% shiftdir  fftshift directions (zyx;111=shift 3D;000=no shift) (111)
%
%     spec  Spectrum        [#s,#x,#y,#z,#t]
%       hz  Frequency axis  [Hz]
%      s1d  abs(spec) summed spatially+coil -> 1D spectral
%      ppm  Frequency axis  [ppm]
%
% 8/2020 Rolf Schulte
if nargin<1, help(mfilename); return; end

timerVal = tic;               % record execution time


%% misc input
export_dcm = true;           % export max-image + plots as dicom
if ~exist('zf','var'),     zf = []; end
if ~exist('lb','var'),     lb = []; end
if ~exist('do_ecc','var'), do_ecc = []; end
if ~exist('ws','var'),     ws = []; end
if isempty(ws),            ws = 0;  end
if length(ws)~=1,          warning('length(ws)~=1'); end
if ~exist('fname','var'),  fname = []; end
if islogical(fname)
    if fname, error('provide fname not true');
    else, fname = [];
    end
end
if ~isempty(fname)
    if ~islogical(fname)
        if ~isempty(regexpi(fname,'\.7$')), fname = fname(1:end-2); end
        if ~isempty(regexpi(fname,'\.h5$')), fname = fname(1:end-3); end
        if ~isempty(regexpi(fname,'\.mat$')), fname = fname(1:end-4); end
    end
end
if ~exist('epsi','var'),   epsi = []; end
if ~exist('verb','var'),   verb = []; end
if isempty(verb),          verb = 2; end
if length(verb)~=1,        warning('length(verb)~=1'); end
if ~exist('do_coco','var'),do_coco = []; end
if isempty(do_coco),       do_coco = true; end
if length(do_coco)~=1,     warning('length(co_coco)~=1'); end
if ~exist('abs_fid','var'),abs_fid = []; end
if isempty(abs_fid),       abs_fid = false; end
if length(abs_fid)~=1,     warning('length(abs_fid)~=1'); end
if ~exist('freq','var'),   freq = []; end
if ~isempty(freq), if ~isnumeric(freq), warning('~isnumeric(freq)'); end; end
if ~exist('lcm_basis','var'), lcm_basis = []; end
if ~isempty(lcm_basis) 
    if ~islogical(lcm_basis) || ~ischar(lcm_basis)
        error('lcm_basis not logical or char');
    end
end
if ~exist('d_ref','var'),  d_ref = []; end
if ~exist('cor_intlv','var'), cor_intlv = []; end
if isempty(cor_intlv),        cor_intlv = true; end
if length(cor_intlv)~=1,   warning('length(cor_intlv)~=1'); end
if ~exist('do_ori','var'), do_ori = []; end
if isempty(do_ori),        do_ori = true; end
if length(do_ori)~=1,      warning('length(do_ori)~=1'); end
if ~exist('shiftdir','var'), shiftdir = []; end
if isempty(shiftdir),      shiftdir = 111; end
if length(shiftdir)~=1,    warning('length(shiftdir)~=1'); end


%% reading in data, if pfile name is given
if ~isnumeric(d)
    if exist(d,'file')
        [d,h] = read_p(d);
    else
        warning('strange input d/file not existing');
    end
end


%% set default for do_ecc
if isempty(do_ecc)
    do_ecc = false;
    if h.rdb_hdr.user23==1
        if h.image.numecho==2
            do_ecc = true;
        else
            warning('h.rdb_hdr.user23==1; nt(=%g)~=2',h.image.numecho);
        end
    end
    if ~isempty(d_ref), do_ecc = true; end 
end


%% if separate reference -> add as second time step
if ~isempty(d_ref)
    if any(size(d)~=size(d_ref)), error('size(d)~=size(d_ref)'); end
    if size(d,4)~=1,   warning('size(d,4)(=%g) ~= 1',size(d,4)); end
    tmp = zeros(size(d,1),size(d,2),size(d,3),2,size(d,5),size(d,6));
    tmp(:,:,:,1,:,:) = d;
    tmp(:,:,:,2,:,:) = d_ref;
    d = tmp;
    clear tmp d_ref
    nt = 2;
end


%% miscelaneous parameters
ws_bounds = [-100 50];                      % bounds for HSVD water suppression
si = size(d); si = [si ones(1,6-length(si))];
nx = h.rdb_hdr.xcsi;
ny = h.rdb_hdr.ycsi;
nz = h.rdb_hdr.zcsi;
nslices = h.rdb_hdr.nslices;
if isempty(epsi)
    ns = h.rdb_hdr.user1;
else
    ns = h.image.user1;                     % number of sampling points
end
if ~exist('nt','var'), nt = h.image.numecho; end  % number of time steps
nc = (h.rdb_hdr.dab(2)-h.rdb_hdr.dab(1))+1; % number of Rx coil elements
nref = h.rdb_hdr.user19;                    % #reference frames (w/o WS)

if si(3)~=1, warning('si(3)(=%g) ~= 1',si(3)); end
if nt~=si(4)
    warning('nt(=%g) ~= size(d,4)(=%g); setting nt=%g',nt,si(4),si(4));
    nt = si(4);
end
if si(5)~=nslices, warning('si(5)(=%g) ~= nslices(=%g)',si(5),nslices); end
if nc~=si(6)
    warning('nc(=%g) ~= size(d,6)(=%g); setting nc=%g',nc,si(6),si(6));
    nc = si(6);
end
if ((nz>1) && (nslices>1))
    warning('nz(=%g)>1 && nslices(=%g)>1; setting nslices=1',nz,nslices);
    nslices = 1;
end
if nz>1, dim = 3;       % 3D spatially
else,    dim = 2;       % 2D spatially incl. multi-slice
end

switch dim
    case 2, nn = [ns nx ny nslices nt nc];
    case 3, nn = [ns nx ny nz nt nc];
    otherwise, error('dim(=%g) not 2 or 3',dim);
end
bw = h.rdb_hdr.spectral_width;              % sampling bandwidth [Hz]
t  = (0:ns-1).'/bw;                         % time [s]

if ((do_ecc) && (nref>0))
    if nt~=1, warning('nt(=%g)~=1',nt); end
    if verb>0, fprintf('Using reference scan for ECC\n'); end
    tmp = zeros(si(1)-nref,si(2),si(3),2,si(5),si(6));
    tmp(:,:,:,1,:,:) = d(1:(si(1)-nref),:,:,1,:,:);
    % tmp(:,:,:,2,:,:) = repmat(d((si(1)-nref+1),:,:,1,:,:),[(si(1)-nref),1,1,1,1,1]);
    tmp(1,:,:,2,:,:) = d((si(1)-nref+1),:,:,1,:,:);
    d = tmp;
    nt = 2;
    nn(5) = 2;
end



%% load or generate epsi structure
if isempty(epsi)
    % if ~isempty(regexpi(h.image.psd_iname,'epsi'))
    if h.rdb_hdr.user22>0, epsi = header2epsi(h,verb); end
    % end
else
    if ~isstruct(epsi)
        if exist(epsi,'file')
            epsi = load(epsi);
        else
            disp(epsi);
            warning('strage input for epsi/file not found');
        end
    end
end


%% misc info + warnings
if verb>0
    fprintf('psd_iname=%s\n',h.image.psd_iname);
    fprintf('ns=%g, nx=%g, ny=%g, nz=%g, nslices=%g, nt=%g, nc=%g\n',...
        ns,nx,ny,nz,nslices,nt,nc);
    fprintf('nucleus=%g\n',h.rdb_hdr.user2);
end
if isempty(epsi)
    if nx*ny*nz~=si(1), warning('nx*ny*nz~=si(1)'); end
    if ns~=si(2),       warning('ns(=%g)~=si(2)(=%g)',ns,si(2)); end
end
if (do_ecc && (nt~=2))
    warning('do_ecc && (nt(=%g)~=2',nt);
    fprintf('reference required; setting ''do_ecc = false''\n');
    do_ecc = false;
end
if ((ws==2) && (nt~=2))
    warning('ws==2 && (nt(=%g)~=2; reference required',nt);
    fprintf('reference required; setting ''ws = 0''\n');
    ws = 0;
end


%% unchop data
if ischop(h)
    if isempty(epsi)
        d(2:2:end,:) = -d(2:2:end,:);
    else
        ilv = epsi.intlv;
        for ii=1:ilv, d((ilv+ii):(2*ilv):end,:) = -d((ilv+ii):(2*ilv):end,:); end
    end
end


%% reorder data
if verb>0, fprintf('Reordering data\n'); end
% size(d)  = [#x*#y*#z, #s,  1, #t,  1, #c]
% size(dd) = [#s,       #x, #y, #z, #t, #c]
if isempty(epsi)
    % if verb>0, fprintf('Reshaping data\n'); end
    dd = zeros(nn);
    for l2=1:nn(2)
        for l3=1:nn(3)
            for l4=1:nz
                for l5=1:nn(5)
                    ll1 = l2+(l3-1)*nn(2)+(l4-1)*nn(3)*nn(2);
                    for l6=1:nn(6)
                        for lsli=1:nslices
                            dd(:,l2,l3,l4+lsli-1,l5,l6) = d(ll1,:,1,l5,lsli,l6).';
                        end
                    end
                end
            end
        end
    end
else
    %% reshaping for EPSI (optional)
    if verb>0, fprintf('EPSI preprocessing\n'); end
    if ~isstruct(epsi), disp(epsi); error('epsi not a structure'); end
    if isfield(epsi,'fidall')
        nn(1) = epsi.npts(2);
        nn(epsi.dir+2) = epsi.npts(1);
        bw = epsi.bws(2);
    end
    dd = epsi_reshape(d,nn,epsi.ind,epsi.intlv,epsi.dir+1,epsi.phafu,...
        cor_intlv,epsi.mode,verb-2);
end


%% (Gaussian) line broadening - both spatial and spectral
if ~isempty(lb)
    if verb>0, fprintf('Gaussian apodisation\n'); end
    lb = [lb(:) ; zeros(4,1)];
    lb_opts = {'','h','h','h'};
    lb_bw = [bw nn(2:4)];
    for l=1:(dim+1)
        if lb(l)~=0
            lbf = lb_fun(nn(l),lb_bw(l),[0 lb(l)],lb_opts{l});
            dd = bsxfun(@times,dd,reshape(lbf,[ones(1,l-1),nn(l),1]));
        end
    end
end


%% spatial zero filling (Fourier interpolation)
if ~isempty(zf)
    if verb>0, fprintf('Zero-filling (spatial)\n'); end
    zf = [zf(:).', nn((length(zf)+1):6)];
    if ~isempty(freq) && (zf(1)~=ns)
        warning('zf(1)(=%g) ~= ns(=%g)',zf(1),ns);
    end
    if any(zf<nn)
        warning('zf<nn'); 
        fprintf('zf = '); disp(zf)
        fprintf('nn = '); disp(nn)
        zf = max([nn;zf],[],1);
        fprintf('Adjusting zf values\n');
        fprintf('zf = '); disp(zf)
    end
    if nslices>1
        if zf(4)>nslices
            warning('zf(4)(=%g) > nslices(=%g); setting zf(4)=%g',...
                zf(4),nslices,nslices);
            zf(4) = nslices;
        end
    end
    if ((nslices==1) && (nz==1) && (zf(4)>1))
        warning('zf(4)(=%g)>1 and 2D single slice (nz=%g && nslices=%g)',...
            zf(4),nz,nslices);
        fprintf('  setting zf(4) to 1\n');
        zf(4) = 1;
    end
    
    % spatial zf; manually, because spatial data symmetric
    tmp = zeros(nn(1),zf(2),zf(3),zf(4),nt,nc);
    for l=2:4
        if isodd(zf(l)-nn(l))
            warning('zf(%g)(=%g) - nn(%g)(=%g) = %g; rounding',l,zf(l),l,nn(l),zf(l)-nn(l)); 
        end
        x{l} = ceil(((zf(l)-nn(l))/2+1):((zf(l)+nn(l))/2)); 
    end
    tmp(:,x{2},x{3},x{4},:,:) = dd;
    dd = tmp;
else
    zf = nn;
end

hz  = -(-zf(1)/2:zf(1)/2-1)/zf(1)*bw;        % frequency axis
ppm = hz/h.rdb_hdr.ps_mps_freq/1d-7;


%%  spatial reconstruction
for l=2:(dim+1)
    % fprintf('l=%g; 4-l=%g, 10^(4-l)=%g; zf(l)=%g\n',l,4-l,10^(4-l),zf(l));
    if zf(l)>1
        if verb>0, fprintf('Spatial reconstruction: dim=%g',l); end
        dd = ifft(fftshift(dd,l),[],l);
        if isodd(floor(shiftdir*10^(2-l)))
            if verb>0, fprintf('; ifftshift'); end
            dd = ifftshift(dd,l);
        end
        if verb>0, fprintf('\n'); end
    end
end


%% coil combination
if do_coco && nc>1
    if verb>0, fprintf('Coil combination\n'); end
    dd = csi_coil_combine(dd);
end


%% eddy current correction - phase deconvolution
if do_ecc
    if verb>0
        fprintf('Eddy current correction by phase deconvolution\n');
    end
    dd = csi_ecc(dd);
end


%% water subtraction via HSVD
if ws==1
    if verb>0
        fprintf('HSVD water subtraction: bounds=(%g %g) [Hz]\n',...
            ws_bounds(1),ws_bounds(2));
    end
    [~,~,~,dd(:,:,:,:,1,:)] = hsvd(dd(:,:,:,:,1,:),25,ws_bounds,bw,1);
end


%% spectral reconstruction
% size(spec) = [#s,#x,#y,#z,#t]
if isempty(freq)
    % fft recon
    if zf(1)>1
        if verb>0, fprintf('Spectral reconstruction\n'); end
        if abs_fid
            if verb>0, fprintf('\ttaking abs(dd) only\n'); end
            dd = abs(dd);
        end
        spec = ifftshift(fft(dd,zf(1),1),1);
    else
        if verb>0, fprintf('Skipping spectral reconstruction\n'); end
    end
else
    % chemical shift modelling
    if verb>0
        fprintf('Spectral reconstruction by chemical shift inversion\n');
    end
    A = exp(-1i*2*pi*t*freq(:).');
    pinv_A = pinv(A);
    zf(1) = length(freq);
    spec = complex(zeros(zf));
    for l2=1:zf(2)
        for l3=1:zf(3)
            for l4=1:zf(4)
                for l5=1:zf(5)
                    for l6=1:zf(6)
                        spec(:,l2,l3,l4,l5,l6) = pinv_A*dd(:,l2,l3,l4,l5,l6);
                    end
                end
            end
        end
    end
    hz = freq;
end


%% Water subtraction
if ws==2
    if verb>0, fprintf('Water subtraction\n'); end
    indws = abs(hz)<30;
    spec = csi_h2osubtraction(spec,indws);
end


%% image orientation: xyz -> AP-RL-SI
if do_ori
    spec = csi_orientation(spec,h);
end


%% appr. 1D spectrum
if ((nargout>2) || (verb>1)) || ~isempty(fname)
    s1d = zeros(zf(1),2,nt);
    for lt=1:nt
        if abs_fid
            sss = sum(real(spec(:,:,:,:,lt,:)),6);
        else
            sss = abs(sqrt(sum(spec(:,:,:,:,lt,:).*conj(spec(:,:,:,:,lt,:)),6)));
        end
        s1d(:,1,lt) = squeeze(sum(sum(sum(sum(sss,2),3),4),5));
        s1d(:,2,lt) = squeeze(max(max(max(max(sss,[],2),[],3),[],4),[],5));
    end
    s1d = s1d./repmat(max(s1d,[],1),[zf(1) 1]);
end
figstr = sprintf('P%05d Exam%d Series%d ',...
    h.image.rawrunnum,h.exam.ex_no,h.series.se_no);
par.fax = hz; par.fig_str = figstr; par.lb = lb; % par structure for plot_mrsi


%% saving result as mat file
if ~isempty(fname)
    if verb>0, fprintf('Saving results (fname=''%s'')\n',fname); end
    tmp = whos('spec');
    if verb>0, fprintf('Memory size of spec = %g [MB] (SI=Base10)\n',tmp.bytes*1d-6); end
    if tmp.bytes>100d6   % reduce size for data bigger than 100MB
        if verb>0 
            fprintf('\t>100MB: converting spec (64bit) into single (32bit))\n'); 
        end
        spec = single(spec); % convert: double (64bit) to single (32bit) precision
    end
    
    save([fname '.mat'],'-mat','-v7.3',...
        'spec','h','zf','lb','do_ecc','ws','fname','epsi','verb',...
        'do_coco','abs_fid','freq','lcm_basis','cor_intlv','do_ori',...
        'shiftdir','hz','s1d','nn','t','dim',...
        'si','nx','ny','nz','nslices','ns','nt','nc','par');
end


%% plotting
if ((verb>1) && (ns>1))
    if (length(size(spec))>5) || ((verb==3)&&(isempty(fname)))
        warning('(length(size(spec))>5) || ((verb==3)&&(isempty(fname))); skipping plotting');
    else
        fprintf('Plotting\n');
        if abs_fid, pc = 0; else, pc = []; end
        if verb==3
            fid = figure('Visible','off');
        else
            figure(20); clf
        end
        imagesc_row(shiftdim(max(abs(spec),[],1),1),[],'row','',true);
        set(gcf,'name',[figstr '- image of max']);
        colormap gray
        if ~isempty(fname)
            print([fname '_maximage'],'-dpng','-r600','-painters');
        end
        if verb==3, close(fid); end
        
        if ~isempty(s1d)
            if verb==3
                fid = figure('Visible','off');
            else
                figure(21); clf
            end
            plot(hz,s1d(:,:,1).');
            set(gca,'XDir','reverse');
            legend('rms','max');
            set(gcf,'name',[figstr '- pseudo 1D spectra']);
            grid on
            if ~isempty(fname)
                print([fname '_s1d'],'-dpng','-r600','-painters');
                if export_dcm>0
                    write_scdicom([fname '_s1d.dcm'],gcf,h);
                end
            end
            if verb==3, close(fid); end
        end
        
        if ~isempty(fname)
            fname_plt_csi = [fname '_plot_csi'];
        else
            fname_plt_csi = [];
        end
        if export_dcm>0
            plot_csi(spec,pc,[],[figstr '- CSI'],fname_plt_csi,11,h);
        else
            plot_csi(spec,pc,[],[figstr '- CSI'],fname_plt_csi,11);
        end
        
        if verb~=3, plot_mrsi(spec,par); end
    end
end


%% exporting images of max as dicom into scanner database
if ~isempty(fname) && export_dcm
    if verb>0, fprintf('Storing dicom images\n'); end
    inp.SeriesDescription = sprintf('CSI max:S%g;X%g;Y%g;Z%g,T%g,C%g', ...
        ns,nx,ny,nz+nslices-1,nt,nc);
    if dim==3, inp.MRAcquisitionType = '3D'; end
    write_dicom(fname,shiftdim(max(abs(spec),[],1),1),h,inp);
end


%% write out to LCModel
if ~isempty(lcm_basis)
    if verb>0, fprintf('Writing LCModel files\n'); end
    if islogical(lcm_basis)
        if ~lcm_basis, return; end
        te = h.image.te*1d-3;
        switch te
            case 35,  lcm_basis = '/home/schulte/data/basis/lcmodel_simulated/press35/gamma_basis_TE35_B0_30_TE1_18_BW5000_N4096.basis';
            case 144, lcm_basis = '/home/schulte/data/basis/lcmodel_simulated/press144/gamma_basis_TE144_B0_30_TE1_20_BW5000_N4096.basis';
            otherwise
                lcm_basis = '/home/schulte/data/basis/lcmodel_simulated/press144/gamma_basis_TE144_B0_30_TE1_20_BW5000_N4096.basis';
                warning('TE(=%g)~=144[ms] (taking TE=144ms as basis)',h.image.te*1d-3); 
        end
    end
    if ~ischar(lcm_basis), error('~ischar(lcm_basis)'); end
    if verb>0, fprintf('Writing out raw data and control file for LCModel\n'); end
    if ~isempty(lb), warning('Data line-broadened'); end
    if (zf~=ns),     warning('Data zero-filled');    end
    
    for l2=1:zf(2)
        for l3=1:zf(3)
            for l4=1:zf(4)
                for l5=1:zf(5)
                    lcm_fname = sprintf('lcm_matlab/P%.5d/P%.5d_l2_%g_l3_%g_l4_%g_l5_%g',...
                        h.image.rawrunnum,h.image.rawrunnum,l2,l3,l4,l5);
                    if verb>0, fprintf('%s\n',lcm_fname); end
                    write_lcmodel(dd(:,l2,l3,l4,l5),h,lcm_basis,[],lcm_fname);
                end
            end
        end
    end
end

if verb>0
    fprintf('recon_mrsi.m finished: ');
    fprintf('total reconstruction time = %.2f s\n',toc(timerVal));
end

end      % main function recon_csi.m
