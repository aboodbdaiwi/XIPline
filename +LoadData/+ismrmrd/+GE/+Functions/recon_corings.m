function [spec,hz,s1d,ppm] = recon_corings(d,h,wfn,zf,lb,do_ecc,ws,fname,verb,...
    do_coco,abs_fid,d_ref,do_ori,do_single)
%RECON_CORINGS  Reconstruct MRSI data acquired with concentric rings
%[spec,hz,s1d,ppm]=recon_concentric_rings(d,h,wfn,zf,lb,do_ecc,ws,fname,verb,...
%    do_coco,abs_fid,freq,lcm_basis,d_ref,do_ori)
%                                                               (default)
%        d  Raw p-file data 
%           [(n_ring*n_intlv*n_zpe),ns,1,n_time,n_slices,n_coils] 
%           Interleaved reference scan: nt=1 (1=mrs,2=h2o)
%        h  Header structure
%      wfn  Location of waveform .mat file
%       zf  Zero fill to:   [#spec,#x,#y,#z]
%       lb  Line broadening [spec,x,y,z]         [Hz,mtx]       (0)
%   do_ecc  Phase deconvolution eddy current correction         (false)
%           requires H2O reference scan: 
%           if ((nt==2)&&(h.rdb_hdr.user23==1)) -> default=true
%           if ~isempty(d_ref) -> default=true; 
%       ws  Water suppression: 0=off; 1=HSVD; 2=subtraction     (0)
%    fname  Print <fname>.png and save reco as <fname>.mat (opt,[]); 
%           also export image of max as dicom if template.dcm is found
%     verb  Verbose mode + plotting 
%           (0=off;1=verb;2=verb+plot;3=ext verb+plot)          (2)
%  do_coco  Do coil combination                                 (true)
%  abs_fid  Spectral Fourier transform of abs(fid)              (false)
%    d_ref  Reference scan (w/o water suppression)
%   do_ori  Image orientation: x-y-z to AP-RL-SI                (true)
%do_single  Convert and reconstruct data in single precision    ([])
%           (false=double; true=single; []->single if size d>1GB)
%
%     spec  Spectrum        [#s,#x,#y,#z,#t]
%       hz  Frequency axis  [Hz]
%      s1d  abs(spec) summed spatially+coil -> 1D spectral
%      ppm  Frequency axis  [ppm]
%
% 6/2020 Rolf Schulte
if nargin<1, help(mfilename); return; end

timerVal = tic;               % record execution time


%% misc input
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
if ~exist('verb','var'),   verb = []; end
if isempty(verb),          verb = 2; end
if length(verb)~=1,        warning('length(verb)~=1'); end
if ~exist('do_coco','var'),do_coco = []; end
if isempty(do_coco),       do_coco = true; end
if length(do_coco)~=1,     warning('length(co_coco)~=1'); end
if ~exist('abs_fid','var'),abs_fid = []; end
if isempty(abs_fid),       abs_fid = false; end
if length(abs_fid)~=1,     warning('length(abs_fid)~=1'); end
if ~exist('d_ref','var'),  d_ref = []; end
if ~exist('cor_intlv','var'), cor_intlv = []; end
if isempty(cor_intlv),     cor_intlv = true; end
if length(cor_intlv)~=1,   warning('length(cor_intlv)~=1'); end
if ~exist('do_ori','var'), do_ori = []; end
if isempty(do_ori),        do_ori = true; end
if length(do_ori)~=1,      warning('length(do_ori)~=1'); end
if ~exist('do_single','var'), do_single = []; end


%% reading in data, if pfile name is given
if ~isnumeric(d)
    if exist(d,'file')
        [d,h] = read_p(d);
    else
        if ischar(d)
            error('file ''%s'' not found',d);
        else
            warning('strange input: d');
        end
    end
end


%% use single precision for large data
if isempty(do_single)
    if prod(size(d))*8*2>1d9
        do_single = true;
    else
        do_single = false;
    end
end
if do_single
    fprintf('Converting input data to single\n');
    d = single(d); 
end


%% reading waveform file
if isempty(wfn), error('wfn empty'); end
if isstruct(wfn)
    wf = wfn;
else
    if ~isempty(regexpi(wfn,'\.k$'))
        if ~exist(wfn,'file'), error('file not found: wfn=%s',wfn); end
        wf = create_corings_wf(h,wfn);
    else
        if ~isempty(regexpi(wfn,'\.wav$')), wfn = wfn(1:end-4); end
        if isempty(regexpi(wfn,'\.mat$')),  wfn = [wfn '.mat']; end
        if ~exist(wfn,'file'), error('file not found: wfn=%s',wfn); end
        wf = load(wfn);          % load waveform
    end
end


%% set default for do_ecc
if isempty(do_ecc)
    do_ecc = false;
    if h.rdb_hdr.user23==1
        if h.image.numecho==2
            do_ecc = true;
        else
            warning('h.rdb_hdr.user23==1; n_time(=%g)~=2',h.image.numecho);
        end
    end
    if ~isempty(d_ref), do_ecc = true; end 
end


%% if separate reference -> add as second time step
if ~isempty(d_ref)
    if any(size(d)~=size(d_ref)), error('size(d)~=size(d_ref)'); end
    if size(d,4)~=1,   warning('size(d,4)(=%g) ~= 1',size(d,4)); end
    tmp = complex(zeros(size(d,1),size(d,2),size(d,3),2,size(d,5),size(d,6)));
    tmp(:,:,:,1,:,:) = d;
    tmp(:,:,:,2,:,:) = d_ref;
    d = tmp;
    clear tmp d_ref
    n_time = 2;
end


%% miscelaneous parameters
ws_bounds = [-100 100];                      % bounds for HSVD water suppression
si = size(d); si = [si ones(1,6-length(si))];
mtx = wf.mtx(:).'; mtx = [mtx ones(1,3-length(mtx))];
nx = mtx(1);
ny = mtx(2);
nz = mtx(3);
n_slices = h.rdb_hdr.nslices;
if ~exist('nt','var'), n_time = h.image.numecho; end      % #time steps
n_coils = (h.rdb_hdr.dab(2)-h.rdb_hdr.dab(1))+1;           % #Rx coil elements
nref = h.rdb_hdr.user19;               % #reference frames (w/o WS)
ns = wf.npts;                          % #sampling points
ind = wf.ind;                          % index to select acquired data
n_spec = wf.n_spec;                    % #spectral points
nk_circ = wf.nk_circ;                  % #sample points per 1 ring
bw_spec = wf.bw_spec;                  % spectral bandwidth [Hz]
n_ring = wf.n_ring;                    % #rings
n_intlv = wf.n_intlv;                  % #spectral interleaves
n_zpe = wf.n_zpe;                      % #z-phase encodings
if abs(rem(nk_circ,n_intlv))>1d-10
    warning('nk_circ(=%g) not multiple of n_intlv(=%g)',nk_circ,n_intlv);
end
if abs(rem(n_spec,n_intlv))>1d-10
    warning('n_spec(=%g) not multiple of n_intlv(=%g)',n_spec,n_intlv);
end

if size(d,1)<(n_ring*n_zpe*n_intlv)
    error('size(d,1)(=%d)<n_ring(=%d)*n_zpe(=%d)*n_intlv(=%d)',...
        size(d,1),n_ring,n_zpe,n_intlv);    
end
if size(d,1)>(n_ring*n_zpe*n_intlv)
    warning('size(d,1)(=%d)>n_ring(=%d)*n_zpe(=%d)*n_intlv(=%d): truncating d',...
        size(d,1),n_ring,n_zpe,n_intlv);
    ii = 1:(n_ring*n_zpe*n_intlv);
    d = d(ii,:,:,:,:,:);
end
if si(2)~=ns, warning('si(2)(=%g)~=ns(=%g)',si(2),ns); end
if si(3)~=1, warning('si(3)(=%g) ~= 1',si(3)); end
if n_time~=si(4)
    warning('nt(=%g) ~= size(d,4)(=%g); setting nt=%g',n_time,si(4),si(4));
    n_time = si(4);
end
if si(5)~=n_slices
    warning('si(5)(=%g) ~= nslices(=%g); setting n_slices=%d',...
        si(5),n_slices,si(5));
    n_slices = si(5);
end
if n_coils~=si(6)
    warning('nc(=%g) ~= size(d,6)(=%g); setting nc=%g',n_coils,si(6),si(6));
    n_coils = si(6);
end
if ((nz>1) && (n_slices>1))
    warning('nz(=%g)>1 && nslices(=%g)>1; truncating slices',nz,n_slices);
    n_slices = 1;
    d = d(:,:,:,:,1,:);
end
nz = max([nz n_slices]);               % use nz for 2D multi-slice and 3D
dim = wf.dim;                          % spatially 2D or 3D
nn = [n_spec nx ny nz n_time n_coils]; % final data size (6D)




if ((do_ecc) && (nref>0))
    if n_time~=1, warning('nt(=%g)~=1',n_time); end
    if verb>0, fprintf('Using reference scan for ECC\n'); end
    tmp = zeros(si(1)-nref,si(2),si(3),2,si(5),si(6));
    tmp(:,:,:,1,:,:) = d(1:(si(1)-nref),:,:,1,:,:);
    tmp(1,:,:,2,:,:) = d((si(1)-nref+1),:,:,1,:,:);
    d = tmp;
    n_time = 2;
    nn(5) = 2;
end
xloc = h.rdb_hdr.user26;
yloc = h.rdb_hdr.user27;
zloc = h.rdb_hdr.user28;
fov = h.rdb_hdr.fov; 
if ((abs(xloc)>0.01) || (abs(yloc)>0.01) || (abs(zloc)>0.01))
    shft = mtx.*[xloc yloc zloc]./repmat(fov,[1 3]);
else
    shft = [];
end


%% misc info + warnings
if verb>0
    fprintf('psd_iname=%s; nucleus=%g\n',h.image.psd_iname,h.rdb_hdr.user2);
    fprintf('ns=%g, nx=%g, ny=%g, nz=%g, nslices=%g, nt=%g, nc=%g\n',...
        ns,nx,ny,nz,n_slices,n_time,n_coils);
    fprintf('n_spec=%g; nk_circ=%g; n_ring=%g; n_intlv=%g; n_zpe=%g\n',... 
        n_spec,nk_circ,n_ring,n_intlv,n_zpe);
    fprintf('bw_spec=%g[Hz]; bw_spat=%g[kHz]\n',... 
        bw_spec,2*h.rdb_hdr.bw);
end
if wf.nexc~=si(1), warning('wf.nexc(=%g)~=si(1)(=%g)',wf.nexc,si(1)); end
if (do_ecc && (n_time~=2))
    warning('do_ecc && (nt(=%g)~=2',n_time);
    fprintf('reference required; setting ''do_ecc = false''\n');
    do_ecc = false;
end
if ((ws==2) && (n_time~=2))
    warning('ws==2 && (nt(=%g)~=2; reference required',n_time);
    fprintf('reference required; setting ''ws = 0''\n');
    ws = 0;
end


%% unchop data
if ischop(h), d(2:2:end,:) = -d(2:2:end,:); end


%% selecting indexed data
if verb>0, fprintf('Selecting indexed data\n'); end
d = d(:,ind(1,:),:,:,:,:);


%% spatial shift
% shft = []; warning('shft=[]');
if ~isempty(shft)
    kk = wf.k;
    if dim>2, kk = repmat(kk,[1 n_zpe 1]); end
    kk = -1i*2*pi*permute(reshape(kk,[size(d,2) size(d,1) 2]),[2 1 3]);
    if verb>0, fprintf('Shifting data by (%g,%g,%g) pixel\n',shft); end
    shft2 = reshape(shft(1,1:2),[1 1 2]);
    phafu = exp(-sum(bsxfun(@times,kk,shft2),3));
    d = bsxfun(@times,d,phafu);
    if ((dim>2) && (abs(shft(1,3))>0.1))
        warning('shft(3) not yet implemented');
    end
end


%% reshaping of data
% size(d)  = [(#rings*#zpe*#intlv), #acq, 1, #time, #slices, #coils]
% size(dd) = [#spec,  #k-space one ring, 1, #z-phase/slices, #time, #coils]
if verb>0, fprintf('Reshaping data\n'); end

dd = complex(zeros(n_spec,n_ring*nk_circ,1,max([n_zpe n_slices]),...
    n_time,n_coils,class(d)));

if verb>2, fprintf('l_spec  i2(1) i2(end) ii1(1) ii2(1) ii2(end)\n'); end
for lc=1:n_coils
    for lt=1:n_time
        for l_slice=1:n_slices
            for lz=1:n_zpe
                for l_spec=1:n_spec
                    if n_zpe>1
                        l4 = lz;
                    else
                        l4 = l_slice;
                    end
                    if n_intlv==1
                        ii1 = (1:n_ring) + (lz-1)*n_ring;
                        ii2 = (1:nk_circ) + (l_spec-1)*nk_circ;
                        dd(l_spec,:,1,l4,lt,lc) = ...
                            reshape(d(ii1,ii2,1,lt,l_slice,lc).',...
                            [1 n_ring*nk_circ]);
                    else
                        for li=1:n_intlv
                            i2 = bsxfun(@plus,(1:(nk_circ/n_intlv)).',...
                                (0:(n_ring-1))*nk_circ+(li-1)*nk_circ/n_intlv);
                            i2 = i2(:).';
                            ss1 = mod(li-l_spec,n_intlv)+1;
                            ii1 = ss1:n_intlv:(n_intlv*n_ring);
                            ii2 = (1:(nk_circ/n_intlv)) + ...
                                (l_spec-1)*nk_circ/n_intlv;
                            if verb>2
                                fprintf('%6d %6d %6d %6d %6d %6d\n',...
                                    l_spec,i2(1),i2(end),ii1(1),ii2(1),ii2(end));
                            end
                            tmp = reshape(d(ii1,ii2,1,lt,l_slice,lc).',...
                                [1 n_ring*nk_circ/n_intlv]);
                            dd(l_spec,i2,1,l4,lt,lc) = tmp;
                        end
                    end
                end
            end
        end
    end
end
clear d


%% correct for time evolution during one circle
if verb>0, fprintf('Correct for time evolution\n'); end
% phafu = repmat(exp(2*pi*1i*wf.t1*bw_spec),[1 n_ring*n_intlv]);
% dd = fft(bsxfun(@times,ifft(dd,[],2),phafu),[],2);
% dd = fft(bsxfun(@times,ifft(dd,[],1),phafu),[],1);
% dd = fft(bsxfun(@times,ifft(dd,256,1),phafu),[],1); dd = dd(1:10,:);
% dd = ifft(bsxfun(@times,fft(dd,256,1),phafu),[],1); dd = dd(1:10,:);
% dd = bsxfun(@times,dd,phafu);
warning('not yet implemented');


%%  spatial reconstruction along z
if n_zpe>1
    ntmp = size(dd); ntmp(4) = mtx(3);
    ddz = complex(zeros(ntmp,class(dd)));
    if verb>0, fprintf('Spatial reconstruction along z\n'); end
    kz = -1i*2*pi*wf.kz;
    xx = reshape((-mtx(3)/2):(mtx(3)/2-1),[1 1 1 mtx(3)]);
    for lkz=1:length(kz)
        ddz = ddz + bsxfun(@times,dd(:,:,:,lkz,:,:),...
            exp(kz(1,lkz)*xx));
    end
    dd = ddz;
    clear ddz
end


%%  spatial reconstruction along xy
% size(dd) = [#spec,  #k-space one ring, 1, #z, #time, #coils]
% size(spec) = [#spec,  #x, #y, #z, #time, #coils]
spec = complex(zeros(nn));
if verb>0, fprintf('Spatial reconstruction along xy\n'); end
for lz=1:max([nz n_slices])
    for lt=1:n_time
        for lc=1:n_coils
            spec(:,:,:,lz,lt,lc) = permute(gridding(dd(:,:,1,lz,lt,lc),...
                wf.k1,wf.dcf1,mtx(1)),[3 1 2]);
        end
    end
end
clear dd


%% (Gaussian) line broadening in spectral domain
if ~isempty(lb)
    if verb>0, fprintf('Gaussian apodisation\n'); end
    lb = [lb(:) ; zeros(4,1)];
    if any(lb(2:end,1)>0), warning('Only spectral apodisation possible'); end
    if abs(lb(1))>1d-10
        lbf = lb_fun(nn(1),bw_spec,[0 lb(1)]);
        spec = bsxfun(@times,spec,lbf.');
    end
end


%% zero filling (Fourier interpolation)
if ~isempty(zf)
    if verb>0, fprintf('Zero-filling\n'); end
    zf = [zf(:).', nn((length(zf)+1):6)];
    if any(zf<nn)
        warning('zf<nn'); 
        fprintf('zf = '); disp(zf)
        fprintf('nn = '); disp(nn)
        zf = max([nn;zf],[],1);
        fprintf('Adjusting zf values\n');
        fprintf('zf = '); disp(zf)
    end
    if ((n_slices==1) && (nz==1) && (zf(4)>1))
        warning('zf(4)(=%g)>1 and 2D single slice (nz=%g && nslices=%g)',...
            zf(4),nz,n_slices);
        fprintf('  setting zf(4) to 1\n');
        zf(4) = 1;
    end
    
    % spatial zf; manually, because spatial data symmetric
    zf_tmp = zf; zf_tmp(1) = nn(1);
    spec = truma(spec,true,zf_tmp);
else
    zf = nn;
end
hz  = -(-zf(1)/2:zf(1)/2-1)/zf(1)*bw_spec;        % frequency axis
ppm = hz/h.rdb_hdr.ps_mps_freq/1d-7;


%% coil combination
if do_coco && n_coils>1
    if verb>0, fprintf('Coil combination\n'); end
    spec = csi_coil_combine(spec);
end


%% eddy current correction - phase deconvolution
if do_ecc
    if verb>0, fprintf('Eddy current correction by phase deconvolution\n'); end
    warning('must be adapted');
    spec = csi_ecc(spec);
    % dd = csi_ecc(dd,2);
end


%% HSVD water removal
if ws==1
    if verb>0
        fprintf('HSVD water subtraction: bounds=(%g %g) [Hz]\n',...
            ws_bounds(1),ws_bounds(2));
    end
    % spec = ifft(fftshift(spec,1),[],1);
    [~,~,~,spec] = hsvd(spec,25,ws_bounds,bw_spec,1);
    % spec = ifftshift(fft(spec,[],1),1);
end

%% spectral reconstruction
if verb>0, fprintf('Spectral reconstruction\n'); end
spec = ifftshift(fft(spec,zf(1),1),1);


%% Water subtraction
if ws==2
    if verb>0, fprintf('Water subtraction\n'); end
    indws = abs(hz)<30;
    spec = csi_h2osubtraction(spec,indws);
end


%% image orientation: xyz -> AP-RL-SI
if do_ori
    if verb>0, fprintf('Orienting MRSI spatially\n'); end
    spec = spec(:,end:-1:1,end:-1:1,:,:,:);
    spec = csi_orientation(spec,h);
end


%% appr. 1D spectrum
if ((nargout>2) || (verb>1)) || ~isempty(fname)
    s1d = zeros(zf(1),2,n_time);
    for lt=1:n_time
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


%% saving result as mat file
if ~isempty(fname)
    if verb>0, fprintf('Saving results (fname=''%s'')\n',fname); end
    tmp = whos('spec');
    if verb>0, fprintf('Memory size of spec = %g [MB] (SI=Base10)\n',tmp.bytes*1d-6); end
    if tmp.bytes>100d6   % reduce size for data bigger than 100MB
        if verb>0 
            fprintf('\t>100MB: converting spec (64bit) into single (32bit))\n'); 
        end
        spec = single(spec); % convert from double (64bit) to single (32bit) precision
    end
    save([fname '.mat'],'-mat','-v7.3',...
        'spec','h','zf','lb','do_ecc','ws','fname','verb',...
        'do_coco','abs_fid','cor_intlv','do_ori',...
        'hz','s1d','nn','dim',...
        'si','nx','ny','nz','n_slices','ns','n_time','n_coils');
end


%% plotting
if ((verb>1) && (ns>1)) 
    fprintf('Plotting\n');
    if length(size(spec))>5
        warning('length(size(spec))>5; skipping plotting');
    else
        if abs_fid, pc = 0; else, pc = []; end
        figstr = sprintf('P%05d Exam%d Series%d ',...
                h.image.rawrunnum,h.exam.ex_no,h.series.se_no);
        switch 3
            case 1
                clf
                plot(hz,s1d(:,:,1).');
                set(gca,'XDir','reverse');
                legend('rms','max');
            case 2
                clf
                plot_csi(spec,pc);
            case 3
                figure(20); clf
                % imagesc_row(shiftdim(max(abs(spec),[],1),1),[],'row','',true);
                imagesc_row(shiftdim(max(abs(spec),[],1),1),[],'',true,true);
                set(gcf,'name',[figstr '- image of max']);
                colormap gray
                if ~isempty(fname)
                    print('-dpng','-r600',[fname '_maximage']);
                end
                
                if ~isempty(s1d)
                    figure(21); clf
                    plot(hz,s1d(:,:,1).');
                    set(gca,'XDir','reverse');
                    legend('rms','max');
                    set(gcf,'name',[figstr '- pseudo 1D spectra']);
                    grid on
                    if ~isempty(fname)
                        print('-dpng',[fname '_s1d']);
                    end
                end
                
                if ((nx==1)||(ny==1)||((nz==1)&&(n_slices==1)))&&(n_time==1)
                    figure(22); 
                end
                if ~isempty(fname)
                    fname_plt_csi = [fname '_plot_csi'];
                else
                    fname_plt_csi = [];
                end
                plot_csi(spec,pc,[],[figstr '- CSI'],fname_plt_csi);
        end
    end
end


%% exporting images of max as dicom into scanner database
if ~isempty(fname)
    if verb, fprintf('Storing dicom images\n'); end
    inp.SeriesDescription = sprintf('CORINGS max:S%g;X%g;Y%g;Z%g,T%g,C%g', ...
        ns,nx,ny,nz+nslices-1,nt,nc);
    if dim==3, inp.MRAcquisitionType = '3D'; end
    write_dicom(fname,shiftdim(max(abs(spec),[],1),1),h,inp);
end
if false % ~isempty(fname)
    % export as dicom, if fname present fname_template.dcm is found
    template_fname = [fname '_template.dcm'];
    if ~exist(template_fname,'file')
        if verb
            fprintf('No template dicom file found\n\t%s\n',template_fname);
        end
    else
        if verb, fprintf('Storing dicom images\n'); end
        inp.SeriesNumber = h.exam.ex_numseries;
        inp.SeriesDescription = sprintf('CSI max:S%g;X%g;Y%g;Z%g,T%g,C%g', ...
            ns,nx,ny,nz+n_slices-1,n_time,n_coils);
        % inp.MRAcquisitionType = '2D';
        inp.SliceThickness = h.image.slthick;
        inp.RepetitionTime = round(h.image.tr*1d-3);
        % inp.EchoTime = round(h.rdb_hdr.user25*1d-2)*1d-1;
        % inp.ImagedNucleus = num2str(h.image.user2);
        write_dicom(fname,shiftdim(max(abs(spec),[],1),1),...
            template_fname,inp);
    end
end


%% finishing
if verb>0
    fprintf('recon_corings.m finished\n');
    fprintf('Total reconstruction time = %.2f s\n',toc(timerVal));
end

end      % main function recon_corings.m
