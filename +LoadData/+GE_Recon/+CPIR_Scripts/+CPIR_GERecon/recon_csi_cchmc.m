function [spec,hz,s1d,ppm,par,h,bb] = recon_csi(d,h,zf,lb,fname,plt,coco,...
    ws,do_ecc,d_ref)
%RECON_CSI  Reconstruct Cartesian CSI product data
%[spec,hz,s1d,ppm,par,h]=recon_csi(d,h,zf,lb,fname,plt,coco,ws,do_ecc,d_ref)
%      d  Raw p-file data [(nx*ny*nz),nspec,1,nt,1,nc]
%         Interleaved reference scan: nt=1 (1=mrs,2=h2o)
%      h  Header structure
%     zf  Zero fill to:   [#s,#x,#y,#z]; if<0->multiplier
%     lb  Line broadening [s,x,y,z]
%  fname  Print <fname>.png and save reco as <fname>.mat               ([]) 
%         also export max-image and plots as dicom
%    plt  Plotting: 0=off, 1=normal, 2+=plot_mrsi                       (2)
%   coco  Coil combination [do percent_noise4decor]                 ([1 0])
%     ws  Water suppression: 0=off; 1=HSVD; 2=subtraction               (0)
% do_ecc  Phase deconvolution eddy current correction               (false)
%         requires H2O reference scan: 
%         if ((nt==2)&&(h.rdb_hdr.user23==1)) -> default=true
%         if ~isempty(d_ref) -> default=true; 
%  d_ref  Reference scan (w/o water suppression)
%
%   spec  Spectrum        [#s,#x,#y,#z,#t]
%     hz  Frequency axis  [Hz]
%    s1d  pseudo spectra  [#s,3]
%    ppm  Frequency axis  [ppm]
%    par  Parameter structure for plot_mrsi
%
%  6/2024 Rolf Schulte
if nargin<1, help(mfilename); return; end

timerVal = tic;              % record execution time

%% misc input
if ~exist('zf','var'),     zf = []; end
if ~exist('lb','var'),     lb = []; end
if ~exist('fname','var'),  fname = []; end
if islogical(fname)
    if fname, error('provide fname not true');
    else, fname = [];
    end
end
if ~isempty(fname)
    if ~islogical(fname)
        if ~isempty(regexpi(fname,'\.7$',  'once')), fname = fname(1:end-2); end
        if ~isempty(regexpi(fname,'\.h5$', 'once')), fname = fname(1:end-3); end
        if ~isempty(regexpi(fname,'\.mat$','once')), fname = fname(1:end-4); end
    end
end
if ~exist('plt','var'),    plt = []; end
if isempty(plt),           plt = 2; end
if length(plt)~=1,         warning('length(plt)~=1'); end
if ~exist('coco','var'),   coco = []; end
if isempty(coco),          coco = 1; end
if length(coco)<2,         coco = [coco 0]; end
if ~exist('ws','var'),     ws = []; end
if isempty(ws),            ws = 0;  end
if length(ws)~=1,          warning('length(ws)~=1'); end
if ~exist('do_ecc','var'), do_ecc = []; end
if ~exist('d_ref','var'),  d_ref = []; end
do_ori = 1;          % correct image orientation


%% reading in data, if pfile name is given
if ~isnumeric(d) || isempty(d)
    if ~exist(d,'file') && ~isempty(d)
        warning('strange input for d (=''%s''); file not existing?',d); 
    end
    [d,h] = read_p(d);
end
d = conj(d);                 % changing dir of axis (keep MRS convention)


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
    tmp = complex(zeros(size(d,1),size(d,2),size(d,3),2,size(d,5),...
        size(d,6),class(d)));
    tmp(:,:,:,1,:,:) = d;
    tmp(:,:,:,2,:,:) = d_ref;
    d = tmp;
    clear tmp d_ref
    nt = 2;
end


%% misc parameters
ws_bounds = [-100 50];                 % bounds for HSVD water suppression
si = size(d); si = [si ones(1,6-length(si))];
nx = h.rdb_hdr.xcsi;
ny = h.rdb_hdr.ycsi;
nz = h.rdb_hdr.zcsi;
nspec = h.rdb_hdr.user1;               % #sampled points (FID)
nslices = h.rdb_hdr.nslices;
if ~exist('nt','var'), nt = h.image.numecho; end  % number of time steps
nc = (h.rdb_hdr.dab(2)-h.rdb_hdr.dab(1))+1; % number of Rx coil elements
nref = h.rdb_hdr.user19;               % #reference frames (w/o WS)


%% misc checks
if si(2)~=nspec, warning('si(2)(=%g)~=nspec(=%g)',si(2),nspec); end
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
    case 2, nn = [nspec nx ny nslices nt nc];
    case 3, nn = [nspec nx ny nz nt nc];
    otherwise, error('dim(=%g) not 2 or 3',dim);
end
bw = h.rdb_hdr.spectral_width;              % sampling bandwidth [Hz]
t  = (0:nspec-1).'/bw;                      % time [s]
zf_pseuim = zf(5:end);
zf = zf(1:min(4,length(zf)));
[zf,do_zf_spat] = sub_zf(zf,nn,nz,nslices);


%% reference scan for eddy current compensation
if ((do_ecc) && (nref>0))
    if nt~=1, warning('nt(=%g)~=1',nt); end
    fprintf('Using reference scan for Eddy Current Compensation\n');
    tmp = complex(zeros(si(1)-nref,si(2),si(3),2,si(5),si(6),class(d)));
    tmp(:,:,:,1,:,:) = d(1:(si(1)-nref),:,:,1,:,:);
    tmp(1,:,:,2,:,:) = d((si(1)-nref+1),:,:,1,:,:);
    d = tmp;
    nt = 2;
    nn(5) = 2;
end


%% misc info + warnings
fprintf('psdname=%s/%s\n',h.image.psdname,h.image.psd_iname);
fprintf('nspec=%g, nx=%g, ny=%g, nz=%g, nslices=%g, nt=%g, nc=%g\n',...
    nspec,nx,ny,nz,nslices,nt,nc);
fprintf('nucleus=%g\n',h.image.specnuc);
if nx*ny*nz~=si(1), warning('nx*ny*nz~=si(1)'); end
if nspec~=si(2),    warning('nspec(=%g)~=si(2)(=%g)',nspec,si(2)); end
if nspec<64,        warning('nspec(=%g)<64',nspec); end

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
if ischop(h), d(2:2:end,:) = -d(2:2:end,:); end


%% reorder data
fprintf('Reordering data\n'); 
% size(d)  = [#x*#y*#z, #s,  1, #t,  1, #c]
% size(dd) = [#s,       #x, #y, #z, #t, #c]
dd = complex(zeros(nn,class(d)));
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


%% image orientation: xyz -> AP-RL-SI
if do_ori==2
    dd = csi_orientation(dd,h.data_acq_tab.rotate,...
        h.data_acq_tab.transpose,false);
end


%% (Gaussian) line broadening - both spatial and spectral
if ~isempty(lb)
    fprintf('Gaussian apodisation\n');
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
if do_zf_spat
    fprintf('Zero-filling (spatial)\n');
    % spatial zf; manually, because spatial data symmetric
    zf_tmp = zf; zf_tmp(1) = nn(1);
    dd = truma(dd,false,zf_tmp);
end


%% frequency axes
hz  = (-zf(1)/2:zf(1)/2-1)/zf(1)*bw;
ppm = hz/h.rdb_hdr.ps_mps_freq/1d-7;


%%  spatial reconstruction
for l=2:(dim+1)
    if zf(l)>1
        fprintf('Spatial reconstruction: dim=%g\n',l);
        % dd = ifftshift(ifft(fftshift(dd,l),[],l),l);
        dd = ifftshift(fft(fftshift(dd,l),[],l),l);
    end
end


%% coil combination
if coco(1)>0 && nc>1
    fprintf('Coil combination\n');
    dd = csi_coil_combine(dd,[],[],[],coco(2));
end


%% eddy current correction - phase deconvolution
if do_ecc
    fprintf('Eddy current correction by phase deconvolution\n');
    dd = csi_ecc(dd);
end


%% water subtraction via HSVD
if ws==1
    fprintf('HSVD water subtraction: bounds=(%g %g) [Hz]\n',...
        ws_bounds(1),ws_bounds(2));
    [~,~,~,dd(:,:,:,:,1,:)] = hsvd(dd(:,:,:,:,1,:),25,ws_bounds,bw,1);
end


%% spectral reconstruction
% size(spec) = [#s,#x,#y,#z,#t]
% fft recon
if zf(1)>1
    fprintf('Spectral reconstruction\n');
    spec = ifftshift(fft(dd,zf(1),1),1);
else
    warning('Skipping spectral reconstruction\n');
end


%% Water subtraction
if ws==2
    fprintf('Water subtraction\n');
    indws = abs(hz)<30;
    spec = csi_h2osubtraction(spec,indws,true);
end


%% image orientation: xyz -> AP-RL-SI
if do_ori==1
    spec = csi_orientation(spec,h.data_acq_tab.rotate,...
        h.data_acq_tab.transpose,true);
end


%% appr. 1D spectrum
if ((nargout>2) || (plt>0)) || ~isempty(fname)
    % size(spec)=[#s,#x,#y,#z,#t,#c]
    s1d = zeros(zf(1),3);
    sss = sum(spec,5);
    s1d(:,3) = sqrt(sum(sum(sum(sum(abs(sss.*conj(sss)),6),4),3),2));
    sss = sqrt(sum(abs(sss.*conj(sss)),6));
    s1d(:,1) = sum(sum(sum(sss,2),3),4);
    s1d(:,2) = max(max(max(sss,[],2),[],3),[],4);
    s1d = bsxfun(@rdivide,s1d,max(s1d,[],1));
end
figstr = sprintf('P%05d Exam%d Series%d ',...
    h.image.rawrunnum,h.exam.ex_no,h.series.se_no);
par.fax = hz; par.fig_str = figstr; par.lb = lb; % par structure for plot_mrsi
par.te = h.rdb_hdr.te*1d-6;  % echo time [s]


%% saving result as mat file
bb = mrsi2pseuimage(spec,zf_pseuim);
if ~isempty(fname)
    fprintf('Saving results to fname=''%s''\n',fname);
    spec = single(spec); % convert: double (64bit) to single (32bit) precision
    
    save([fname '.mat'],'-mat','-v7.3',...
        'spec','h','zf','lb','do_ecc','ws','fname','plt',...
        'coco','do_ori',...
        'hz','s1d','nn','t','dim','bb',...
        'si','nx','ny','nz','nslices','nspec','nt','nc','par');
end



%% plot check
if length(size(spec))>5
    warning('(length(size(spec))>5); skipping plotting');
    plt = 0;
    fname = [];
end


%% plotting pseudo image
if (plt>0) || ~isempty(fname)
    if plt>0
        fprintf('Plotting pseudo image (visible)\n');
        fid = figure(20); clf
    else
        fprintf('Plotting pseudo image (invisible; saving pix)\n');
        fid = figure('Visible','off');
    end
    % bs = sort(bb(:),'descend');
    % mm = 0.97*mean(bs(1:ceil(0.02*length(bs))));
    % imagesc_row(bb,[0 mm],'',size(bb,4)==1,true);
    mini=prctile(bb(:),1);
    maxi=prctile(bb(:),99.9);
    imagesc_row(bb,[mini, maxi],'',size(bb,4)==1,true);
    set(fid,'name',[figstr '- pseudo image']);
    colormap gray
    if ~isempty(fname)
        if isdeployed
            print(fid,[fname '_pseuimage'],'-dpng','-r300','-painters');
            % export to dicom database
            inp.SeriesDescription = sprintf('CSI pseudo:S%g;X%g;Y%g;Z%g,T%g,C%g', ...
                nspec,nx,ny,nz+nslices-1,nt,nc);
            if dim==3
                inp.MRAcquisitionType = '3D';
                inp.slthick = h.image.slthick/zf(4);
                inp.SpacingBetweenSlices = inp.slthick;
                if abs(h.series.start_loc-h.series.end_loc)<0.01
                    loc = h.series.start_loc;
                    slthick = h.image.slthick;
                    h.series.start_loc = loc - slthick/2;
                    h.series.end_loc = loc + slthick/2;
                end
            end
            write_dicom_cchmc(fname,bb,h,inp,[],[],[],3);
        else
            print(fid,[fname '_pseuimage'],'-dpng','-r300');
        end
    end
    if plt<=0, close(fid); end
end


%% pseudo spectrum
if ((plt>0) || ~isempty(fname)) && ~isempty(s1d)
    if plt>0
        fprintf('Plotting pseudo spectrum (visible)\n');
        fid = figure(21); clf
    else
        fprintf('Plotting pseudo spectrum (invisible; saving pix)\n');
        fid = figure('Visible','off');
    end
    plot(hz,s1d(:,:,1).');
    set(gca,'XDir','reverse');
    legend('sum','max','rms');
    set(fid,'name',[figstr '- pseudo 1D spectra']);
    axis([min(hz) max(hz) -0.01 1.02]);
    grid on
    xlabel('freq [Hz]');
    if ~isempty(fname)
        if isdeployed
            print(fid,[fname '_s1d'],'-dpng','-r300','-painters');
            % export to dicom database
            inpscd.SeriesNumber = h.series.se_no+100;
            inpscd.SeriesDescription = ...
                sprintf('CSI Se%d',h.series.se_no);
            inpscd.InstanceNumber = 0;
            write_scdicom([fname '_s1d.dcm'],fid,h,inpscd);
        else
            print(fid,[fname '_s1d'],'-dpng','-r300');
        end
    end
    if plt<=0, close(fid); end
end


%% plotting csi grid (invisible + dicom)
if ~isempty(fname)
    if isdeployed, hh = h;
    else, hh = []; end
    plot_csi(spec,[],[],[figstr '- CSI'],[fname '_plot_csi'],11,hh);
end


%% plotting 3-plane of pseudo images
if ((dim>2) && ((plt>0) || ~isempty(fname)))
    % bb3 = mrsi2pseuimage(spec,[128 128 128]);
    if isempty(zf_pseuim), zf_pseuim = [128 128 128]; end
    bb3 = mrsi2pseuimage(spec,zf_pseuim);
    if plt>0
        fprintf('Plotting pseudo 3-plane (visible)\n');
        fid = figure(22); clf
    else
        fprintf('Plotting pseudo 3-plane (invisible; saving pix)\n');
        fid = figure('Visible','off');
    end
    % bs = sort(bb3(:),'descend');
    % mm = 0.96*mean(bs(1:20d3));
    mini=prctile(bb(:),1);
    maxi=prctile(bb(:),99.9);
    % imagesc_ind3d(bb3,[],[0 mm],[],false,true);
    imagesc_ind3d(bb3,[],[mini, maxi],[],false,true);
    set(fid,'name',[figstr '- pseudo 3-planee']);
    colormap gray
    if ~isempty(fname)
        if isdeployed
            print([fname '_pseu3plane'],'-dpng','-r300','-painters');
        else
            print([fname '_pseu3plane'],'-dpng','-r300');
        end
    end
    if plt<=0, close(fid); end
end


%% plotting interactive MRSI viewer
if plt>1
    fprintf('Starting plot_mrsi\n');
    plot_mrsi(spec,par);
end


%% finishing up
fprintf('recon_csi: runtime = %.2f s\n',toc(timerVal));

end      % recon_csi.m


%% sub-functions
% zero filling
function [zf,do_zf_spat] = sub_zf(zf,nn,nz,nslices)
do_zf_spat = false;
if ~isempty(zf)
    zf = [zf(:).', nn((length(zf)+1):6)];
    if any(zf<0)             % negative numbers -> factor
        ii = (zf<0) & (nn>1);
        zf(ii) = -zf(ii).*nn(ii);
        zf(nn==1) = 1;
    end
    if any(abs(zf)<1d-10)    % 0 -> no zero-filling
        ii = abs(zf)<1d-10;
        zf(ii) = nn(ii);
    end
    if any(zf<nn)            % no cropping
        warning('zf<nn');
        fprintf('zf = %s\n',num2str(zf));
        fprintf('nn = %s\n',num2str(nn));
        zf = max([nn;zf],[],1);
        fprintf('Adjusting zf values\n');
        fprintf('zf = %s\n',num2str(zf));
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
    zf = round(zf);
    
    if any(zf(2:4)~=nn(2:4)), do_zf_spat = true; end

%     % spatial zf; manually, because spatial data symmetric
%     tmp = zeros(nn(1),zf(2),zf(3),zf(4),nt,nc);
%     for l=2:4
%         if isodd(zf(l)-nn(l))
%             warning('zf(%g)(=%g) - nn(%g)(=%g) = %g; rounding',l,zf(l),l,nn(l),zf(l)-nn(l)); 
%         end
%         x{l} = ceil(((zf(l)-nn(l))/2+1):((zf(l)+nn(l))/2)); 
%     end
%     tmp(:,x{2},x{3},x{4},:,:) = dd;
%     dd = tmp;
else
    zf = nn;
end

end      % sub_zf
