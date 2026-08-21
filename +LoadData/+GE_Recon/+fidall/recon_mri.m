function [bb,bbabs] = recon_mri(d,h,mtx,fname,plt,coco,lb,grdwrp)
%RECON_MRI Reconstruct standard Cartesian MRI data
%[bb,bbabs] = recon_epi(d,h,mtx,fname,plt,coco,lb,grdwrp)
%         d   Raw (p-file) data  (or pfile fname)
%         h   Header from p-file (or empty)
%       mtx   Reconstruction matrix size (Fourier interpolation)        (0)
%             1 value  => isotropic
%             2 values => 1=freq+2=phase
%             mtx(1)==0:  do nothing
%             mtx(1)==-1: round up to 2^n, minimum=128
%             mtx(1)<-1:  multiplier
%     fname   Print <fname>.png and save reco as <fname>.mat           ('')
%       plt   (1) plotting: 0=off, 1=image, 2=k-space  
%             (2) export as dicom                                   ([1 0])
%      coco   Coil combination (see mri_coil_combine.m)                 (2)
%        lb   Spatial apodisation                                       (0)
%    grdwrp   Apply gradwarp distortion correction (or gradcoil)     (true)
%
%        bb   Reconstructed data (mtx,mtx,#exc,#metabolites,#slice,#coils)
%     bbabs   RMS-coil combined images  (mtx,mtx,#exc,#metabolites,#slice)
%
%  2/2025 Rolf Schulte
if (nargin<1), help(mfilename); return; end


%% input variables
if ~exist('mtx','var'),   mtx = []; end
if isempty(mtx),          mtx = 0; end
if ~exist('fname','var'), fname = []; end
if ~exist('plt','var'),   plt = []; end
if isempty(plt),          plt = 1; end
if length(plt)<2,         plt(2) = 0; end
if ~exist('coco','var'),  coco = []; end
if isempty(coco),         coco = 2; end
if ~exist('lb','var'),    lb = []; end
if isempty(lb),           lb = 0; end
if ~exist('grdwrp','var'),grdwrp = []; end
if isempty(grdwrp),       grdwrp = true; end


%% reading in data, if pfile name is given
if ~isnumeric(d) || isempty(d)
    if ~exist(d,'file') && ~isempty(d)
        warning('strange input for d (=''%s''); file not existing?',d); 
    end
    fprintf('Reading data ''%s''\n',d);
    fname_tmp = d;
    [d,h] = read_p(d,true);
end


%% output file name
if ~isempty(fname)
    if islogical(fname)
        if fname && exist('fname_tmp','var')
            fname = fname_tmp;
        else
            fname = [];
        end
    end
end
if ~isempty(fname)
    if ~isempty(regexpi(fname,'\.7$', 'once')),  fname = fname(1:end-2); end
    if ~isempty(regexpi(fname,'\.h5$', 'once')), fname = fname(1:end-3); end
    if ~isempty(regexpi(fname,'\.mat$', 'once')),fname = fname(1:end-4); end
end


%% misc variables and checks
fprintf('size(d) = %s\n',num2str(size(d)));
[ny,nx,n3,nt,nslices,ncoils] = size(d);
if n3~=1
    warning('size(d,3)(=%g)~=1',n3);
    if nt==1
        fprintf('Swapping dim 3 and 4\n');
        d = permute(d,[1 2 4 3 5 6]);
        nt = n3;
    else
        warning('discarding data in dim 3');
        d = d(:,:,1,:,:,:);
    end
end

rot90fac = h.data_acq_tab.rotate(1);      % image rotate orientation
trnsps = h.data_acq_tab.transpose(1);     % image translate orientation
if ~((trnsps==0) || (trnsps==3))
    warning('h.data_acq_tab.transpose (=%g) not 0 or 3',trnsps);
    fprintf('\tunknown image orientation');
end
pfov = h.rdb_hdr.phase_scale;             % phase fov
if pfov>1,   warning('pfov(=%g)>1',pfov); end
if pfov<0.2, warning('pfov(=%g)<0.2',pfov); end


%% unchop
chp = ischop(h);
if ~isempty(regexpi(h.image.psdname,'3dgrass', 'once'))
    warning('3dgrass: chp = ~chp');
    chp = ~chp;
end
if chp, d(2:2:end,:) = -d(2:2:end,:); end


% 3D
if isodd(bitget(h.rdb_hdr.data_collect_type,7))
    % RHTYP3D = 64
    do3D = true;
else
    do3D = false;
end


%% reshaping data
% [YRES XRES PHASES ECHOES SLICES RECEIVERS]
% -> [#y,#x,#slices,#timesteps,#coils]
d = permute(d,[2 1 5 4 6 3]);


%% apodisation
if abs(lb)>1d-10
    lbf1 = lb_fun(ny,ny,[0 lb],'h').';
    lbf2 = lb_fun(nx,nx,[0 lb],'h');
    d = bsxfun(@times,d,lbf1);
    d = bsxfun(@times,d,lbf2);
end


%% zero-filling
if mtx(1)~=0
    if mtx(1)==-1        % rouding to 2^n, 128 minimum
        mtx = max([nx,ny]);
        mtx = 2^ceil(log2(mtx));
        if mtx<128, mtx = 128; end
    end
    if mtx(1)<-1        % multiplier
        mtx = ceil(abs(mtx(1))*max([nx,ny])/2)*2;
    end
    if isscalar(mtx), mtx(2) = mtx(1); end
    if length(mtx)>2
        warning('length(mtx)(=%d)>2; ignoring',length(mtx));
    end
    if mtx(1)<ny
        warning('mtx(1)(=%d)<ny(=%d); ignoring',mtx(1),ny); 
        mtx(1) = ny;
    end
    if mtx(2)<nx
        warning('mtx(2)(=%d)<nx(=%d)',mtx(1),nx);
        mtx(2) = nx;
    end
    if pfov<1
        mtx(1) = ceil(mtx(1)*pfov/2)*2;
    end
    mtx = [mtx(1),mtx(2),nslices,nt,ncoils];
    fprintf('Zero-filling to %s\n',num2str(mtx));
    d = truma(d,false,mtx);
else
    if pfov<1
        mtx1 = h.image.dim_Y;
    else
        mtx1 = ny;
    end
    mtx = [mtx1,nx,nslices,nt,ncoils];
end


%% half-voxel shift
d = d(:,end:-1:1,:,:,:,:);
if true
    fprintf('Half-voxel shift in %dD\n',2+do3D)
    for ld=1:(2+do3D)
        phafu = shiftdim(exp(pi*1i*((-mtx(ld)/2:mtx(ld)/2-1))/mtx(ld)).',-ld+1);
        d = bsxfun(@times,d,phafu);
    end
end


%% reconstruct data
fprintf('Reconstructing data\n');
% d = d(:,end:-1:1,:,:,:,:);
bb = fftshift(fftshift(ifft2(ifftshift(ifftshift(d,1),2)),1),2);
bb = fliplr(bb);
clear d
if do3D
    fprintf('Reconstructing 3D\n');
    % empirically found; needs more testing to confirm
    bb = fft(ifftshift(bb,3),[],3);
    bb = bb(:,:,end:-1:1,:,:,:);
end


%% add zeros if phase FOV
if pfov<1
    yy = round(h.image.dim_Y*pfov/2)*2;
    if ny~=yy, warning('ny(=%d) ~= yy(=%g)',ny,yy); end
    % if mtx(1)==0
    %     mtx1 = h.image.dim_Y;
    % else
    %     mtx1 = mtx(1);
    % end
    nn = size(bb);
    nn(1) = mtx(1);
    bb = truma(bb,false,nn);
end


%% coil combination
if ncoils>1
    bb = mri_coil_combine(bb,0,[],coco);
end


%% gradwarp distortion correction
if ischar(grdwrp)
    gradcoil = grdwrp;
    grdwrp = true;
else
    gradcoil = '';
end
if grdwrp
    try
        bb = grad_warp(bb,h,gradcoil,do3D);
    catch ME
        warning(sprintf('grad_warp failed with Error\n\t%s',ME.message));
        grdwrp = -1;
    end
end


%% image orientation
% bb = mri_orientation(bb,rot90fac,trnsps,true);
bb = mri_orientation(bb,rot90fac,trnsps);


%% calculate real images
if size(bb,5)>1
    bbabs = abs(sqrt(sum(bb.*conj(bb),5)));
else
    if ~isreal(bb), bbabs = abs(bb);
    else,           bbabs = bb; end
end


%% saving result as mat file
if ~isempty(fname)
    fprintf('Saving to file ''%s''\n',fname);
    bb = single(bb);
    bbabs = single(bbabs);
    save([fname '.mat'],'-mat','-v7.3',...
        'bb','bbabs','h','mtx','fname','plt',...
        'nx','ny','n3','nt','nslices','ncoils','grdwrp');
end


%% plotting
if (plt(1)>0) || ~isempty(fname)
    fprintf('Plotting\n');
    if (plt(1)>0)
        fid = figure(10); clf;
    else
        fid = figure('visible','off');
    end
    imagesc_row(bbabs,[],'',nt==1,true);
    % imagesc_row(bbabs,[],'gl3',nt==1,true);
    colormap gray;
    figstr = sprintf('P%05d Exam%d Series%d',...
        h.image.rawrunnum,h.exam.ex_no,h.series.se_no);
    set(fid,'name',figstr);
    drawnow
    if ~isempty(fname)
        print(fid,[fname '.png'],'-dpng','-r300');
    end
    if plt(1)==0, close(fid); end
end
if plt(1)>1
    figure(11); clf;
    imagesc_row(abs(d),[],'inm',nt==1,true);
    grid on
    set(gcf,'name',figstr);
end


%% exporting images as dicom into scanner database
if ~isempty(fname) && plt(2)>0
    fprintf('Storing dicom images\n');
    tmp = sprintf('MNSRP;MRI;SL%g,TS%g',nslices,nt);
    tmp = [tmp ';' h.series.se_desc];
    inp.SeriesDescription = tmp;
    write_dicom(fname,bbabs,h,inp,[],[],[],3);
end


end      % recon_mri.m
