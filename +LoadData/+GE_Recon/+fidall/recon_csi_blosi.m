function b1map=recon_csi_blosi(d,h,zf,lb,fname,plt,n_pts,std_thresh)
%RECON_CSI_BLOSI Calculate B1+ Map from Cartesian CSI data
%  b1map = recon_csi_blosi(d,h,zf,lb,fname,plt,n_pts)
%      d   Raw p-file data [(nx*ny*nz),ns,1,nt,1,nc]
%      h   Header structure
%     zf   Zero fill to [#s,#x,#y,#z]; if<0->multiplier      ([0 -2 -2 -2])
%     lb   Line broadening [s,x,y,z]                                   ([])
%  fname   Print <fname>.png and save reco as <fname>.mat              ([])
%          also export maps and ratios to dicom
%    plt   Plotting: 0=off,1=B1+ map,2+=CSI                             (1)
%  n_pts   #largest FID points for Bloch-Siegert phase                 (30)
%std_thresh Mask B1+ map to std(blosi phase)<std_thresh         [deg]  (60)
%
% 11/2024 Rolf Schulte
if (nargin<1), help(mfilename); return; end


%% input parameters
if ~exist('fname','var'), fname = []; end
if islogical(fname)
    if fname, error('provide fname not true');
    else, fname = [];
    end
end
if ~isempty(fname)
    if ~isempty(regexpi(fname,'\.7$',  'once')), fname = fname(1:end-2); end
    if ~isempty(regexpi(fname,'\.h5$', 'once')), fname = fname(1:end-3); end
    if ~isempty(regexpi(fname,'\.mat$','once')), fname = fname(1:end-4); end
end
if ~exist('zf','var'),    zf = []; end
if isempty(zf),           zf = [0 -2 -2 -2]; end
if ~exist('lb','var'),    lb = []; end
if ~exist('plt','var'),   plt = []; end
if isempty(plt),          plt = 1; end
if ~exist('n_pts','var'), n_pts = []; end
if isempty(n_pts),        n_pts = 30; end
if ~exist('std_thresh','var'), std_thresh = []; end
if isempty(std_thresh),   std_thresh = 60; end


%% misc parameters
if isdeployed      % dicom export
    do_dcm = true;
else
    do_dcm = false;
end
if (plt>1) && ~isempty(fname)
    fname_csi = [fname '_csi'];
else
    fname_csi = [];         % don't save MRSI
end


%% reading in data, if pfile name is given
if ~isnumeric(d) || isempty(d)
    if ~exist(d,'file') && ~isempty(d)
        warning('strange input for d (=''%s''); file not existing?',d); 
    end
    if isempty(d) || ~ischar(d), d = ''; end
    [d,h] = read_p(d);
end


%% more parameters
nblosi = h.rdb_hdr.user19;             % # of blosi off-freq pulses/sign
xcsi = h.rdb_hdr.xcsi;
ycsi = h.rdb_hdr.ycsi;
zcsi = h.rdb_hdr.zcsi;
if nblosi~=(xcsi*ycsi*zcsi)
    warning('nblosi(=%d) ~= (xcsi(=%d)*ycsi(=%d)*zcsi(=%d)',...
        nblosi,xcsi,ycsi,zcsi);
end
necho = h.image.numecho;
if (necho<2), error('necho(=%d)<2',necho); end


%% reconstruct Cartesian CSI
[spec,~,s1d,~,~,h] = recon_csi(d,h,zf,lb,fname_csi,plt>1);
clear d
[nspec,nx,ny,nz,nt] = size(spec);
fprintf('nspec=%d, nx=%d, ny=%d, nz=%d, nt=%d\n',nspec,nx,ny,nz,nt);
% nspat = nx*ny*nz;
if nt~=necho, warning('nt(=%d)~=necho(=%d)',nt,necho); end
if size(spec,6)>1,  error('size(spec,6)(=%d)>1',size(spec,6)); end


%% IFFT back to FID
fid = ifft(fftshift(spec,1),[],1);
fid = conj(fid);


%% index to largest signal
f1d = ifft(fftshift(s1d(:,2,1)));
% f1d = s1d(:,2,1);
ii = [true(floor(nspec/2),1) ; false(ceil(nspec/2),1)];
ii(1:5,1) = false;
f1d(~ii) = 0;
[~,ind] = sort(f1d,'descend');
ind = ind(1:n_pts,1);


%% blosi parameters
KBS = header2kbs(h);
if nt~=2
    fprintf('nt(=%d)~=2: averaging blosi\n',nt);
    tmp = complex(zeros(nspec,nx,ny,nz,2,class(spec)));
    tmp(:,:,:,:,1) = mean(fid(:,:,:,:,1:2:end),5);
    tmp(:,:,:,:,2) = mean(fid(:,:,:,:,2:2:end),5);
    fid = tmp;
    clear tmp;
end
blosi_phase = angle(conj(fid(:,:,:,:,2)).*fid(:,:,:,:,1));
if ~isreal(blosi_phase), warning('~isreal(blosi_phase)'); end
ii = blosi_phase<-pi/2;
blosi_phase(ii) = blosi_phase(ii)+2*pi;
bs_mean = reshape(mean(blosi_phase(ind,:,:,:),1),[nx ny nz]);
b1map = real(sqrt(bs_mean/2/KBS));
bs_std = reshape(std(blosi_phase(ind,:,:,:),[],1),[nx ny nz]);
if std_thresh>0
    % mask = (bs_std*180/pi<std_thresh) & (bs_mean>0);
    mask = (bs_std*180/pi<std_thresh);
    b1map = b1map.*mask;
end
if any(isnan(b1map(:)))
    warning('isnan(b1map): setting to zero');
    b1map(isnan(b1map)) = 0;
end
if any(isinf(b1map(:)))
    warning('isinf(b1map): setting to zero');
    b1map(isinf(b1map)) = 0;
end


%% saving data
if ~isempty(fname)
    save([fname '_b1map.mat'],'-mat','-v7.3','b1map','h','zf','n_pts');
end


%% exporting as dicom
if ~isempty(fname) && do_dcm
    fprintf('Storing dicom images\n');
    inp = struct;
    if nz>1, inp.MRAcquisitionType = '3D'; end
    write_dicom([fname '_b1map'],b1map,h,inp);
end


%% plotting ratios
if plt>0
    figstr = sprintf('P%05d Exam%d Series%d - ',...
        h.image.rawrunnum,h.exam.ex_no,h.series.se_no);

    if nz==1
        figure(22); clf; set(gcf,'name',[figstr 'B1+ map']);
        imagesc(b1map);
        title('B1+ [uT]','FontSize',14);
        colorbar; colormap jet
        if ~isempty(fname)
            print([fname '_b1map.png'],'-dpng','-r300');
            if do_dcm>0
                inp.SeriesNumber = h.series.se_no+100;
                write_scdicom([fname '_b1map.dcm'],gcf,h,inp);
            end
        end

        figure(23); clf; set(gcf,'name',[figstr 'std(blosi phase)']);
        imagesc(bs_std*180/pi,[0 180]);
        title('std(blosi phase) [deg]','FontSize',14);
        colorbar; colormap jet

    else
        figure(22); clf; set(gcf,'name',[figstr 'B1+ map']);
        imagesc_ind3d(b1map);
        colorbar; colormap jet
        if ~isempty(fname)
            print([fname '_b1map.png'],'-dpng','-r300');
            if do_dcm>0, write_scdicom([fname '_b1map.dcm'],gcf,h,inp); end
        end

        figure(23); clf; set(gcf,'name',[figstr 'std(blosi phase)']);
        imagesc_ind3d(bs_std*180/pi,[],[0 180]);
        colorbar; colormap jet

        figure(24); clf; set(gcf,'name',[figstr 'mean(blosi phase)']);
        imagesc_ind3d(bs_mean*180/pi,[],[-1 1]*mean(bs_mean(:))*180/pi);
        colorbar; colormap jet
    end
end


end      % recon_csi_blosi.m
