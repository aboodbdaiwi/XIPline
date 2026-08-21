function snr_map=recon_csi_snr(d,h,zf,fname,plt,bw_signal,per_noise)
%RECON_CSI_SNR Calculate SNR Map from Cartesian CSI data
%   snr_map = recon_csi_snr(d,h,zf,fname,plt,bw_signal,per_noise)
%         d   Raw p-file data [(nx*ny*nz),ns,1,nt,1,nc]
%         h   Header structure
%        zf   Zero fill to [#s,#x,#y,#z]; if<0->multiplier   ([0 -2 -2 -2])
%     fname   Print <fname>.png and save reco as <fname>.mat           ([])
%             also export maps and ratios to dicom
%       plt   Plotting: 0=off,1=SNR map,2+=CSI                          (1)
% bw_signal   Integration area [Hz] default=100, 13C=600 
% per_noise   Noise percentage [spec sides,spatial]        [%]    ([30,30])
%
%  5/2024 Rolf Schulte
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
if ~exist('zf','var'),       zf = []; end
if isempty(zf),              zf = [0 -2 -2 -2]; end
if ~exist('bw_signal','var'),bw_signal = []; end
if ~exist('per_noise','var'), per_noise = []; end
if isempty(per_noise),       per_noise = 30; end
if length(per_noise)<2,      per_noise = [per_noise 30]; end
if ~exist('plt','var'),      plt = []; end
if isempty(plt),             plt = 1; end
n_maxsig = 5;                % #largest signal to average for max
lb = 10;                     % spectral apodisation [Hz]


%% misc parameters
if isdeployed      % dicom export
    do_dcm = true;
    renderer = '-painters';
else
    do_dcm = false;
    renderer = '-opengl';
end
if (plt>0) && ~isempty(fname)
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
nt = h.image.numecho;
if (nt~=1), warning('nt(=%d)~=1',nt); end


%% get nuclei specific integration boundary
if isempty(bw_signal)
    switch h.image.specnuc
        case 13         % integrate over quartett
            bw_signal = 600;
        otherwise       % integrate over singlett
            bw_signal = 100;
    end
end


%% reconstruct Cartesian CSI
[spec,hz,s1d,~,~,h] = recon_csi(d,h,zf,lb,fname_csi,plt>1);
spec0 = recon_csi(d,h,zf,0,'',false);   % for noise

clear d
[nspec,nx,ny,nz] = size(spec);
nspat = nx*ny*nz;
if size(spec,5)>1, error('size(spec,5)(=%d)>1',size(spec,5)); end
if size(spec,6)>1, error('size(spec,6)(=%d)>1',size(spec,6)); end


%% determine integration boundaries
itmp = round(per_noise(1)/200*nspec);
ind_noise = [(1:itmp),((1:itmp)+nspec-itmp)];
[~,imax] = max(s1d(:,3));
ind_signal = abs(hz-hz(imax))<bw_signal/2;


%% scaling relative to 1x1x1 cm^3 voxel size and 90° flip
slthick = h.image.slthick;        % slice thickness [mm]
fov  = h.rdb_hdr.fov;              % FOV [mm]
xcsi = h.rdb_hdr.xcsi;
ycsi = h.rdb_hdr.ycsi;
zcsi = h.rdb_hdr.zcsi;
fa   = h.image.mr_flip;
nex  = h.image.nex;
if zcsi>1     % 3D
    res = fov./[xcsi ycsi zcsi];
else          % 2D
    res = [fov./xcsi fov./ycsi slthick];
end
scale_vox = 1000/prod(res);                      % (1cm)^3 voxel size
scale_fa  = 1/sind(fa);                          % 90° flip angle
scale_ave = sqrt(1000/(xcsi*ycsi*zcsi*nex*nt));  % 10x10x10
scale = scale_vox*scale_fa*scale_ave;            % overall scaling
fprintf('scaling = %g (spatial) * %g (flip) * %g (#exc) = %g\n',...
    scale_vox,scale_fa,scale_ave,scale);


%% quantify: using root-sum-squares: better for SNR
% signal = sum(abs(spec(ind_signal,:,:,:)),1);
signal = sqrt(sum(abs(spec(ind_signal,:,:,:).*conj(spec(ind_signal,:,:,:))),1));

svec = signal(:).';
[~,ii] = sort(svec,'descend');
ii = ii(1,1:n_maxsig);
signal_max = mean(svec(1,ii));
% signal_max = max(signal(:));

noise = spec0(ind_noise,:,:,:);
noise_val = abs(std(noise(:)));
fprintf('max(signal) = %g; noise_val = %g; max(snr) = %g\n',...
    signal_max,noise_val,signal_max/noise_val);
fprintf('scaled max(snr) = %g\n',scale*signal_max/noise_val);
snr_map = scale*reshape(signal./noise_val,[nx,ny,nz]);


%% saving data
if ~isempty(fname)
    save([fname '_snr_map.mat'],'-mat','-v7.3','snr_map','h',...
        'zf','bw_signal','per_noise','scale');
end


%% exporting as dicom
if ~isempty(fname) && do_dcm
    fprintf('Storing dicom images\n');
    inp = struct;
    if nz>1, inp.MRAcquisitionType = '3D'; end
    write_dicom([fname '_snr_map'],snr_map,h,inp);
end


%% plotting ratios
if plt>0
    figstr = sprintf('P%05d Exam%d Series%d - SNR map',...
        h.image.rawrunnum,h.exam.ex_no,h.series.se_no);

    figure(22); clf
    set(gcf,'name',figstr);
    if zcsi==1
        imagesc(snr_map);
        colorbar; colormap jet
        if ~isempty(fname)
            print([fname '_snr.png'],'-dpng','-r300',renderer);
            if do_dcm>0
                inp.SeriesNumber = h.series.se_no+100;
                write_scdicom([fname '_snr.dcm'],gcf,h,inp);
            end
        end
    else
        imagesc_row(snr_map,'','',true,true);
        colorbar; colormap jet
        if ~isempty(fname)
            print([fname '_snr_row.png'],'-dpng','-r300',renderer);
            if do_dcm>0
                inp.SeriesNumber = h.series.se_no+100;
                write_scdicom([fname '_snr_row.dcm'],gcf,h,inp);
            end
        end

    
        figure(23); clf
        set(gcf,'name',figstr);
        imagesc_ind3d(snr_map);
        colorbar; colormap jet
        if ~isempty(fname)
            print([fname '_snr_ind3d.png'],'-dpng','-r300',renderer);
            if do_dcm>0, write_scdicom([fname '_snr_ind3d.dcm'],gcf,h,inp); end
        end
    end

    figure(24); clf
    set(gcf,'name',figstr);
    subplot(2,1,1);
    plot(hz(ind_signal),reshape(abs(spec(ind_signal,:,:,:)),[sum(ind_signal) nspat]));
    grid on; title('signal'); axis tight;

    subplot(2,1,2); 
    plot(reshape(abs(noise),[size(noise,1) nspat]));
    grid on; title('noise'); axis tight;

end


end      % recon_csi_snr.m
