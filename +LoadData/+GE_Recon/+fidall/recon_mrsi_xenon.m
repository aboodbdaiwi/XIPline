function bb=recon_mrsi_xenon(d,h,wfn,fname,plt_mrsi,zf,freqs,nint,...
    mean_frac,thresh_fac)
%RECON_MRSI_XENON Reconstruct and quantify dissolved-phase xenon MRSI
%      bb = recon_mrsi_xenon(d,h,wfn,fname,plt_mrsi,zf,freqs,nint,...
%              mean_frac,thresh_fac)
%       d   Raw p-file data [(nx*ny*nz),ns,1,nt,1,nc]
%       h   Header structure
%     wfn   Location of waveform .mat file
%   fname   Print <fname>.png and save reco as <fname>.mat             ([])
%           also export maps and ratios to dicom
%plt_mrsi   Plot and save MRSI data                                     (0)
%      zf   Zero fill to [#s,#x,#y,#z]; if<0->multiplier   ([256 -2 -2 -2])
%   freqs   Xenon frequencies (f0 relative to gas) [ppm]      ([0 197 217])
%    nint   Integration boundaries                 [pts] ([-5 5;-1 2;-2 1])
%mean_frac  Fraction of largest gas signal for threshold              (0.2)
%thresh_fac Factor for threshold                                      (0.6)
%
%  9/2024 Rolf Schulte
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
if ~exist('wfn','var'),      wfn = []; end
if ~exist('plt_mrsi','var'), plt_mrsi = []; end
if isempty(plt_mrsi),        plt_mrsi = 0; end
if ~exist('zf','var'),       zf = []; end
if isempty(zf),              zf = [256 -2 -2 -2]; end
if ~exist('freqs','var'),    freqs = []; end
if isempty(freqs),           freqs = [0 197 217]; end
if ~exist('nint','var'),     nint = []; end
if isempty(nint),            nint = [-4 4;-2 4;-4 2]; end
if ~exist('mean_frac','var'),mean_frac = []; end
if isempty(mean_frac),       mean_frac = 0.2; end
if ~exist('thresh_fac','var'), thresh_fac = []; end
if isempty(thresh_fac),      thresh_fac = 0.6; end
slowft = false;         % use FFT


%% misc parameters
plt = 1;
bb_str = {'gas','tissue','blood','tissue/gas ratio','blood/gas ratio',...
    'blood/tissue ratio'};
rat = [2 1;3 1;3 2];         % calculate ratios
ratsca = [0.015 7d-3 0.8];
fff = {'tg_rat','bg_rat','bt_rat'};

if (plt_mrsi>0) && ~isempty(fname)
    fname_mrsi = [fname '_mrsi'];
else
    fname_mrsi = [];         % don't save MRSI
end


%% reading in data, if pfile name is given
if ~isnumeric(d)
    if exist(d,'file')
        [d,h] = read_p(d);
    else
        warning('strange input d/file not existing');
    end
end


%% plot first point of FID
plt_fid = 0;
if plt>0
    if ~isempty(regexpi(h.image.psd_iname,'fidall',  'once'))
        plt_fid = 1;
    end
    if ~isempty(regexpi(h.image.psd_iname,'fidcsi',  'once'))
        if h.image.user19==0
            plt_fid = 1;
        else
            if isfield(h.data_acq_tab,'view_acq')
                plt_fid = 2;
            end
        end
    end
end
if plt_fid>0
    figstr = sprintf('P%05d Exam%d Series%d ',...
        h.image.rawrunnum,h.exam.ex_no,h.series.se_no);
    figure(29); clf
    set(gcf,'name',[figstr '- FID: first point']);
    fid = sqrt(mean(mean(mean(mean(abs(d.*conj(d)),3),4),5),6));
    fid = mean(fid(:,1:2),2);
    if plt_fid>1   % data not acquired sequentially; reorder
        fid = fid(h.data_acq_tab.view_acq,1);
    end
    fid0 = fid(1,1);
    if abs(fid0)<0.1
        warning('fid(1,1)(=%g)==0 -> using max for normalisation',fid0);
        fid0 = max(fid);
    end
    fid = fid/fid0*100;
    plot(fid)
    % title('FID: first point');
    xlabel('excitations');
    ylabel('signal [%]');
    axis([0 size(fid,1) 0 max(fid)])
    grid on
    drawnow
    if ~isempty(fname)
        if isdeployed
            print([fname '_fid.png'],'-dpng','-r300','-painters');
            write_scdicom([fname '_fid.dcm'],gcf,h);
        else
            print([fname '_fid.png'],'-dpng','-r300');
        end
    end
end


%% reconstruct product CSI or fidall MRSI
if (h.rdb_hdr.xcsi==1) && (h.rdb_hdr.ycsi==1) && (h.rdb_hdr.zcsi==1)
    % [spec,hz,s1d,ppm,par,h] = recon_mrsi(d,h,wfn,zf,'','','',fname_mrsi,...
    %     1+plt_mrsi,[],[],[],slowft);
    [spec,hz,s1d,ppm,par,h] = recon_mrsi(d,h,wfn,zf,'',fname_mrsi,...
        plt_mrsi,[],slowft);
else
    if ~isempty(wfn)
        warning('x,y,zcsi(=%d,%d,%d)>1 && ~isempty(wfn)',...
            h.rdb_hdr.xcsi,h.rdb_hdr.ycsi,h.rdb_hdr.zcsi);
    end
    [spec,hz,s1d,ppm,par,h] = recon_csi(d,h,zf,'',fname_mrsi,plt_mrsi);
end
[~,nx,ny,nz] = size(spec);
nn = nx*ny*nz;


%% quantify
[area,iint] = integrate_mrsi(spec,ppm,freqs,nint);


%% calculate mask
bs = area(1,:,:,:); bs = sort(bs(:));
ii = round(nn*(1-mean_frac)):nn;
thresh = thresh_fac*mean(bs(ii));
mask = area(1,:,:,:)>thresh;
mc = circular_mask([nx,ny,nz],[],max([nx,ny,nz])/2);
% mask = mask & permute(mc,[4 1 2 3]);
mask = permute(mask,[2 3 4 1]) & mc;


%% adjust signal scaling for gas image
te = par.te;
if h.rdb_hdr.ps_mps_freq<20d7
    b0str = '1.5T';
    T2s = [17 2.2 2]*1d-3;        % gas-tissue-blood
else
    b0str = '3T';
    T2s = [8.7 1.3 1.2]*1d-3;     % gas-tissue-blood
end
bb_scale = exp(te./T2s);
mr_flip = h.image.mr_flip;
if mr_flip<1
    warning('mr_flip(=%g)<1; setting to 10',mr_flip);
    mr_flip = 10;
end
bb_scale(1,1) = sind(mr_flip)/sind(mr_flip/100)*...
    bb_scale(1,1);
fprintf('%s: T2*=%g;%g;%g[ms]; flip=%g;%g;%g[deg] - gas;tissue;blood\n',...
    b0str,T2s*1d3,mr_flip*[0.01 1 1]);
fprintf('    bb_scale=%g;%g;%g\n',bb_scale);


%% calculate ratio maps
nfreq = length(freqs);
nrat = size(rat,1);
bb = zeros(nx,ny,nz,nfreq+nrat);
bb(:,:,:,1:nfreq) = permute(area,[2 3 4 1]);
for l=1:3
    bb(:,:,:,l) = bb_scale(l)*bb(:,:,:,l);   % scale images
end
for lr=1:nrat
    bb(:,:,:,lr+nfreq) = bb(:,:,:,rat(lr,1))./bb(:,:,:,rat(lr,2)).*mask;
end


%% plotting ratios
if plt>0
    figstr = sprintf('P%05d Exam%d Series%d ',...
        h.image.rawrunnum,h.exam.ex_no,h.series.se_no);
    
    figure(24); clf
    set(gcf,'name',[figstr '- gas,tissue,blood,t/g,b/g,b/t']);
    imagesc_row(bb,'','row');
    colormap jet
    if ~isempty(fname)
        if isdeployed
            print([fname '_bb.png'],'-dpng','-r300','-painters');
            write_scdicom([fname '_bb.dcm'],gcf,h);
        else
            print([fname '_bb.png'],'-dpng','-r300');
        end
    end
    
    figure(25); clf
    set(gcf,'name',[figstr '- spec with integ-boundaries']);
    x1 = ppm(iint(:,1));
    x2 = ppm(iint(:,2));
    xx = [x1;x1;x2;x2;x1];
    yy = [0;1;1;0;0]*ones(1,size(iint,1));
    plot(ppm,s1d,'',xx,yy,'--')
    xlabel('freq [ppm]');
    legend('rms','max');
    axis([min(ppm) max(ppm) 0 1]);
    set(gca,'xdir','reverse'); grid on
    if ~isempty(fname)
        if isdeployed
            print([fname '_spec.png'],'-dpng','-r300','-painters');
            write_scdicom([fname '_spec.dcm'],gcf,h);
        else
            print([fname '_spec.png'],'-dpng','-r300');
        end
    end
end
if plt>1
    for ll=1:3
        figure(25+ll); clf
        imagesc_row(bb(:,:,:,3+ll),[0 ratsca(ll)],[],true);
        colorbar
        set(gcf,'name',[figstr '- ' bb_str{3+ll}]);
        if ~isempty(fname)
            if isdeployed
                print([fname '_' fff{ll} '.png'],'-dpng','-r300','-painters');
            else
                print([fname '_' fff{ll} '.png'],'-dpng','-r300');
            end
        end
    end
end


%% saving data
if ~isempty(fname)
    spec = single(spec);
    save([fname '.mat'],'-mat','-v7.3',...
        'bb','bb_str','mask','spec','hz','s1d','ppm','par','h','fid',...
        'zf','freqs','nint','rat','bb_scale','mean_frac','thresh_fac');
end


%% exporting as dicom
if ~isempty(fname) && isdeployed
    fprintf('Storing dicom images\n');
    if nz>1
        inp.MRAcquisitionType = '3D';
        if isfield(par,'niso')
            if any(abs(diff(par.niso))>1d-5)
                spacing = round(h.rdb_hdr.fov/nz/par.niso(3)*10)/10;
                inp.SpacingBetweenSlices = spacing;
            end
        else
            warning('no field par.niso');
        end
    end
    for lser=1:(nfreq+nrat)
        if lser<=nfreq
            desc = bb_str{lser};
        else
            desc = sprintf('%s/%s ratio',bb_str{rat(lser-nfreq,1)},...
                bb_str{rat(lser-nfreq,2)});
        end
        inp.SeriesDescription = desc;
        inp.SeriesNumber = 100*h.series.se_no+lser;
        write_dicom([fname '_' num2str(lser)],bb(:,:,:,lser),h,inp);
    end
end

end      % recon_mrsi_xenon.m
