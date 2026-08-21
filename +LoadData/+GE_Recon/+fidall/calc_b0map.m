function df0 = calc_b0map(bb,TE,mask,do_fiex,phasewrap,fname,h,plt)
%CALC_B0MAP Extract B0 map from 2&3 echo GRE images
%    df0 = calc_b0map(bb,TE,mask,do_fiex,phasewrap,fname,h,plt)
%     bb   Complex-valued images acquired at different TEs
%          size: 3D=[x,y,z,te,coils]; 2D=[x,y,te,1,slice,coils]
%     TE   Echo times                                           [s]
%   mask   Masking parameters or logical                             (true)
%          [#noise,threshold,smoothing]: if true -> mpar = [0.9 0.2 1]
%do_fiex   Apply polynomial smoothing                               (false)
%phasewrap Position of phase wrap: -pi:pi                       [rad]  ([])
%  fname   Save mat+png+dcm to fname                                   ('')
%      h   Header structure (required for dicom_export)                ([])
%    plt   Plotting                                                  (true)
%    df0   B0 map                                               [Hz]
%
% 11/2024 Rolf Schulte
if (nargin<2), help(mfilename); return; end


%% input parameters
if ~exist('mask','var'),    mask = []; end
if ~exist('do_fiex','var'), do_fiex = []; end
if isempty(do_fiex),        do_fiex = false; end
if ~exist('phasewrap','var'), phasewrap = []; end
if ~exist('plt','var'),     plt = []; end
if isempty(plt),            plt = true; end
if ~exist('fname','var'),   fname = ''; end
if ~exist('h','var'),       h = []; end
niso = [];


%% load bb if mat file
if ischar(bb)
    fnmat = bb;
    if ~isempty(regexpi(fnmat,'\.7$',  'once')), fnmat = fnmat(1:end-2); end
    if ~isempty(regexpi(fnmat,'\.h5$', 'once')), fnmat = fnmat(1:end-3); end
    if ~isempty(regexpi(fnmat,'\.mat$','once')), fnmat = fnmat(1:end-4); end
    if exist([fnmat '.mat'],'file')
        xx = load([fnmat '.mat']);
        bb = xx.bb;
        if isempty(h), h = xx.h; end        
        if isfield(xx,'niso'), niso = xx.niso; end
        clear xx
    else
        error('fnmat (=''%s.mat'') not found',fnmat);
    end
    if ~ischar(fname) && ~isempty(fname)
        if fname
            fname = fnmat;
        else
            fname = '';
        end
    end
end


%% inverting sign of input images
fprintf('inverting sign: bb = conj(bb);\n');
bb = conj(bb);


%% parameters
if size(bb,3)<=3
    do_2D = true;
    fprintf('Assuming 2D data\n');
    if size(bb,4)~=1, error('dim4: #metabolites(=%d)~=1',size(bb,4)); end
    if size(bb,6)>1
        pind = [1 2 5 3 6 4];
    else
        pind = [1 2 5 3 4];
    end
    bb = permute(bb,pind);
    if (size(mask,1)==size(bb,1)) && (size(mask,2)==size(bb,2))
        mask = permute(mask,pind);
    end
else
    do_2D = false;
end
nn = size(bb);
dim = length(nn);
if dim<4, error('dim(=%d)<4',dim); end
if dim>5, error('dim(=%d)>5',dim); end
if isreal(bb), error('isreal(bb)'); end
TE = TE(:).';
nTE = size(bb,4);
% if length(TE)==(nTE-1), TE = [0 TE]; end
if length(TE)~=nTE, error('length(TE)(=%d)~=nTE(=%d)',length(TE),nTE); end
if max(TE)>0.1, warning('max(TE)(=%g)>0.1[s]',max(TE)); end


%% create/check fname
if ~isempty(fname)
    if ~isempty(regexpi(fname,'\.7$',  'once')), fname = fname(1:end-2); end
    if ~isempty(regexpi(fname,'\.h5$', 'once')), fname = fname(1:end-3); end
    if ~isempty(regexpi(fname,'\.mat$','once')), fname = fname(1:end-4); end
end


%% calculate phase difference
switch nTE
    case 2
        phi = angle(bb(:,:,:,2).*conj(bb(:,:,:,1)));
        dte = TE(2)-TE(1);
    case 3    % mrm72_80_2014_robinson.pdf
        phi = angle(bb(:,:,:,3).*conj(bb(:,:,:,2).^2).*bb(:,:,:,1));
        dte = TE(3)-2*TE(2)+TE(1);
        
        if false
            % rewind phases
            % pha = angle(bsxfun(@times,bb(:,:,:,2:end),conj(bb(:,:,:,1))));
            pha = angle(bb(:,:,:,2:end).*conj(bb(:,:,:,1:end-1)));
            figure(10); imagesc_ind3d(pha);colorbar
            % pha = bsxfun(@rdivide,pha,reshape(diffte,[1 1 1 nTE-1]));
            % df0 = phi/(2*pi*dte);
            dpha = pha(:,:,:,2) - (TE(3)-TE(2))/(TE(2)-TE(1))*pha(:,:,:,1);
            dpha = angle(exp(1i*dpha));
            figure(11); imagesc_ind3d(dpha);colorbar
            phi = cat(4,pha(:,:,:,1),pha(:,:,:,1) + dpha,phi/3);
            dte = TE(2)-TE(1);
        end

    otherwise
        error('nTE(=%d)~=2 or 3',nTE);
end
fprintf('dte=%g[ms] -> max(df0)=+-%g[Hz]\n',dte*1d3,abs(1/dte/2));


%% weighted coil combination
if dim>4
    fprintf('Coil combination of phase images\n');
    ww = abs(bb(:,:,:,1,:));
    ww = bsxfun(@rdivide,ww,sum(ww,5));
    phi = sum(phi.*ww,5);
end


%% move phase wrap
if ~isempty(phasewrap)
    if abs(phasewrap)>pi, warning('|phasewrap(=%g)|>pi',phasewrap); end
    if phasewrap>0
        phi(phi>phasewrap) = phi(phi>phasewrap) - 2*pi;
    else
        phi(phi<phasewrap) = phi(phi<phasewrap) + 2*pi;
    end
end


%% calculate B0 map
fprintf('Calculating df0\n');
df0 = phi/(2*pi*dte);
bbabs = sqrt(mean(mean(bb.*conj(bb),5),4));


%% masking
do_mask = false;
if isempty(mask), mask = false; end
if all(size(mask,1)==size(df0))
    do_mask = true;
else
    if mask(1)
        do_mask = true;
        if length(mask)<2, mask = [mask 2]; end
        if length(mask)<3, mask = [mask 0.3]; end
        method = mask(2);
        par = mask(3);
        mask = calc_mask(bbabs,method,par);
    end
end
if do_mask, df0 = bsxfun(@times,df0,mask); end


%% polynomial smoothing
if do_fiex
    fprintf('Polynomial smoothing\n');
    fmax = max(df0(:));
    fmin = min(df0(:));
    if do_2D
        for l3=1:size(df0,3)
            df0(:,:,l3) = fiex3d(df0(:,:,l3));
        end
    else
        df0 = fiex3d(df0);
    end


    % crop data
    if any(df0(:)>fmax)
        fprintf('cropping %d voxel to fmax (=%g)\n',sum(df0(:)>fmax),fmax);
        df0(df0>fmax) = fmax;
    end
    if any(df0(:)<fmin)
        fprintf('cropping %d voxel to fmin (=%g)\n',sum(df0(:)<fmin),fmin);
        df0(df0<fmin) = fmin;
    end
end


%% final checks
if any(isnan(df0(:)))
    warning('df0 contains NaN: setting to zero'); 
    df0(isnan(df0)) = 0;
end
if any(isinf(df0(:)))
    warning('df0 contains inf: setting to zero'); 
    df0(isinf(df0)) = 0;
end


%% saving results
if ~isempty(fname)
   save([fname '_df0.mat'],'-mat','-v7.3',...
        'df0','h','nn','dim','TE','do_fiex','phasewrap','niso');
end


%% plotting
if (plt>0) || ~isempty(fname)
    fprintf('Plotting results\n');
    if (plt>0)
        fid = figure(100); clf;
    else
        fid = figure('visible','off');
    end
    set(fid,'DefaultAxesFontSize',14);
    if do_2D
        imagesc_row(squeeze(df0));
    else
        imagesc_ind3d(df0);
    end
    colormap jet
    if ~isempty(h)
        figstr = sprintf('P%05d Exam%d Series%d - B0map',...
            h.image.rawrunnum,h.exam.ex_no,h.series.se_no);
    else
        figstr = 'B0map';
    end
    set(fid,'name',figstr);
    if ~isempty(phasewrap)
        caxis(([-0.5 0.5]-phasewrap/2/pi)/dte);
    else
        caxis([-0.5 0.5]/dte);
    end
    colorbar
    drawnow    
    if ~isempty(fname)
        print([fname '_df0.png'],'-dpng','-r300');
    end
    if ~(plt>0), close(fid); end    
end


%% exporting images as dicom into scanner database
export_dcm = isdeployed;
if ~isempty(fname) && ~isempty(h) && export_dcm
    if isempty(niso), niso = [1,1,1]; end
    fprintf('Exporting images to dicom\n');
    inp.SeriesDescription  = sprintf('MNSRP;B0map;%dD%g;%s',...
        dim,nn(1),h.series.se_desc);
    % inp.MRAcquisitionType = sprintf('%dD',dim);
    if do_2D
        dim_str = '2D';
    else
        dim_str = '3D';
    end
    inp.MRAcquisitionType = dim_str;
    if any(abs(diff(niso))>1d-5)
        spacing = round(h.rdb_hdr.fov/nn(3)/niso(3)*10)/10;
        inp.SpacingBetweenSlices = spacing;
    end
    inp.SeriesNumber = h.series.se_no+100;
    % matlab-based dicomwrite.m without template
    % dicom header generated from p-file header
    write_dicom([fname '_df0'],df0,h,inp,[],[],export_dcm,3);
end


%% reshaping
if do_2D
    df0 = ipermute(df0,pind);
end

end      % calc_b0map.m
