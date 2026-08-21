function [spec,hz,ddd,time,outp] = recon_nutation(d,h,zf,lb,plt,verb)
%RECON_nutation  Reconstruction of nutation acquisition 
%
%[spec,hz,ddd,time,outp] = recon_nutation(d,h,zf,lb,plt,verb)
%
%         d  Raw p-file data [nf,ns,1,1,1,nc]    (or pfile fname)
%            nf=#frames; ns=#spectral pts; nc=#coils
%         h  p-file header                       (or empty)
%        zf  zero-filling                                (default=[])
%        lb  line-broading (exponential Gaussian) [Hz]   (default=[])
%       plt  plotting                                    (default=true)
%      verb  verbose                                     (default=true)
%
%      spec  Spectrum        [ns,nfa]
%        hz  Frequency axis  [Hz]
%       ddd  FID: coil-combined and averaged
%      time  Time axis       [s]
%      outp  misc output structure
%
%  2/2026 Matt Willmering
if nargin<1, help(mfilename); return; end

%% input/default parameters
if ~exist('zf','var'),       zf = []; end
if ~exist('lb','var'),       lb = []; end
if ~isempty(lb), if length(lb)==1, lb(2) = 0; end; end
if ~exist('plt','var'),      plt = []; end
if isempty(plt),             plt = true; end
if ~exist('verb','var'),     verb = []; end
if isempty(verb),            verb = true; end

%% reading in data, if pfile name is given
if ~isnumeric(d) || isempty(d)
    if ~exist(d,'file') && ~isempty(d)
        warning('strange input for d (=''%s''); file not existing?',d); 
    end
    [d,h] = read_p(d,false);
end

%% miscelaneous parameters
csi_dims = h.rdb_hdr.csi_dims;         % #CSI dims (must be 0)
nf   = h.rdb_hdr.nframes;              % acquired frames (size(d,1)
ns   = h.rdb_hdr.frame_size;           % acq pts (user1)
bw   = h.rdb_hdr.spectral_width;       % BW [Hz] (user0)
nc = size(d,6);                        % #coils
time = (0:ns-1)/bw;                    % time axis
if isempty(zf), zf = ns; end           % zero filling; default=off

signal_bound = -2000;                   % cut-off freq for SNR calc [Hz]

%% print info
if verb
    fprintf('TE=%g [ms]; TR=%g [ms]; nacqs(cv4)=%g; BW(cv0)=%g [Hz]; ns(cv1)=%g\n',...
        h.image.te*1d-3,h.image.tr*1d-3,h.image.user4,h.rdb_hdr.user0,ns);
    fprintf('nf=%g, nc=%g',nf,nc);
    fprintf('size(d) = %s\n',num2str(size(d)));
end

%% data checks
nrcvr = (h.rdb_hdr.dab(2)-h.rdb_hdr.dab(1))+1;   % #receiver channels
if nf~=size(d,1),  warning('nf(=%g) ~= size(d,1)(=%g)',nf,size(d,1)); end
if ns~=size(d,2),  warning('ns(=%g) ~= size(d,2)(=%g)',ns,size(d,2)); end
if nc~=nrcvr,      warning('nc(=%g) ~= nrcvr(=%g)',nc,nrcvr); end
if csi_dims~=0,    warning('csi_dim(=%g)~=0: use recon_csi.m for CSI data',csi_dim); end
if (size(d,3)~=1), warning('size(d,3)(=%g)~=1',size(d,3)); end 
if (size(d,4)~=1), warning('size(d,4)(=%g)~=1',size(d,4)); end 
if (size(d,5)~=1), warning('size(d,5)(=%g)~=1',size(d,5)); end 

% unchop data (in case opnex=1)
if ischop(h) && (nf>1), d(2:2:nf,:) = -d(2:2:nf,:); end

%% coil combination
ind1 = (1:nf);
dd = mrs_coil_combine(d,ind1,[],[],verb);      % signal-weighted
% dd = mrs_coil_combine(d,ind1,[],100,verb);     % SNR-weighted
% dd = csi_coil_combine(d,ind1);                   % signal-weighted svds
% warning('check coil combination');

%% line broadening
ddd = dd;
if ~isempty(lb)
    fprintf('Apodising: exponential lb=%g [Hz]; Gaussian lb=%g [Hz]\n',lb(1),lb(2));
    apo = lb_fun(ns,bw,lb);
    ddd = dd.*repmat(apo,[size(ddd,1) 1]);
end

%% spectral reconstruction
hz  = -(-zf/2:zf/2-1)*bw/zf;             % frequency axis
spec = ifftshift(fft(ddd,zf,2),2);

[~,read_max] = max(abs(mean(ddd,1)));
[~,acq_max] = max(abs(ddd(:,read_max))); %read_max to avoid switching

%% SNR
[snr,signal,noise,f_max,ind_noise] = mrs_snr(ddd(acq_max,:),spec(acq_max,:),bw,[],[],false);
[~,spec_max] = min(abs(hz-f_max));

%% line-width
[lw,out] = fwhm(spec,bw,[],false,false);

%% Nutation
FAs = (10*(1:h.rdb_hdr.user4))';%assuming always steps of 10
Intensity = abs(ddd(:,read_max))./max(abs(ddd(:)));

%obtain fit to variables
fit_options = fitoptions('Method','NonlinearLeastSquares',...
           'Lower',[0,0],...
           'Upper',[Inf,Inf],...
           'StartPoint',[1,1]);
fit_type = fittype('abs(a*sind(b*x))','options',fit_options);
[NutationFit,gof] = fit(FAs,Intensity,fit_type);

%create fit to plot
FAsFit = (0:1:max(FAs))';
NutationFitPlot = abs(NutationFit.a*sind(NutationFit.b*FAsFit));

%% print info
if verb
    vol = h.rdb_hdr.roilenx*h.rdb_hdr.roileny*h.rdb_hdr.roilenz*1d-3;
    acq_dur = h.image.tr*1d-6*h.image.user4;
    fprintf('Voxel size = (%g, %g, %g) [mm^3] = %g [ml]\n',...
        h.rdb_hdr.roilenx,h.rdb_hdr.roileny,h.rdb_hdr.roilenz,vol);
    fprintf('Acqusition time = %g [s]\n',acq_dur);
    fprintf('SNR = %g/%g = %.3g\n',signal,noise,snr);
    fprintf('SNR/acq_dur/vol = %.3g [ml*min]^-1\n',snr/vol/acq_dur*60);
    fprintf('linewidth (FWHM) = %.3g [Hz]',lw(1));
    fprintf('\n');
    fprintf('**************************************************\n');
    outp.label = ...
        sprintf('P-file\tdate\tTE\tTR\tcv4\tBW\tns\tnc\tdur\tvol\tSNR\tunitSNR\tLW\n');
    outp.str = ...
        sprintf('P%.5d.7\t%s\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%.3g\t%.3g\t%.3g\n',...
        h.image.rawrunnum,h.rdb_hdr.scan_date,...
        h.image.te*1d-3,h.image.tr*1d-3,h.image.user4,h.rdb_hdr.user0,ns,nc,...
        acq_dur,vol,snr,snr/vol/acq_dur*60,lw(1));
    fprintf('%s',outp.label);
    fprintf('%s',outp.str);    
end

%% saving result as mat file
[fpath,dataname,ext] = fileparts(h.fname);
fname = fullfile(fpath,h.series.se_desc);
if ~isempty(fname)
    fprintf('Saving results to ''%s.mat''\n',fname);
    save([fname '.mat'],'-mat','-v7.3',...
        'hz','time','FAs','NutationFit','snr','spec','ddd','h','lb',...
        'lw','fname');
end


%% plotting
if plt
    ss = abs(spec)./signal;
    figure('Units','inches','Position',[1,1,5,5]);
    hold on
    plot(hz,ss(acq_max,:),'k',...
        hz(ind_noise),ss(acq_max,ind_noise),'cx',...
        f_max(1),abs(ss(acq_max,spec_max)),'cs',...
        out.freqs(acq_max,:),abs(ss(acq_max,out.ind_max(acq_max,:))),'c*',...
        [out.fp(acq_max,:) out.fm(acq_max,:)],...
        repmat(abs(ss(acq_max,out.ind_max(acq_max,:)))/2,[2 1]),'c*-',...
        'MarkerSize',8,'LineWidth',1,'MarkerFaceColor','r');
    set(gca,'XDir','reverse');
    xlim([min(hz) max(hz)]);
    xlabel('Frequency (Hz)');
    %info box
    dim = [0.15 0.7 0.75 0.2];
    str = ['SNR = ',num2str(snr,4),newline,...
        'Linewidth = ',num2str(lw(1)),'Hz',...
        ];
    annotation('textbox',dim,'String',str,'FitBoxToText','on');
    figstr = sprintf('P%05d Exam%d Series%d',...
        h.image.rawrunnum,h.exam.ex_no,h.series.se_no);
    set(gcf,'name',figstr);
    hold off

    if ~isempty(fname)
        print([fname '_spec.png'],'-dpng','-r600');
        saveas(gcf,[fname '_spec.fig'],'fig');
        inp.InstanceNumber = h.rdb_hdr.user4;
        write_scdicom([fname '_spec.dcm'],gcf,h,inp);
    end

    figure('Units','inches','Position',[1,1,5,5]);
    hold on
    plot(FAs,abs(ddd(:,read_max))./max(abs(ddd(:))),'ko-','LineWidth',1,'MarkerSize',8,MarkerFaceColor='c');
    plot(FAsFit,NutationFitPlot,'Color',[0,1,1,0.5],'LineWidth',4)
    xlabel('Flip Angle (Degrees)')
    %info box
    dim = [0.15 0.7 0.75 0.2];
    str = ['True 90deg = ',num2str(90/NutationFit.b,4),'deg',newline,...
        'R^2 = ',num2str(gof.rsquare,4),newline,...
        'S_{max} = ',num2str(NutationFit.a,4),...
        ];
    annotation('textbox',dim,'String',str,'FitBoxToText','on');
    figstr = sprintf('P%05d Exam%d Series%d',...
        h.image.rawrunnum,h.exam.ex_no,h.series.se_no);
    set(gcf,'name',figstr);
    hold off

    if ~isempty(fname)
        print([fname '_nutation.png'],'-dpng','-r600');
        saveas(gcf,[fname '_nutation.fig'],'fig');
        inp.InstanceNumber = h.rdb_hdr.user4+1;
        write_scdicom([fname '_nutation.dcm'],gcf,h,inp);
    end

end

end      % recon_nutation.m
