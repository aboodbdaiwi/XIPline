function [rf,g,gg]= design_spsp2d(what,sv,do_blosi,shift4rfstat,par_add)
%DESIGN_SPSP2D  SPSP pulse design by direct 2D matrix inversion
% [rf,g,gg] = design_spsp2d(what,sv,do_blosi,shft4rfstat)
%      what   Parameter selection (1-9) or:
%             inp structure
%             thick     slice thickness                  [m]
%             nspec     #spectral pts
%             dt_spec   duration of flyback sub-block    [s]
%             ftw       fractional transition width
%             ff(2)     spectral bw                      [Hz]
%             mtx_fac   mtx = mtx_fac.*[nspat nspec]
%             overfov   factor beyond Nyquist                    ([1 1.25])
%             samp_fac  spatial downsampling factor              (4)
%             ax_fac    axis factor for sinh function            ([2 0.5])
%             ref       Isocentre of RF pulse                    (0)
%                       0=symmetric; -1=self-refocused; 1=max phase
%             flyback   Flyback (true) or bi-direcional (false) Gz
%             halfshift shift excitation to first n/2 position   (false)
%             gmax      max grad                         [T/m]   (40d-3)
%             smax      max slew rate                    [T/m/s] (150)
%             gmax_fac factor to derate flyback                  (1)
%        sv   Save waveforms + figures
%  do_blosi
%shift4rfstat Centre half-shifted pulse for correct (old) rfstat calc
%   par_add   Parameter structure with selected fields for overwriting inp
%
% 1/2021  Rolf Schulte
if nargin<1, help(mfilename); return; end


%% default parameters
if ~exist('sv','var'),       sv = []; end
if isempty(sv),              sv = false; end
if ~exist('do_blosi','var'), do_blosi = []; end
if isempty(do_blosi),        do_blosi = true; end
if ~exist('shift4rfstat','var'), shift4rfstat = []; end
if isempty(shift4rfstat),    shift4rfstat = false; end
if ~exist('par_add','var'),  par_add = []; end
if isempty(par_add),         par_add = struct; end
hwdt = 4d-6;                 % sequencer sampling time [s]
b1max_ref = 0.128128;        % fidall B1max reference [G]


%% design parameters
if isstruct(what)
    inp = what;
    clear what;
else
    if isnumeric(what)
        inp = sub_what2par(what);     % stored pulses
    else
        error('strange input: what');
    end
end

%% default values
def = {{'gmax',40d-3},{'smax',150},{'gmax_fac',1},{'halfshift',false},...
    {'nucleus','13C'},{'ref',0},{'ax_fac',[2 0.5]},{'overfov',[1 1.25]},...
    {'samp_fac',4},{'reg',0},{'rampsamp',0},{'do_real',false}};
for l=1:length(def)
    if ~isfield(inp,def{l}{1})
        inp.(def{l}{1}) = def{l}{2};
    end
end


%% overwrite with par_add
fldn = fieldnames(par_add);
for l=1:length(fldn)
    inp.(fldn{l}) = par_add.(fldn{l});
end


%% parameters
thick = inp.thick;
nspec = inp.nspec;
dt_spec = inp.dt_spec;
ftw = inp.ftw;
mtx_fac = inp.mtx_fac;
overfov = inp.overfov;
ax_fac = inp.ax_fac;
ref = inp.ref;
flyback = inp.flyback;
halfshift = inp.halfshift;
gmax = inp.gmax;
smax = inp.smax;
gmax_fac = inp.gmax_fac;
nucleus = inp.nucleus;
samp_fac = inp.samp_fac;
reg = inp.reg;
rampsamp = inp.rampsamp;
do_real = inp.do_real;


%% checks
% if halfshift, ref = -1; end
if ~flyback && gmax_fac~=1, error('~flyback && gmax_fac~=1'); end
% if flyback && halfshift, warning('cannot design shifted flyback'); end
if shift4rfstat && ~halfshift
    warning('options shift4rfstat does not make sense for non-halfshift pulse'); 
end


%% physical parameters
if abs(dt_spec/hwdt-round(dt_spec/hwdt)) > 1d-12
    dt_spec = round(dt_spec/hwdt)*hwdt;
    warning('dt_spec not multiple of hwdt; rounding dt_spec to %g[ms]',...
        dt_spec*1d3);
end
ff(1) = gyrogamma(nucleus)/2/pi*thick*gmax; % bw of spatial dim
ff(2) = inp.ff(2);


%% design gradient trajectory
[g,ind_rf,nwf,nplat,nramp] = sub_design_grad(nspec,dt_spec,hwdt,flyback,gmax,smax,...
    gmax_fac,rampsamp);


%% filename for saving
t_total = length(g)*hwdt;
if halfshift
    type_str = 'hshi';
else
    type_str = 'cent';
end
% if flyback, type_str = 'flyb'; end
if flyback, type_str = [type_str '_fb']; end

fname = sprintf('spsp2D_%s_%s_nspec%g_dt%g_ff%g_%g_ftw%g_%g_ref%g_thick%g_dur%g',...
    type_str,nucleus,nspec,round(dt_spec*1d5)/100,...
    round(ff(1)),round(ff(2)),round(ftw(1)*100)/100,round(ftw(2)*100)/100,...
    round(ref*10)/10,round(thick*1d4)/10,...
    round(t_total*1d4)/10);
fname(strfind(fname,'.'))='p';
fname(strfind(fname,'-'))='m';
if rampsamp>0, fname = [fname '_rs']; end
fprintf('fname=%s\n',fname);
if sv
    fileid = fopen([fname 'design_spsp2d_output.txt'],'w+t');
else
    fileid = 1;
end
fprintf(fileid,'fname=%s\n',fname);


%% print grad info
fprintf(fileid,'Prescribed gmax=%g [mT/m]  smax=%g [T/m/s]\n',...
    gmax*1d3,smax);
smax_act = max(abs(diff(g)/hwdt));
gmax_act = max(abs(g));
fprintf(fileid,'Actual     gmax=%g [mT/m]  smax=%g [T/m/s]\n',...
    gmax_act*1d3,smax_act);
if (gmax_act>gmax), warning('gmax_act(=%g)>gmax(=%g)',gmax_act,gmax); end
if (smax_act>smax), warning('smax_act(=%g)>smax(=%g)',smax_act,smax); end
fprintf(fileid,'Total duration = %g [us]\n',t_total*1d6);


%% actual RF pulse design
% convert grad into k-space trajectory
k = cumsum(g(end:-1:1));               % spatial k-space        [s]
k = k(end:-1:1)*hwdt/gmax;
kf = ((nwf-1):-1:0)*hwdt;              % spectral k-space       [s]
% kf = ((nwf/2-1):-1:(-nwf/2))*hwdt;     % spectral k-space       [s]

bw = [1/(samp_fac*hwdt) 1/dt_spec];    % Nyquist BW for fitting [Hz Hz]
if any(ff>bw), error('ff(=%g,%g)>bw(=%g,%g)',ff(1),ff(2),bw(1),bw(2)); end
ind = ind_rf;                          % complete index list for rf pts
tmp = logical((0:nwf-1)/samp_fac-floor((0:nwf-1)/samp_fac) ~= 0);
ind(tmp) = false;                      % index for rf downsampled
ni = sum(ind);                         % # rf pts for fitting
indf = find(ind);                      % index list
nspat = round(sum(ind)/nspec);
mtx = round(mtx_fac.*[nspat nspec]);   % fit size: (#spatial, #spectral)

fprintf(fileid,'nspat=%g, nspec=%g; ',nspat,nspec);
fprintf(fileid,'mtx=(%g,%g); samp_fac=%g\n',mtx(1),mtx(2),samp_fac);
fprintf(fileid,'bw(spat,spec)=(%g,%g) [Hz]\n',bw(1),bw(2));
fprintf(fileid,'ff(spat,spec)=(%g,%g) [Hz]\n',ff(1),ff(2));


%% generate grid
% spectral dim > Nyquist, as dt_spec only average; minimum is ramp time
F = zeros(2,6);
for l=1:2          % 1=spat,2=spec dim: F=normalised freq [-0.5..0.5]
    tmp1 = ff(l)/bw(l)*(1-ftw(l));
    tmp2 = ff(l)/bw(l)*(1+ftw(l));
    F(l,:) = [-1 -tmp2 -tmp1 tmp1 tmp2 1]/2;

    switch 1
        case 1
            tmp = sinh(linspace(-1,1,mtx(l)+1)*ax_fac(l));
            tmp = tmp(1:end-1);
            ax{l} = tmp/max(tmp)/2*overfov(l);
        case 2
            ax{l} = (-mtx(l)/2:mtx(l)/2-1)/mtx(l)*overfov(l);
    end
end
fprintf(fileid,'spatial:  F=(%g,%g,%g,%g,%g,%g) [Hz]\n', ...
    round(F(1,:)*bw(1)));
fprintf(fileid,'spectral: F=(%g,%g,%g,%g,%g,%g) [Hz]\n', ...
    round(F(2,:)*bw(2)));

[gri(:,:,1),gri(:,:,2)] = ndgrid(ax{1},ax{2});   % generate 2D grid
% generate mask and exclude transition bands
if halfshift
    % excitation profile (abs)
    b2d = gri(:,:,1)>F(1,3) & gri(:,:,1)<F(1,4) & ...
        ((gri(:,:,2)>(F(2,3)-0.5) & gri(:,:,2)<(F(2,4))-0.5) | ...
         (gri(:,:,2)>(F(2,3)+0.5) & gri(:,:,2)<(F(2,4))+0.5)); 
    mask = gri(:,:,1)>F(1,5) | gri(:,:,1)<F(1,2) | ...
        gri(:,:,2)>(F(2,5)+0.5) | gri(:,:,2)<(F(2,2)-0.5) | ...
        (gri(:,:,2)>(F(2,5)-0.5) & gri(:,:,2)<(F(2,2)+0.5)) | ...
        b2d; % actual mask
else
    b2d = (gri(:,:,1)>F(1,3) & gri(:,:,1)<F(1,4) & ...
        gri(:,:,2)>F(2,3) & gri(:,:,2)<F(2,4));  % excitation profile (abs)
    mask = gri(:,:,1)>F(1,5) | gri(:,:,1)<F(1,2) | ...
        gri(:,:,2)>F(2,5) | gri(:,:,2)<F(2,2) | b2d; % actual mask
end
ngri = sum(mask(:));                             % # of fit points
griV = reshape(gri(repmat(mask,[1 1 2])),ngri,2); % actual exclusion

% A: excitation matrix
xx = bsxfun(@times,griV,-2*pi*1i*bw);
% A = exp(xx(:,1)*k(ind)+xx(:,2)*kf(ind));
ktmp = kf(ind); ktmp = ktmp-mean(ktmp);
A = exp(xx(:,1)*k(ind)+xx(:,2)*ktmp);


% b: target profile (incl linear spectral phase)
% isodelay = (t_total+(-1+ref)*nspec*dt_spec/2)/t_total;
isodelay = ((ref)*nspec*dt_spec/2)/t_total;
if false % ref==-1    % linear phase does not properly work for halfshift, use -1
    b = b2d(mask);
else
    b = b2d(mask).*exp(-2*pi*nspec*isodelay*1i*griV(:,2));
end
fprintf(fileid,'size(A)=(%gx%g), size(b)=(%gx%g)\n',...
    size(A,1),size(A,2),size(b,1),size(b,2));
fprintf(fileid,'isodelay = %g\n',isodelay);

% W: weighting matrix
switch 2
    case 1
        W = ones(mtx(1),mtx(2)) + 2*b2d + ...
            2*(gri(:,:,1)>F(1,2)*1.2 & gri(:,:,1)<F(1,5)*1.2 & ...
            ((gri(:,:,2)>(0.5-F(2,4))) | ((gri(:,:,2)<(-0.5-F(2,3))))));
        W = W(mask).';
    case 2
        W = [];
end

% figure(1);x=(1:length(k));plot(x,k,'b',x(ind_rf),k(ind_rf),'rx')


%% plot design
figure(11);
subplot(2,1,1);   % spatial
[ln1,ln2] = find(b2d,1);
plot(F(1,:)*bw(1),[0 0 1 1 0 0],'b',...
    ax{1}(mask(:,ln2))*bw(1),b2d(mask(:,ln2),ln2),'c.');
axis([-bw(1)/2 bw(1)/2 -0.1 1.15]); grid on; title('spatial');

subplot(2,1,2);   % spectral: plot twice bw
if halfshift
    plot((F(2,:)+0.5)*bw(2),[0 0 1 1 0 0],'b',...
        (F(2,:)-0.5)*bw(2), [0 0 1 1 0 0],'b',...
        (F(2,:))*bw(2),     [0 0 1 1 0 0],'g:',...
        (F(2,:)+1)*bw(2),   [0 0 1 1 0 0],'g:',...
        (F(2,:)-1)*bw(2),   [0 0 1 1 0 0],'g:',...
        ax{2}(mask(ln1,:))*bw(2),b2d(ln1,mask(ln1,:)),'c.');
    if strcmpi(nucleus,'13c')
        hold on
        [f13c,f_str] = freqs13c;
        off_freq = bw(2)/2*[1 0 1 -1 -1];
        for l=[1 3 4 5]
            yy = (l-1)/8+0.1;
            spsp_freq = f13c(l) - off_freq(l);
            plot(f13c - spsp_freq,ones(5,1)*yy,'rx');
            plot(2*bw(2)*[-1 1],[1 1]*yy,'y');
            for l2=1:5
                text(freqs13c(l2)-spsp_freq,yy,f_str{l2},'fontsize',8)
            end
        end
    end
else
    plot(F(2,:)*bw(2),[0 0 1 1 0 0],'b',...
        (F(2,:)+1)*bw(2),[0 0 1 1 0 0],'b:',...
        (F(2,:)-1)*bw(2),[0 0 1 1 0 0],'b:',...
        ...%(F(2,:)+0.5)*bw(2),[0 0 1 1 0 0]/2,'g:',...
        ...%(F(2,:)-0.5)*bw(2),[0 0 1 1 0 0]/2,'g:',...
        ax{2}(mask(ln1,:))*bw(2),b2d(ln1,mask(ln1,:)),'c.');
    if strcmpi(nucleus,'13c')
        hold on
        for l=1:5
            yy = (l-1)/8+0.1;
            [f13c,f_str] = freqs13c;
            plot(f13c - f13c(l),ones(5,1)*yy,'rx');
            plot(2*bw(2)*[-1 1],[1 1]*yy,'y');
            for l2=1:5
                text(freqs13c(l2)-freqs13c(l),yy,f_str{l2},'fontsize',8)
            end
        end
    end
end
hold off; grid on; axis([-bw(2) bw(2) -0.1 1.15]);
title('spectral');
% drawnow


%% actual fit: A=excitation matrix; b=target profile
% differences in rf_fit tiny, hence use fastest
% linear least quares fit
if isempty(W)
    if reg>0
        pinvA = pinv(A'*A+reg*eye(size(A,2)));
        rf_fit = pinvA*(A'*b);
    else
        rf_fit = (A'*A)\(A'*b);            % 2.1s (example runtime in profiler)
        % rf_fit = inv(A'*A)*(A'*b);       % 2.6s
        % rf_fit = A\b;                    % 13.2s
        % rf_fit = pinv(A)*b;              % 28.5s
    end
else
    % rf_fit = ((A'*W*A)\(A'*W*b));    % weighted lin lsq as with diag(W)
    % tmp = bsxfun(@times,A',W.');
    % rf_fit = ((tmp*A)\(tmp*b));

    rf_fit = ((bsxfun(@times,A',W)*A)\(bsxfun(@times,A',W)*b));
    % rf_fit = lscov(A,b,W);
    % similar in speed as storing --> better for memory
end
% rf_fit = (A'*A+10*eye(ni))\(A'*b);   % regularised lin lsq


%% shift half-shifted pulse to on-resonant excitation
if shift4rfstat
    fm = ones(1,nwf)*bw(2)/2; % frequency modulation
    pm = -cumsum(2*pi*fm*hwdt); % phase modulation
    rf_fit = rf_fit.*exp(1i*pm); % add phase modulation to rf pulse
end


%% checking + setting imag(rf_fit) to zero
sarf = sum(abs(real(rf_fit)));
saif = sum(abs(imag(rf_fit)));
fprintf(fileid,...
    'sum(abs(real(rf_fit)))=%g, sum(abs(imag(rf_fit)))=%g; ratio=%g\n',...
    sarf,saif,sarf/saif);
if sarf>saif*100 || do_real
    fprintf(fileid,'setting imag(rf_fit)=0\n');
    rf_fit = real(rf_fit);
else
    % fprintf(fileid,'Warning: large imag(rf_fit)!?\n');
end


%% interpolate and scale rf to physical full waveform in [T]
rf = interp1(indf(1:ni),rf_fit,(1:nwf))/(gyrogamma(nucleus)*hwdt*samp_fac);
rf(~ind_rf) = 0; rf(isnan(rf)) = 0;


%% calculate actual isodelay
tmp = linspace(1,0,nwf);
if ~isreal(rf), warning('using real for centre of mass for isodelay'); end
rf_cog = real(rf);
if halfshift
    if flyback
        rf_cog = abs(rf_cog);        
    else
        rf_cog = rf_cog.*sign(g);
    end
end
isodelay_act = sum(tmp.*real(rf_cog))/sum(real(rf_cog));   % centre-of-mass 
fprintf(fileid,'actual isodelay=%g; TE+=%.4g[ms]\n',...
    isodelay_act,isodelay_act*t_total*1d3);


%% misc plotting, output, conversion to scanner waveforms
% plotting of fit
figure(12);
subplot(2,2,1); imagesc(abs(b2d),[0 1]); axis off; title('b target');
b_fit = zeros(mtx);
b_fit(mask) = A*rf_fit;                % actual fit response
tmp = max(abs(b_fit(:)));
b_fit = b_fit/tmp;                     % scale fit response
rf_fit = rf_fit/tmp;                   % scale fit coefficients
res = abs(b_fit-b2d);
subplot(2,2,2); imagesc(abs(b_fit),[0 1]); axis off; title('b\_fit');
subplot(2,2,3); imagesc(res,[-0.2 0.2]); axis off; 
title('residual [-0.2 0.2]');
subplot(2,2,4); 
if isreal(rf_fit)
    plot(rf_fit);
else
    plot_complex(rf_fit); 
end
title('rf\_fit'); axis tight; drawnow;


%% calculate area, max B1, and print RF info
area     = sum(real(rf))/length(rf)/max(abs(rf));
abswidth = sum(abs(rf))/length(rf)/max(abs(rf));
nom_fa = 90;
% maximum B1 for proton
if halfshift && ~shift4rfstat
    % use abswidth
    b1max = nom_fa * 1e6 / (360 * t_total*1e6 * abswidth * gyrogamma(1)/2/pi/1e4);
else
    % use area
    b1max = nom_fa * 1e6 / (360 * t_total*1e6 * area * gyrogamma(1)/2/pi/1e4);
end
if b1max < max(abs(rf))/gyrogamma(1)*gyrogamma(nucleus)*1e4
    warning('b1max = %d is unexpected small (< %d)',...
        b1max,max(abs(rf))/gyrogamma(1)*gyrogamma(nucleus)*1e4);
end
fprintf(fileid,'1H B1max=%g[G]=%g[uT];  f1max = %g [Hz]\n',...
    b1max,b1max*100,b1max/1e4*gyrogamma(1)/2/pi);
tmp = b1max*gyrogamma(1)/gyrogamma(nucleus);
if ~strcmpi(nucleus,'1h')
    fprintf(fileid,'%s B1max = %g [G]\n',nucleus,tmp);
end
if b1max>b1max_ref
    warning('B1max(=%g)>%g (fidall reference)',b1max,b1max_ref);
end
fprintf(fileid,'error: RMS = %g; ',sqrt(mean(res(mask).^2)));
fprintf(fileid,'max = %g\n',max(res(mask)));


%% plotting overall waveform
rf = rf*b1max*1d-4/max(abs(rf));
figure(13);
xt = (0:nwf-1)*hwdt*1d3;
subplot(2,1,1);
if isreal(rf)
    plot(xt,rf*1d6);
else
    plot_complex(xt,rf*1d6);
end
hold on
plot([0 xt(end)],[1 1]*b1max_ref*100,'m:');
hold off
tmp = max(abs(rf))*1d6;
axis([0 xt(end) (min(real(rf))*1d6-0.1*tmp) 1.1*tmp]);
xlabel('time [ms]'); ylabel('B1 [uT]');
drawnow
subplot(2,1,2);
plot(xt,g*1d3,'b',xt(ind_rf),g(ind_rf)*1d3,'c.','MarkerSize',2);
axis([0 xt(end) 1.2d3*min(g) 1.2d3*max(g)]);
xlabel('time [ms]'); ylabel('G [mT/m]');
drawnow


%% EPIC RF wave form parameters
% minimum slice thickness [mm] for protons
thick = thick*gyrogamma(nucleus)/gyrogamma('1h');
gg = [g ; ones(size(g))*2*pi/gyrogamma(nucleus)];


%% adding profile from fit to target plot
% kk = cumsum(gg(:,end:-1:1),2);
% kk = gyrogamma(nucleus)/2/pi*kk(:,end:-1:1)*hwdt;
nn = 1000; xx = -(-nn/2:(nn/2-1)).'/nn;
A1 = exp(-2*pi*1i*xx*k*bw(1));
A2 = exp(-2*pi*1i*xx*kf*bw(2)*4);
rfX = rf.'*gyrogamma(nucleus)*hwdt*2/pi;
figure(11);
subplot(2,1,1); hold on;
plot(bw(1)*xx,real(A1*rfX),'r');
hold off;
subplot(2,1,2); hold on;
plot(bw(2)*4*xx,abs(A2*rfX),'r');
hold off;


%% calculating actual excitation profile
if do_blosi
    figure(14); clf;
    % blosi(rf,gg,[ceil(thick*300)*0.01 ceil(ff(2)*0.3)*10],[100 120],[],nucleus);
    mtx_blosi = [80 100];
    if sv, mtx_blosi = mtx_blosi*2; end
    [Mxy,Mz,ax]=blosi(rf,gg,[thick*4 bw(2)*2],mtx_blosi,[],nucleus);
    xlabel('freq [Hz]'); ylabel('z [m]');
    drawnow
end


%% save waveform + figures
if sv
    % plot figures and write into file
    % figure(11); print_all([fname '_fig11']);
    % figure(12); print_all([fname '_fig12']);
    % figure(13); print_all([fname '_fig13']);
    % if do_blosi, figure(14); print_all([fname '_fig14']); end
    figure(11); print('-dpng','-r600',[fname '_fig11']);
    figure(12); print('-dpng','-r600',[fname '_fig12']);
    figure(13); print('-dpng','-r600',[fname '_fig13']);
    if do_blosi, figure(14); print('-dpng','-r600',[fname '_fig14']); end
    
    % put parameters into structure
    par.rftype = 'spsp';
    par.pw = t_total*1d6;
    par.flip = nom_fa;
    par.bwcs = ff(2);
    par.nlobe = nspec;
    par.area = nplat/(nplat+2*nramp);
    par.bwsp = ff(1);
    par.res = length(rf);
    par.minwidth = thick*1d3; % for proton [mm]
    par.pwlobe = nplat*hwdt*1d6;
    par.pwlobea = nramp*hwdt*1d6;
    par.pwlobed = nramp*hwdt*1d6;
    par.isodelay = isodelay_act;
    par.maxB1 = b1max; % max B1 1H [G]
    
    if isreal(rf)
        Create5xxRfExtFile([fname '.rho'], rf/max(abs(rf)), par);
    else
        Create5xxRfExtFile([fname '.rho'], abs(rf)/max(abs(rf)), par);
        Create5xxRfExtFile([fname '.theta'], angle(rf)/pi, par);
    end
    Create5xxRfExtFile([fname '.gz'], g/max(abs(g)), par);
    
    save(fname,'F','Mxy','Mz','W','abswidth','area','ax','ax_fac',...
        'b','b1max','b1max_ref','b2d','b_fit','bw','def','do_blosi',...
        'dt_spec','ff','flyback','fname','ftw','g','gg','gmax',...
        'gmax_act','gmax_fac','gri','griV','halfshift','hwdt','ind',...
        'ind_rf','inp','isodelay','isodelay_act','k','kf','mask',...
        'mtx','mtx_blosi','mtx_fac','nn','nom_fa','nplat','nramp',...
        'nspat','nspec','nucleus','nwf','overfov','par','rampsamp',...
        'ref','reg','res','rf','rf_fit','samp_fac','shift4rfstat',...
        'smax','smax_act','t_total','thick','type_str','what','xt');
    fclose(fileid);
end
% whos

% to simulate with blosi:
% [Mxy,Mz,ax]=blosi(rf,gg,[60d-3 2d3],[120 180],[],13);
% xlabel('freq [Hz]'); ylabel('z [m]');

end % function design_spsp2d


%% sub-functions
function [g,ind_rf,nwf,nplat,nramp] = sub_design_grad(nspec,dt_spec,hwdt,flyback,...
    gmax,smax,gmax_fac,rampsamp)
npts = ceil(nspec*dt_spec/hwdt);       % total # pts excl rewinder
g = zeros(1,npts);                     % init
nsub = floor(dt_spec/hwdt);            % # pts for one grad lobe

if ~flyback   % bi-directional trajectory    
    nramp = ceil(gmax/smax/hwdt)+1;    % # pts for ramp
    nplat = nsub - 2*nramp + 1;        % # pts for plateau
    gsub = [linspace(0,gmax,nramp) , ...
        ones(1,nplat)*gmax , ...
        linspace(gmax,0,nramp)];
    nn = floor(nplat/2-nramp/2);       % # pts of rewinder plateau
    if nn<0, error('rewinder plateau nn(=%g)<0',nn); end
    nwf = npts+nn+2*nramp;             % total # waveform pts 
    ind_rf = false(1,nwf);             % index to pts for rf (plateau)

    % combine sub-lobes for main trajectory
    if rampsamp>0
        nover = floor(nramp/4)*2;
    else
        nover = 0;                     % >0 ->ramp sampling
    end
    for l=1:nspec
        ind_g = (l-1)*nsub + (1:nsub);
        g(1,ind_g) = (isodd(l)*2-1)*gsub(1:end-1);
    
        ind_i = (l-1)*nsub + (1:(nplat+2*nover)) + nramp - nover;
        ind_rf(ind_i) = true;
    end
    % add rewinder
    grefo = [linspace(0,gmax,nramp) ones(1,nn)*gmax linspace(gmax,0,nramp)];
    g = [g -(isodd(nspec)*2-1)*grefo];
else          % flyback trajectory
    gmax2 = gmax/gmax_fac;                  % reduced gmax for RF
    nramp2 = ceil(gmax2/smax/hwdt);         % # pts for RF ramp
    
    gmax1 = gmax;                           % gmax for return
    nramp1 = ceil(gmax1/smax/hwdt);         % # pts for return ramp
    nplat = nsub-2*nramp1-2*nramp2;         % # pts for both plateaus
    if nplat<0, error('nplat<0'); end
    if nramp1*gmax1>(nramp2+nplat)*gmax2
        fprintf('no return nplateau required; shorten return lobe\n');
        nplat1 = 0;
        nplat2 = nplat;
        while true
            if (nramp1-2)*gmax1<(nramp2+nplat2+4)*gmax2
                break;
            end
            nramp1 = nramp1-1;
            nplat2 = nplat2+2;
            gmax1 = nramp1*smax*hwdt;
        end
    else
        fprintf('return nplateau required\n');
        nplat1 = 0;
        nplat2 = nplat;
        while true
            if (nramp1+nplat1)*gmax1>(nramp2+nplat2)*gmax2
                break;
            end
            nplat1 = nplat1+1;
            nplat2 = nplat2-1;
        end
    end
    % fine tuning to remove errors from time rounding
    gmax1 = gmax2*(nramp2+nplat2)/(nramp1+nplat1);

    % actual construction of flyback sublobe
    tmp1 = linspace(0,-gmax1,nramp1+1);
    tmp2 = -ones(1,nplat1)*gmax1;
    tmp3 = linspace(-gmax1,0,nramp1+1);
    tmp4 = linspace(0,gmax2,nramp2+1);
    tmp5 = ones(1,nplat2)*gmax2;
    tmp6 = linspace(gmax2,0,nramp2+1);
    gsub = [tmp1(1:end-1) tmp2 tmp3(1:end-1) ...
        tmp4(1:end-1) tmp5 tmp6(1:end-1)];
    if sum(gsub)>1d-6, warning('sum(gsub)>1d-6'); end
    
    nwf = npts+2*nramp1+nplat1+1;           % total # waveform pts 
    ind_rf = false(1,nwf);                  % index for RF plateau
    
    % combine sub-lobes for main trajectory    
    if rampsamp>0
        nover = floor(nramp2/3)*2;
    else
        nover = 0;                     % >0 ->ramp sampling
    end
    noff = 2*nramp1 + nplat1 + nramp2 - nover;
    
    for l=1:nspec
        ind_g = (l-1)*nsub + (1:nsub);
        g(1,ind_g) = gsub;
    
        % if rampsamp>0
        %     ind_i = (l-1)*nsub + (1:(nplat2+2*nramp2)) + 2*nramp1 + nplat1;
        % else
        %     ind_i = (l-1)*nsub + (2:nplat2) + 2*nramp1 + nplat1 + nramp2;
        % end
        ind_i = (l-1)*nsub + (1:(nplat2+2*nover)) + noff;
        ind_rf(ind_i) = true;
    end
    
    % add rewinder
    grew = [tmp1(1:end-1) tmp2 tmp3];
    g = [g 0.5*grew];
    
    % remove first flyback lobe
    ii = 2*nramp1+nplat1;
    g = [0 g(1,(ii+1):end)];
    ind_rf = ind_rf(1,ii:end);
    nwf = nwf-ii+1;
    
    % use average values for simulator
    nplat = round((nplat1+nplat2)/2);
    nramp = round((nramp1+nramp2)/2);
end

if isodd(size(g,2))
    g = [0 g];
    ind_rf = [false ind_rf];
    nwf = nwf+1;
end
if any(size(g)~=size(ind_rf)), error('size(g)~=size(ind_rf)'); end
if rampsamp==2, ind_rf = true(size(ind_rf)); end

end      % sub-function sub_design_grad


%% example pulses
function par = sub_what2par(what)

switch what
    case 1    % alias outside; CS centered
        par.thick = 8d-3;             % slice thickness               [m]
        par.nspec = 21;               % # of spectral pts
        par.dt_spec = 1.120d-3;       % duration of sub-block         [s]
        par.ftw = [0.5 0.7];
        par.ff(2) = 80;               % spectral bw 
        par.mtx_fac = [4 10];         % mtx = mtx_fac.*[nspat nspec]
        par.overfov = [1 1.15];       % factor beyond Nyquist
        par.ax_fac = [2 2.5];
        par.ref = -1;
        par.flyback = false;
    case 2    % alias outside; CS shifted by half
        par.thick = 5d-3;             % slice thickness               [m]
        par.nspec = 21;               % # of spectral pts
        par.dt_spec = 1.216d-3;       % duration of sub-block         [s]
        par.ftw = [0.5 0.75];
        par.ff(2) = 80;               % spectral bw
        par.mtx_fac = [4 8];          % mtx = mtx_fac.*[nspat nspec]
        par.halfshift = true;
        par.flyback = false;
    case 3    % alias outside; CS shifted by half
        par.thick = 8d-3;
        par.nspec = 21;
        par.dt_spec = 1.088d-3;
        par.ftw = [0.5 0.7];
        par.ff(2) = 80;
        par.mtx_fac = [4 8]*2;
        par.halfshift = true; 
        par.flyback = false;
    case 4    % fly-back design
        par.thick = 14d-3;
        par.nspec = 18;
        par.dt_spec = 1.216d-3;
        par.ftw = [0.8 0.7];
        par.ff(2) = 85;
        par.mtx_fac = [30 20];
        par.overfov = [2 1.25];
        par.ax_fac = [1 0.5];
        par.ref = -1;
        par.flyback = true;
        par.gmax_fac = 2; 
    case 5    % alias outside; CS shifted by half; large BW (no Lac+PyrH)
        par.thick = 8d-3;
        par.nspec = 15;
        par.dt_spec = 1.120d-3;
        par.ftw = [0.7 0.8];
        par.ff(2) = 85;
        par.mtx_fac = [10 30];
        par.halfshift = true;
        par.flyback = false;
    case 6    % 1H flyback pulse for 2D MMRSI: excite metabolites, not H2O+fat
        % stop-freq: H2O=4.7-0.1ppm; lipid=1.3+0.1ppm
        % pass-freq: 4.1 & 1.9ppm
        % set spsp_fatx1=-218
        par.nucleus = '1H';
        par.thick = 7d-3;
        par.nspec = 22;
        par.dt_spec = 1.5d-3;
        par.ftw = [0.35 0.17];
        par.ff(2) = 325;
        par.mtx_fac = [5 10];
        par.overfov = [2 1.37];
        par.ref = -0.85;
        par.flyback = true;
        par.gmax_fac = 2;
        par.samp_fac = 4;
        par.gmax = 30d-3;
        par.smax = 120;
        par.reg = 1d3;
        par.rampsamp = true;
    case 7    % 1H flyback pulse for 3D MRSI: excite metabolites, not H2O+fat
        % stop-freq: H2O=4.7-0.1ppm; lipid=1.3+0.1ppm
        % pass-freq: 4.1 & 1.9ppm
        par.nucleus = '1H';
        par.thick = 40d-3;
        par.nspec = 25;
        par.dt_spec = 1.3d-3;
        par.ftw = [0.1 0.17];
        par.ff(2) = 325;
        par.mtx_fac = [1.5 3];
        par.overfov = [2 1.42];
        par.ref = -0.83;
        par.flyback = true;
        par.gmax_fac = 6;
        par.samp_fac = 2;
        par.gmax = 20d-3;
        par.smax = 120;
        par.reg = 1d3;
        par.rampsamp = true;
    case 8    % 1H pulse for 2D MRI: excite H2O, no fat
        % lipid-freqs relative to water: 
        % (4.7-1.3)*128=435Hz
        % (4.7-0.9)*128=486Hz
        par.nucleus = '1H';
        par.thick = 3.2d-3;
        par.nspec = 10;
        par.dt_spec = 1.55d-3;
        par.ftw = [0.2 0.6];
        par.ff(2) = 170;
        par.mtx_fac = [4 10];
        par.overfov = [1.5 1.55];
        par.ref = -0.55;
        par.flyback = false;
        par.samp_fac = 4;
        par.gmax = 33d-3;
        par.smax = 100;
        par.reg = 1;%5d2;
        par.do_real = true;
        par.rampsamp = 1;
    case 9    % 1H flyback pulse for 2D MRI: excite H2O, no fat
        % lipid-freqs relative to water: 
        % (4.7-1.3)*128=435Hz
        % (4.7-0.9)*128=486Hz
        par.nucleus = '1H';
        par.thick = 3.2d-3;
        par.nspec = 10;
        par.dt_spec = 1.5d-3;
        par.ftw = [0.3 0.4];
        par.ff(2) = 260;
        par.mtx_fac = [4 10];
        % par.overfov = [1.5 1.3];
        par.ref = -0.5;
        par.halfshift = true;
        par.flyback = true;
        par.gmax_fac = 2; 
        par.samp_fac = 4;
        par.gmax = 33d-3;
        par.smax = 120;
        par.reg = 5d2;
        par.do_real = true;
        par.rampsamp = 1;% 2;
    otherwise
        error('what not found');
end
end      % sub-function sub_what2par



