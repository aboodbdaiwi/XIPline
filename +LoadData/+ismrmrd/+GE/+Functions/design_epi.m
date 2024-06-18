function [k,dcf,t,grad,out,ga_burst] = design_epi(fov,mtx,nexc,gmax,smax,...
    rampsamp,ref,do_burst,undsamp,fname,nucleus,samp,bwmax,gdt,verb)
%DESIGN_EPI  Design Cartesian readout: prewinder and freq encode
% [k,dcf,t,grad,out] = design_epi(fov,mtx,nexc,gmax,smax,...
%       rampsamp,ref,do_burst,undsamp,fname,nucleus,samp,bwmax,gdt,verb)
%                                                         unit     default
%       fov  Field of view                                [mm]
%       mtx  Matrix size
%      nexc  Number of excitations
%      gmax  Maximum gradient amplitude                   [mT/m]
%      smax  Maximum slew rate                            [T/M/s]
%     *****  fov,mtx,gmax,smax: 1 value  => isotropic
%                               2 values => 1=freq+2=phase
%rest optional (use [] for default):
%  rampsamp  Ramp sampling (requires regridding)          logical  true
%       ref  EPI reference scan                                    0
%            null phase encode of frame ref(1)
%            total #frames            = ref(2)
%            null complete EPI of frame ref(3)
%            null phase encode of frame ref(4) and invert epsi
%  do_burst  Segment for burst                                     false
%   undsamp  Design undersampled trajectory: NOT YET IMPLEMETED    0 
%            0=off (fully); 1=parallel imaging (ARC) 
%            2=partial k-space; 3=1+2=all
%     fname  filename (if given->write_ak_wav) (true=generate)     false
%   nucleus  Nucleus                                               '1H'
%      samp  Sampling rate: 1=fully, >1=Over, <1=under sampled     [1.1 1]
%     bwmax  maximum receiver BW (full BW (+-))           [Hz]     500d3
%       gdt  Gradient update time                         [s]      4d-6
%      verb  print output + plot gradients & k-space      logical  true
%
%         k  k-space trajectory
%       dcf  density compensation function
%         t  time
%      grad  gradient trajectory    [#pts,#exc,direction]
%            directions: 1=phase; 2=freq
%       out  output structure of write_ak_wav.m
%
%  6/2020 Rolf Schulte
if (nargin<1), help(mfilename); return; end

fac_burst = 0.5;


%% input
fov = 1d-3*fov(:).';
if length(fov)~=1 && length(fov)~=2, error('length(fov)~=1 or 2'); end
if length(fov)==1,           fov = [fov fov]; end
if length(mtx)~=1 && length(mtx)~=2, error('length(mtx)~=1 or 2'); end
if length(mtx)==1,   mtx = [mtx mtx]; end
if any(rem(mtx,1)),          error('mtx not natural number'); end
if length(nexc)~=1,          error('length(nexc)~=1'); end
if rem(nexc,1),              error('nexc not natural number'); end
if nexc<1,                   error('nex(=%g)<1',nexc); end
if any(gmax<1),              warning('gmax<1: input [mT/m]'); end
gmax = 1d-3*gmax(:).';                                 % [T/m]
if length(gmax)~=1 && length(gmax)~=2, error('length(gmax)~=1 or 2'); end
if length(gmax)==1,          gmax = [gmax gmax]; end
if length(smax)~=1 && length(smax)~=2, error('length(smax)~=1 or 2'); end
if length(smax)==1,          smax = [smax smax]; end
if ~exist('rampsamp','var'), rampsamp = []; end
if isempty(rampsamp),        rampsamp = true; end
if ~exist('ref','var'),      ref = []; end
if isempty(ref),             ref = 0; end
if ~exist('do_burst','var'), do_burst = []; end
if isempty(do_burst),        do_burst = false; end
if ~exist('undsamp','var'),  undsamp = []; end
if isempty(undsamp),         undsamp = 0; end
if undsamp~=0, warning('undsamp not yet implemented'); end
if ~exist('fname','var'),    fname = []; end
if ~exist('nucleus','var'),  nucleus = []; end
if isempty(nucleus),         nucleus = '1H'; end
if ~exist('samp','var'),     samp = []; end
if isempty(samp),            samp = [1.1 1]; end
if ~exist('bwmax','var'),    bwmax = []; end
if isempty(bwmax),           bwmax = 500d3; end  % max Rx BW (full BW) [Hz]
if ~exist('gdt','var'),      gdt = []; end
if isempty(gdt),             gdt = 4e-6; end     % gradient update time [s]
if ~exist('verb','var'),     verb = []; end
if isempty(verb),            verb = true; end


%% misc parameters
gpre_derate = 0.85;          % factor to derate pre-phaser to reduce PNS peak
grad_fufa = 0.995;           % reduce gmax+smax to prevent digitisation error
fov_samp = fov.*samp;        % include over-/ undersampling via FOV
kmax = 2*pi*mtx./fov_samp;   % k-space extend [rad/m]
dt = 1/(fov_samp(1)*abs(gyrogamma(nucleus))/2/pi*gmax(1)*grad_fufa);
dki = 2^ceil(log2(dt/gdt));  % round down for oversampling
bw = 1/(gdt*dki);            % full acq BW [Hz]
if bw>bwmax                  % sampling rate limited regime
    fprintf('Attention: sampling rate limited regime\n');
    fprintf('Derating bw from %g [kHz] to %g [kHz]\n',bw*1d-3,bwmax*1d-3);
    bw = bwmax;
end
g_read = bw/(fov_samp(1)*abs(gyrogamma(nucleus))/2/pi); % grad strength [T/m]
n_koff = 4;                  % skip first acq points (artefact)
n_goff = ceil(n_koff*dki);   % #grad pts offset for start of traj

n_ramp = ceil(g_read/(smax(1)*grad_fufa)/gdt)+1;
smax_pre2 = smax(2)*gpre_derate;
smax_phas = smax(2)*grad_fufa;
if do_burst 
    if isodd(mtx(2)/nexc)
        error('mtx(2)(=%g)/nexc(=%g) is odd',mtx(2),nexc);
    end
    if nexc<2, warning('nexc(=%g)<2: increase for burst',nexc); end
end


%% phase encoding blips
ga_phas = 2*pi/fov_samp(2)/gyrogamma(nucleus)*nexc;
% G_phas = gradient_lobe_old(ga_phas,gdt,gmax(2)*grad_fufa,smax_phas,0,true);
% n_phas = size(G_phas,2);     % #pts of phase encoding blip
[G_phas,n_phas] = gradient_lobe(ga_phas,0,gmax(2)*grad_fufa,smax_phas,gdt,true);
ng2_pha = 2*n_ramp-n_phas;   % #pts on ramps w/o phase encoding
if ng2_pha<0, error('ng2_pha (=%g) < 0',ng2_pha); end
if abs(ng2_pha/2-floor(ng2_pha/2))>1d-10, error('ng2_pha not even'); end


%% readout parameters
n_plat = ceil(kmax(1)/(gyrogamma(nucleus)*g_read*gdt)); % pts/readout (plateau)
if (n_plat~=(mtx(2)+1))
    warning('n_plat(=%g)~=(mtx(2)(=%g)+1)',n_plat,mtx(2));
end
if rampsamp
    n_plat = ceil(n_plat - n_ramp + n_ramp/n_phas^2*2 + 1);
end
if n_plat<1, warning('n_plat(=%g)<1; setting to 1',n_plat); n_plat = 1; end
n_lobe = mtx(2)/nexc;        % #epi lobes
if rem(n_lobe,1)>1d-30
    warning('mtx(2)(=%g) not multiple of nexc(=%g); rounding up n_lobe(=%g)',...
        mtx(2),nexc,n_lobe);
    n_lobe = ceil(n_lobe);
end
if n_lobe<1, error('n_lobe(=%g)<1',n_lobe); end
n_pts = 2*n_ramp + n_plat;   % grad pts per lobe


%% print info
if verb
    fprintf('bw (full acq) = %g [kHz]\n',bw*1d-3);
    fprintf('gdt = %g [us]; kdt = %g [us]\n',gdt*1d6,1d6/bw);
    fprintf('g_read = %g [mT/m]\n',g_read*1d3);
    fprintf('n_plat = %g\n',n_plat);
    fprintf('n_lobe = %g\n',n_lobe);
end


%% prewinders (at beginning of trajectory)
% prewinder for frequency encoding direction
ga_pre1 = -g_read*(n_plat+n_ramp)*gdt/2;
[G_pre1,ng_pre] = gradient_lobe(ga_pre1,0,gmax(1)*gpre_derate,...
    smax(1)*gpre_derate,gdt);

% prewinder for phase-encoding direction
for l=1:nexc
    ga_pre2 = -kmax(2)/2/gyrogamma(nucleus) + ga_phas*(l-1)/nexc;
    % G_pre2{l} = gradient_lobe_old(ga_pre2,gdt,gmax(2)*gpre_derate,smax_pre2,0);
    G_pre2{l} = gradient_lobe(ga_pre2,0,gmax(2)*gpre_derate,smax_pre2,gdt);
    if size(G_pre2{l},2)>ng_pre, ng_pre = size(G_pre2{l},2); end
end
if do_burst
    ga_burst2 = fac_burst*kmax(2)/gyrogamma(nucleus) + ga_phas*(nexc+1)/nexc;
    for l=2:nexc
         % G_pre2{l} = gradient_lobe_old(ga_burst2,gdt,...
         %     gmax(2)*gpre_derate,smax_pre2,0);
         G_pre2{l} = gradient_lobe(ga_burst2,0,gmax(2)*gpre_derate,smax_pre2,gdt);
         if size(G_pre2{l},2)>ng_pre, ng_pre = size(G_pre2{l},2); end
    end
end


%% initialise
n1 = n_goff+ng_pre+(2*n_ramp+n_plat)*n_lobe+1; % #grad pts/excitation
n1 = ceil(n1/2)*2;           % round up to even number
grad = zeros(n1,nexc,2);     % initialise
ind = false(n1,1);           % index list for data acq


%% add prewinders
ii = n_goff+ng_pre-size(G_pre1,2)+(1:size(G_pre1,2))+1;
if do_burst
    grad(ii,1,1) = G_pre1.';                     % freq dir
else
    grad(ii,:,1) = repmat(G_pre1.',[1 nexc]);    % freq dir
end
for l=1:nexc, grad(n_goff+(1:size(G_pre2{l},2)),l,2) = G_pre2{l}; end


%% single frequency encoding lobe
g_freq = g_read*[(0:(n_ramp-1)).'/n_ramp; ...
    ones(n_plat,1); (n_ramp:-1:1).'/n_ramp];


%% assemble full gradient waveform
for l=1:n_lobe
    sgn = 2*rem(l,2)-1;                          % sign
    i_start = n_goff + ng_pre + (l-1)*n_pts;     % index to current lobe
    grad(i_start+(1:n_pts),:,1) = sgn*repmat(g_freq,[1 nexc]);  % freq
    
    if l<n_lobe
        ii = i_start+n_pts-n_phas/2+(1:n_phas);
        grad(ii,:,2) = repmat(G_phas.',[1 nexc]);% phase encoding
    end
    % ind(i_start+n_ramp+(1:(n_plat)),1) = true; % index w/o ramp sampling
    % ii = i_start + n_ramp + (1:(n_plat-1)) + ~isodd(l);
    ii = i_start + n_ramp + (1:n_plat);          % sample
    ind(ii,1) = true;                            % index w/o ramp sampling
end
if size(grad,1)~=n1
    warning('size(grad,1)(=%g)~=n1(=%g)',size(grad,1),n1); 
end
if do_burst
    grad = reshape(grad,[n1*nexc 1 2]);
end


%% convert to k-space trajectory + calculate index for reco
kk = cumsum(grad,1)*gdt*abs(gyrogamma(nucleus));
if do_burst
    koff_burst = (fac_burst+1)*kmax(2);
    ga_burst = koff_burst/gyrogamma(nucleus);
    for l=2:nexc
        ii = (1:n1) + n1*(l-1);
        kk(ii,1,2) = kk(ii,1,2) - koff_burst*(l-1);
    end
else
    ga_burst = [];
    
end
if dki>1
    k = kk(1:dki:end,:,:);
    ind = ind(1:dki:end,1);
else
    if abs(dki-0.5)<1d-10
        warning('interpolating k with splines');
        k = interp1((1:size(kk,1)),kk,(1:0.5:size(kk,1)+0.5),'spline');
        tmp = false(size(ind,1)*2,1);
        tmp(1:2:end,:) = ind;
        tmp(2:2:end,:) = ind;
        ind = tmp;
    else
        k = kk;
    end
end



%% misc indexing
[~,ii] = min(cumsum(G_pre1));
read_start = ceil((ii+n_goff)/dki);
if rampsamp        % overwrite previous index list
    ind = true(size(k,1),1);
    ind(abs(k(:,1,1))>kmax(1)/2,1) = false;
    ind(1:read_start,1) = false;
end
t = (0:size(k,1)-1)/bw;
if do_burst
    if ~rampsamp
        ind = repmat(ind,[nexc 1]).'; 
    else
        ind = ind.';
    end
end
k = k(ind,:,:);
t = t(ind).';
if ~do_burst
    ind = repmat(ind,[1 nexc]).';
    t = repmat(t,[1 nexc]);
end
ni2 = size(ind,2);
if abs(ni2/32-floor(ni2/32))>1d-12
    fprintf('size(ind,2)(=%g) not multiple of 32; adding false to end\n',ni2);
    ind = [ind , false(size(ind,1),(ceil(ni2/32)*32-ni2))];    
end


if verb
    kmipma = min(k(:,1,1)) + max(k(:,1,1));
    fprintf('min(k(:,1,1))=%g; max(k(:,1,1))=%g; min+max=%g\n',...
        min(k(:,1,1)),max(k(:,1,1)),kmipma);
    if abs(kmipma)>1d-10, warning('|min(k)+max(k)|(=|%g|)>1d-10',kmipma); end
    fprintf('#ind pts = %g; #ind pts/lobe = %g\n',...
        sum(ind(:)),sum(ind(end,:))/n_lobe);
end


%% print info + plot grad & k-space
gmax_act = [max(max(abs(grad(:,:,1)))) , max(max(abs(grad(:,:,2))))];
smax_act = [max(max(abs((diff(grad(:,:,1),[],1)/gdt)))) , ...
    max(max(abs((diff(grad(:,:,2),[],1)/gdt))))];
if verb
    fprintf('max grad = %g %g [mT/m]; specified gmax = %g %g [mT/m]\n',...
        gmax_act(1)*1d3,gmax_act(2)*1d3,gmax(1)*1d3,gmax(2)*1d3);
    fprintf('max slew = %g %g [mT/m]; specified smax = %g %g [mT/m]\n',...
        smax_act(1),smax_act(2),smax(1),smax(2));
    % fprintf('Discard first %g acq pts\n', read_start);
    fprintf('Duration = %g [ms]\n',n1*gdt*1d3);
    fprintf('size(ind) = %g %g\n',size(ind,1),size(ind,2));

    % plotting
    figure(11); clf;
    if do_burst
        nexc_plt = 1;
    else
        nexc_plt = nexc;
    end
    for l=1:nexc_plt
        subplot(nexc_plt,1,l); 
        tg = (0:size(grad,1)-1)*gdt*1d3;
        plot(tg,squeeze(grad(:,l,:)*1d3));
        axis([0 max(tg) -max(gmax)*1.02d3 max(gmax)*1.02d3]);
        grid on;
        if l==1, title('Gradient trajectory'); end
        if l==nexc, xlabel('time [ms]'); end
    end
    
    figure(12); clf;
    aa = 1.1*max(kk(:));
    plot([-1 -1 1 1 -1]*kmax(1)/2,[-1 1 1 -1 -1]*kmax(2)/2,'m',...
        kk(:,:,1),kk(:,:,2),'-',k(:,:,1),k(:,:,2),'x',...
        [-1 1]*aa,[0 0],'r:',[0 0],[-1 1]*aa,'r:');
    axis([-1 1 -1 1]*aa); grid on;
    title('k-space');
end


%% reshape k: (1,nk,2)
tmp1 = k(:,:,1); tmp2 = k(:,:,2);
k = zeros(1,size(k,1)*size(k,2),2);
k(1,:,1) = tmp1(:)/kmax(1); 
k(1,:,2) = tmp2(:)/kmax(2);


%% density compensation function
dcf = [abs(diff(k(1,:,1)))*mtx(1)/samp(2) 0];
% plot([1:size(k,2)],dcf,'b',[1:size(k,2)],dcf2,'r')


%% add EPI reference scan to front
if ref(1)>0
    ref = ref(:).';
    ref = [ref zeros(1,4-size(ref,2))];
    if length(ref)~=4
        error('length(ref)(=%g)~=2,3 or 4',length(ref)); 
    end
    if verb
        fprintf('Adding reference scan:\n');
        if ref(1)>0
            fprintf('\tnull phase encode=%g (EPI reference scan)\n',ref(1)); 
        end
        fprintf('\trepeat grad %g times\n',ref(2));
        if ref(3)>0
            fprintf('\tnull EPI train=%g (phase reference)\n',ref(3)); 
        end
        if ref(4)>0
            fprintf('\tnull phase encode=%g and invert EPSI\n',ref(4));
        end
    end
    grad = repmat(grad,[1 ref(2) 1]);
    
    % ref(1) -> EPI reference scan
    if ref(1)>0
        if do_burst
            ii2 = ref(1);                       % null ref(1) view
        else
            ii2 = (ref(1)-1)*nexc + (1:nexc);   % null ref(1) view
        end
        grad(:,ii2,2) = 0;                      % null phase encodes -> epsi
        if do_burst                             % add again burst gradients
            G_burst = gradient_lobe(ga_burst,0,gmax(2)*gpre_derate,smax_pre2,gdt);
            for l=2:nexc
                ii1 = n_goff+(1:size(G_burst,2))+(l-1)*n1;
                grad(ii1,ii2,2) = G_burst;
            end
        end
    end

    % ref(3) -> phase reference scan
    if ref(3)>0
        if do_burst
            ii2 = ref(3);                       % null ref(3) view
        else
            ii2 = (ref(3)-1)*nexc + (1:nexc);   % null ref(3) view
        end
        grad(:,ii2,2) = 0;                      % null phase encodes
        grad(:,ii2,1) = 0;                      % null EPSI train
        if do_burst                             % add again burst gradients
            G_burst = gradient_lobe(ga_burst,0,gmax(2)*gpre_derate,smax_pre2,gdt);
            for l=2:nexc
                ii1 = n_goff+(1:size(G_burst,2))+(l-1)*n1;
                grad(ii1,ii2,2) = G_burst;
            end
        end
    end

    % ref(4) -> inverted EPI reference scan
    if ref(4)>0
        if do_burst
            ii2 = ref(4);                       % null ref(4) view
        else
            ii2 = (ref(4)-1)*nexc + (1:nexc);
        end
        grad(:,ii2,2) = 0;                      % null phase encodes -> epsi
        if do_burst                             % add again burst gradients
            G_burst = gradient_lobe(ga_burst,0,gmax(2)*gpre_derate,smax_pre2,gdt);
            for l=2:nexc
                ii1 = n_goff+(1:size(G_burst,2))+(l-1)*n1;
                grad(ii1,ii2,2) = G_burst;
            end
        end
        grad(:,ii2,1) = -grad(:,ii2,1);          % invert EPSI
    end
    ref_str = sprintf('%g_%g',ref(1:2));
    if any(ref(3:4)>0)
        ref_str = sprintf('%s_%g_%g',ref_str,ref(3:4));
    end
else
    ref_str = '0';
end




%% generate filename
if (islogical(fname) || isnumeric(fname))
    if fname
        if isnumeric(nucleus), error('Please enter string for nucleus'); end
        fname = sprintf('epi_%s_fov%g_mtx%gx%g_nexc%g_gmax%g_smax%g_ref%s',...
            nucleus,round(max(fov)*1d3),mtx(1),mtx(2),...
            nexc,round(max(gmax_act)*1d3),round(max(smax_act)),ref_str);
        if rampsamp, fname = [fname '_rs']; end
        if do_burst, fname = [fname '_burst']; end
    else
        fname = [];
    end
end


%% echo-time
if verb
    TE = (n_goff+ng_pre+(2*n_ramp+n_plat)*(n_lobe/2+0.5))*gdt;
    % if do_burst, TE = TE + n1*(nexc-1)*gdt; end
    fprintf('TE += %g [ms]\n',TE*1d3);
end


%% OUTPUT
if ~isempty(fname)
    fprintf('Saving trajectory to %s\n',fname);
    fov1H = max(fov)/abs(gyrogamma('1h'))*abs(gyrogamma(nucleus));
    n_ind = size(ind,2); npix = max(mtx);
    % n_acq = ceil((n_ind+10)/2)*2; % add a few points for gradient delay
    
    
    desc = sprintf('EPI trajectories; fov=[%g %g], mtx=[%g %g], nexc=%g\n',...
        fov(1),fov(2),mtx(1),mtx(2),nexc);
    desc = sprintf('%sgmax=%g, smax=%g, rampsamp=%g, ref=%s, do_burst=%g\n',...
        desc,gmax,smax,rampsamp,ref_str,do_burst);
    desc = sprintf('%snucleus=%s, samp=[%g %g], bwmax=%g\n',...
        desc,nucleus,samp(1),samp(2),bwmax);
    desc = sprintf('%sgdt=%g, fov1H=%g, bw=%g, n_ind=%g\n',...
        desc,gdt,fov1H,bw,n_ind);    
    
    out = write_ak_wav([fname,'.wav'],grad,bw,fov1H,desc,n_ind);
    save(fname,'k','dcf','t','out',...
        'fov','mtx','nexc','gmax','smax','nucleus','fname',...
         'rampsamp','ref','do_burst','samp','bwmax','gdt',...
         'ind','fov1H','dcf','fov_samp','t','bw','desc','n_ind','npix');
else
    out = [];
end
if do_burst
    ga_burst = [0;ga_burst];        % y/phase encode for design_rf4burst
end

end    % main function design_epi.m
