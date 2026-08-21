function [g_fresp,hz,val,esp_time,esp] = check_grad_acoustics(grad,esp,verb,gdt)
%CHECK_GRAD_ACOUSTICS  Check acoustic response of gradient trajectory 
%[g_fresp,hz,val,esp_time,esp] = check_grad_acoustics(grad,esp,verb,gdt)
%     grad  Gradient waveforms with 4us time resolution    [T/m]
%           [#pts/interleave,#interleaves,#groups]
%           with: #groups = 2 for 2d-imaging, =3 for 3d-imaging
%           if #groups=1 and complex grad -> 2D
%      esp  Forbidden echo spacings (cell array or string)
%           string with name of gradient coil: 'xrmw','xrm','hrmw'
%           values from /usr/g/bin/epiesp.dat, epiespW.dat, epiespHRMw.dat
%     verb  Verbose (=1&2) and plotting (=2); figureID=verb(2)    ([2,102])
%    figid  Figure ID for plotting                                    (102)
%      gdt  Gradient update time                           [s]       (4d-6)
%
%  g_fresp  Frequency response of gradient trajectories
%       hz  Axis                                           [Hz]
%      val  Maximum of freq response in forbidden zones
% esp_time  Echo spacing (full period; EPI; =2*esp)        [us]
%
%  8/2024 Rolf Schulte
if nargin<1, help(mfilename); return; end


%% input parameters
if ~exist('esp','var'), esp = []; end
if isempty(esp),        esp{3} = []; end
if ~exist('verb','var'),verb = []; end
if isempty(verb),       verb = 2; end
if length(verb)<2,      verb = [verb,102]; end
if ~exist('gdt','var'), gdt = []; end
if isempty(gdt),        gdt = 4d-6; end


%% load from wav if grad is filename
fname = '';
if ischar(grad)
    if exist(grad,'file')
        fname = grad;
        fprintf('read_ak_wav(''%s'')',fname)
        [grad,~,~,desc,~,params] = read_ak_wav(fname);
        fprintf('desc=%s\n',desc);
        gdt = params(6)*1d-6;
    else
        error('ischar(grad(=''%s''); but file not found',grad);
    end
end


%% default parameters
maxfreq = 3d3;          % maximum frequency range for plotting [Hz]
zf_fac = 7;             % zero-filling factor
threshold = 0.3;        % maximum allowed in frequency band
grad = grad*1d3;        % change from [T/m] to [mT/m]


%% input checking
[n1,n2,n3] = size(grad);
if ~isreal(grad)
    if size(grad,3)~=1, error('grad complex and size(grad,3)~=1'); end
    tmp = zeros(n1,n2,2);
    tmp(:,:,1) = real(grad);
    tmp(:,:,2) = imag(grad);
    grad = tmp;
    clear tmp
end
if any(isnan(grad(:))), warning('grad contains NaN'); end
if any(isinf(grad(:))), warning('grad contains inf'); end


%% bands
% epiespXXX.dat file location: /usr/g/bin (old); /srv/nfs/psd/etc (new)
if ischar(esp)
    switch lower(esp)
        case {'xrmw','vrmw'}      % MR750w: epiespW.dat & epiespVRWw.dat
            tmp{3} = [330 460 0.0];
        case 'xrm'                % MR750:  epiesp.dat
            tmp{1} = [410 510 0.0];
            tmp{2} = [410 510 0.0];
            tmp{3} = [360 440 0.0];
        case 'hrmw'               % Premier: epiespHRMw.dat
            tmp{1} = [420 459 0.0 ; 816 865 1.6];
            tmp{2} = [420 459 0.0 ; 816 876 1.6];
            tmp{3} = [];
        case 'hrmb'               % new 7T: epiespHRMb.dat
            tmp{1} = [393 445 0.0; 476 528 0.0; 952 988 0.0];
            tmp{2} = [393 445 0.0; 476 504 0.0];
            tmp{3} = [393 445 0.0; 476 491 0.0];
        case 'hrmbuhp'            % ??: epiespHRMbUHP.dat
            tmp{1} = [350 445 0.0; 481 488 0.0];
            tmp{2} = [350 445 0.0; 481 488 0.0];
            tmp{3} = [418 426 0.0];
        case 'irmw'               % ??: epiespIRMw.dat
            tmp{1} = [827 896 0.0];
            tmp{2} = [827 896 0.0];
            tmp{3} = [];

        otherwise, error('gradient coil (%s) unkown',esp);
    end
    esp = tmp;
else
    if ~iscell(esp), error('esp must be char or cell'); end
end
bands = esp;
grad_axis = {'x','y','z'};
clr = {'r','g','b'};
for l1=1:length(bands)
    if verb(1)>0, fprintf('%s axis (colour=%s)\n',grad_axis{l1},clr{l1}); end
    for l2=1:size(bands{l1},1)
        bands{l1}(l2,2) = 1/(2d-6*esp{l1}(l2,1));
        bands{l1}(l2,1) = 1/(2d-6*esp{l1}(l2,2));
        if verb(1)>0
            fprintf('esp =   %g %g [us] %g\n',esp{l1}(l2,:));
            fprintf('bands = %.4g %.4g [Hz] %g\n',bands{l1}(l2,:));
        end
    end
end


%% zerofilling for FFT (must decay to zero to avoid wiggles in front)
zf1 = zf_fac*n1;
if isodd(n1+zf1), zf1 = zf1+1; end
grad(n1+zf1,:,:) = 0; 


if length(size(grad))>3, error('length(size(grad))>3'); end
if length(size(grad))<2, error('length(size(grad))<2'); end
[n1,n2,n3] = size(grad);
if n1==1, error('size(grad,1)==1'); end


%% calculating actual frequency response of gradient
g_fresp = ifftshift(ifft(grad,[],1),1);
bw = 1/gdt;
hz = (-n1/2:n1/2-1)/n1*bw;


%% check + print info
[magf,iigf] = max(max(max(abs(g_fresp),[],2),[],3),[],1);
freq_gf = abs(hz(iigf));
esp_time = round(1d6/freq_gf/2)*2;
% iigf = find(abs(g_fresp)==magf);
if verb(1)>0
    fprintf('Maximum complete g_fresp=%g at f=%g [Hz]; esp=%g[us]\n',...
        magf,freq_gf,esp_time/2);
    fprintf('Maximum of g_fresp in bands\n');
end
val = [];
for lg=1:n3
    for l1=1:length(bands)
        for l2=1:size(bands{l1},1)
            ind = (hz>=bands{l1}(l2,1)) & (hz<=bands{l1}(l2,2));
            if sum(ind)==0, warning('sum(ind)(=%g)==0',sum(ind)); end
            agb = abs(g_fresp(ind,:,lg));
            magb = max(max(agb));
            rmsagb = sqrt(mean(agb(:).^2));
            if verb(1)>0
                fprintf('lg=%g l1=%g l2=%g: magb=%g\n',lg,l1,l2,magb);
            end
            if magb>threshold
                warning('magb(=%g)>threshold(=%g)',magb,threshold); 
            end
            val{lg}{l1}{l2} = [magb rmsagb];
        end
    end
end


%% plotting
if verb(1)>1
    figure(verb(2)); clf
    if magf==0, warning('maf(=%g)==0; setting to 1',magf); magf = 1; end
    for lg=1:n3
        subplot(n3,1,lg);
        plot(hz,abs(g_fresp(:,:,lg)))
        hold on
        for l1=1:length(bands)
            ba = bands{l1};
            for l2=1:size(ba,1)
                h = fill([ba(l2,1) ba(l2,2) ba(l2,2) ba(l2,1)],...
                    [threshold threshold magf magf],clr{l1});
                h.FaceAlpha = 0.15;
            end
        end
        hold off
        xlabel('freq [Hz]');
        ylabel([grad_axis{lg} ' axis'])
        grid on;
        axis([-maxfreq maxfreq 0 magf]);
    end
    if ~isempty(fname), set(verb(2),'name',fname); end
end

end      % check_grad_acoustics.m
