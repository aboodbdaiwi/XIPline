function [PThresh,pt,PTmax,gmax,smax,t,f] = pns(grad,chronaxie,rheobase,alpha,gdt,plt)
%PNS Convolution model to calculate peripheral nerve stimulation
%[PThresh,pt,PTmax,gmax,smax,t,f] = pns(grad,chronaxie,rheobase,alpha,gdt,plt)
%                                           (CV par)       [unit] (default)   
%      grad  gradient matrix (dimensions,points,repetitions) [T/m]
%            (for 2D also dim=1 + complex data possible)
%            also possible: (points,repetitions,dimensions)
% chronaxie  Chronaxie time constant        (cfrfact [us]) [s]
%  rheobase  Rheobase scale factor          (cfrinf [T/s]) [T/s]
%     alpha  Effective coil length     (cfdbdtdx/y/z [cm]) [m]
%            SRmin=rheobase/alpha
%       gdt  Gradient update time                          [s] (4d-6)
%       plt  Plotting + print output                           (true)
%
%   PThresh  Threshold of PNS               (cfdbdtper)    [%]
%        pt  Threshold of different axis (uncombined)
%     PTmax  Maximum PNS threshold (combined)
%      gmax  Maximum gradient strength (combined)          [T/m]
%      smax  Maximum slewrate (combined)                   [T/m/s]
%
%alternatively:
%  pns(grad,coil)
%      coil  XRM,XRMW,WHOLE,ZOOM,HRMW,HRMB
%
% Scanner  Gradient coil   chronaxie rheobase alpha  gmax   smax
%                          [s]       [T/s]    [m]    [mT/m] [T/m/s]
% MR750    XRM             334d-6    23.4     0.333  50     200
% MR750w   XRMW            360d-6    20.0     0.324  33     120
% HDx      TRM WHOLE       370d-6    23.7     0.344  23     77
% HDx      TRM ZOOM        354d-6    29.1     0.309  40     150
% Rio      HRMW            642.4d-6  17.9     0.310  70     150
% 7T       HRMB            359d-6    26.5     0.370  100    200
%
% values on scanner from /w/config/Scandbdt.cfg (chronaxie & rheobase) + 
%  /export/home/mx/host/config/current/GRSubsystemHWO.xml (alpha=EffectivedBdTlengthX/100
%  /w/config/GradientConfig.cfg (gmax=xFSAmp*10 & smax=SRMode)
%
% Literature: "Introduction to dB/dt" presentation by Toni Linz; 6/8/06
%             DIN EN 60601-2-33
% 8/2020 Rolf Schulte
if nargin<1, help(mfilename); return; end

model = 2;         % 1=convolution in time domain; 2=multiplication in FD


%% input checking
if length(size(grad))==3
    [n1,n2,n3] = size(grad);
    if (n1>n3) && (n3<4)
        fprintf('Shifting grad dimension:\n');
        fprintf('\tassuming (points,repetitions,dimensions)\n');
        grad = shiftdim(grad,2);
    end
end
if ~isreal(grad)
    if size(grad,1)~=1, error('grad complex and size(grad,1)~=1'); end
    grad = [real(grad) ; imag(grad)];
end
if any(isnan(grad(:))), warning('grad contains NaN'); end
if any(isinf(grad(:))), warning('grad contains inf'); end

if ischar(chronaxie)
    switch lower(chronaxie)
        case 'xrm',   chronaxie=334d-6; rheobase=23.4; alpha=0.333;
        case 'xrmw',  chronaxie=360d-6; rheobase=20.0; alpha=0.324;
        case 'whole', chronaxie=370d-6; rheobase=23.7; alpha=0.344;
        case 'zoom',  chronaxie=354d-6; rheobase=29.1; alpha=0.309;
        case 'hrmw',chronaxie=642.4d-6; rheobase=17.9; alpha=0.310;
        case 'hrmb',  chronaxie=359d-6; rheobase=26.5; alpha=0.370;
        otherwise, error('gradient coil (%s) unkown',chronaxie);
    end
end

if ~exist('rheobase','var'), error('Provide rheobase'); end
if ~exist('alpha','var'), error('Provide alpha'); end
if ~exist('gdt','var'), gdt = []; end
if isempty(gdt), gdt = 4d-6; end
if ~exist('plt','var'), plt = true; end
if (chronaxie>800d-6) || (chronaxie<200d-6)
    warning('chronaxie=%g; typical values in help',chronaxie);
end
if (rheobase>30) || (rheobase<10)
    warning('rheobase=%g; typical values in help',rheobase);
end
if (alpha>0.4) || (alpha<0.3)
    warning('alpha=%g; typical values in help',alpha);
end


%% zerofilling for FFT (must decay to zero to avoid wiggles in front)
if model==2
    if isodd(size(grad,2)), grad(:,size(grad,2)+1,:) = 0; end
    zf2 = ceil(7d-3/gdt/2)*2;
    dozf = false;
    if size(grad,2)<zf2
        dozf = true;
    else
       if any(any(grad(:,end-zf2+1:end,1)~=0)), dozf = true; end
    end
    if dozf, grad(:,size(grad,2)+zf2,:) = 0; end
end

if length(size(grad))>3, error('length(size(grad))>3'); end
if length(size(grad))<2, error('length(size(grad))<2'); end
[n1,n2,n3] = size(grad);
if n1>3, error('size(grad,1)>3'); end
if n2==1, error('size(grad,2)==1'); end
if n3>n2, warning('size(grad,3)>size(grad,2)'); end


%% misc parameters
SRmin = rheobase/alpha;      % min slewrate causing stimulation; Slide 3 Eq 3+4 [T/m/s]
SR = [zeros(n1,1,n3),diff(grad,[],2)]/gdt;  % instantaneous slewrate [T/m/s]

switch model
    case 1    % convolution in time domain (identical to calcpns)
        pt = zeros(n1,n2,n3);
        for l3=1:n3
            for l2=1:n2
                % Slide 34 Eq 30
                pt(:,l2,l3) = 100*gdt*chronaxie/SRmin*...
                    sum(SR(:,(1:l2),l3)./(chronaxie+repmat((l2:-1:1)-0.5,[n1,1])*gdt).^2,2);
            end
        end
    case 2    % multiplication in freq domain
        % small difference exist because time is shifted 1/2*gdt in time

        decay = fftshift(fft(1./(((0:n2-1)+0.5)*gdt+chronaxie).^2));
        % decay = fftshift(fft(1./(((0:n(2)-1))*gdt+chronaxie).^2));
        SR_fd = fftshift(fft(SR,[],2),2);
        pt = 100*gdt*chronaxie/SRmin*...
            ifft(ifftshift(SR_fd.*repmat(decay,[n1,1,n3]),2),[],2);
    otherwise
        error('model (%g) unknown',model);
end

if any(isnan(pt(:))), warning('pt contains NaN'); end
if any(isinf(pt(:))), warning('pt contains inf'); end

PThresh = sqrt(sum(pt.^2,1));          % percent threshold; Slide 35; Eq 31


%% plotting
if plt
    % lp = 30;
    [~,lp] = max(max(PThresh,[],2),[],3);   % index to plot worst trajectory
    t = (0:n2-1)*gdt*1d3;                   % time [ms]
    figure(13); 
    clf;
    rsos_grad = sqrt(sum(grad(:,:,lp).^2));
    subplot(2,2,1); plot(t,grad(:,:,lp)*1d3,'',t,rsos_grad*1d3,'r--');
    xlabel('time [ms]'); ylabel('grad [mT/m]');
    tmp = ceil(max(abs(rsos_grad))*1050);
    grid on; axis([0 t(end) -tmp tmp]);
    
    rsos_SR = sqrt(sum(SR(:,:,lp).^2,1));
    subplot(2,2,3); plot(t,SR(:,:,lp),'',t,rsos_SR,'r--');
    xlabel('time [ms]'); ylabel('slewrate [T/m/s]');
    tmp = ceil(max(abs(rsos_SR))*1.05);
    grid on; axis([0 t(end) -tmp tmp]);
    
    subplot(2,2,4); 
    plot(t,pt(:,:,lp),'',t,PThresh(:,:,lp),'r--',...
        [t(1) t(end)],[100 100],'m:',[t(1) t(end)],[80 80],'m:',...
        [t(1) t(end)],-[100 100],'m:',[t(1) t(end)],-[80 80],'m:');
    xlabel('time [ms]'); ylabel('PNS threshold [%]');
    tmp = 1.05*max([PThresh(:,:,lp) 100]);
    grid on; axis([0 t(end) -tmp tmp]);
    
    if model==2
        f = (-n2/2:(n2-1)/2)/n2/gdt*1d-3;
        subplot(2,2,2);
        plot(f,abs(SR_fd(:,:,lp))/max(max(abs(SR_fd(:,:,lp)))),'',...
            f,abs(decay)/max(abs(decay(:))),'r');
        xlabel('freq [kHz]'); ylabel('FFT of decay + slewrate');
        grid on; axis([-5 5 -0.05 1.05]);
    end
end


%% print output
if plt
    fprintf('\nchronaxie = %g [us]\n',chronaxie*1d6);
    fprintf('rheobase  = %g [T/s]\n',rheobase);
    fprintf('effective coil length: alpha = %g [cm]\n',alpha*1d2);
    fprintf('Stimulation slewrate:  SRmin = %g [T/m/s]\n',SRmin);
    fprintf('Gradient update time:    gdt = %g [us]\n',gdt*1d6);
    fprintf('Max gradient strength:  gmax = %g [mT/m]\n',max(grad(:))*1d3);
    fprintf('Maximum slewrate:       smax = %g [T/m/s]\n',max(SR(:)));

    if n3==1, fprintf('Maximum PThresh = %.4g [%%]\n',max(abs(PThresh)));
    else
        fprintf('Maximum PThresh = \n');
        for l3=1:n3, fprintf('\tl3=%g: %.4g [%%]\n',l3,max(abs(PThresh(:,:,l3)))); end
        fprintf('\tall: %.4g [%%]\n',max(abs(PThresh(:))));
    end
    
    if max(PThresh)>80
        if max(PThresh)>100
            fprintf('Warning: PThresh exceeding first controlled mode (100%%)!!!\n');
        else
            fprintf('Warning: PThresh exceeding normal mode (80%%)!\n');
        end
    end
end


%% misc other output arguments
if nargout>2, PTmax = max(PThresh(:)); end
if nargout>3
    if ~exist('rsos_grad','var'), rsos_grad = sqrt(sum(grad(:,:,1).^2,1)); end
    gmax = max(rsos_grad); 
end
if nargout>4
    if ~exist('rsos_SR','var'), rsos_SR = sqrt(sum(SR(:,:,1).^2,1)); end
    smax = max(rsos_SR); 
end

end           % main function pns.m
