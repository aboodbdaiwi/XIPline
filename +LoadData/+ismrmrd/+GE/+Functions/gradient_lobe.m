function [grad,npts] = gradient_lobe(ga_des,g_off,gmax,smax,gdt,rnd2even,verb)
% GRADIENT_LOBE Calculate gradient lobe with start and end values and 
% specific areaa. Results in triangluar or trapezoidal gradient shape.
%
%[grad,npts] = gradient_lobe(ga_des,g_off,gmax,smax,gdt,rnd2even,verb)
%                                                     [unit]    (default)
%     ga_des   Desired area of gradient lobe          [s*T/m]
%      g_off   Offset gradient values ([start end])   [T/m]     ([0 0])
%       gmax   Maximum gradient strength              [T/m]
%       smax   Maximum slew rate                      [T/m/s]
%        gdt   Sampling dwell time                    [s]       (4d-6)
%   rnd2even   Round up to even #waveform pts                   (false)
%       verb   0=off; 1=verbose; 2=verbose+plotting             (0)
%
%       grad   Gradient waveform [1,#pts]             [T/m]
%       npts   #grad waveform points
%
% 7/2020  Rolf Schulte
if (nargin<1), help(mfilename); return; end


%% input parameters
if ~exist('g_off','var'),    g_off = []; end
if isempty(g_off),           g_off = [0 0]; end
switch length(g_off)
    case 1, g_off = [g_off 0]; 
    case 2
    otherwise, warning('length(g_off)(=%d)~=1 or 2',length(g_off));
end
if ~exist('gmax','var'),     gmax = []; end
if isempty(gmax),            gmax = 33d-3; end
if ~exist('smax','var'),     smax = []; end
if isempty(smax),            smax = 120; end
if ~exist('rnd2even','var'), rnd2even = []; end
if isempty(rnd2even),        rnd2even = false; end
if ~exist('verb','var'),     verb = []; end
if isempty(verb),            verb = 0; end
if ~exist('gdt','var'),      gdt = []; end
if isempty(gdt),             gdt = 4d-6; end


%% input checks
if length(ga_des)~=1, error('length(ga_des)~=1'); end
if ~isreal(ga_des)
    warning('ga_des is complex number; taking abs');
    ga_des = abs(ga_des);
end
if gmax<0, error('gmax(=%g)<0',gmax); end
if smax<0, error('smax(=%g)<0',smax); end
if any(abs(g_off)>gmax)
    warning('abs(g_off)(=[%g %g])>gmax(=%g)',...
        abs(g_off(1)),abs(g_off(2)),gmax);
end
if abs(ga_des/gmax/gdt)>1d10
    error('abs(ga_des/gmax/dt(=%g))>1',ga_des/gmax/gdt);
end
if abs(ga_des/gmax/gdt)>1d6
    warning('abs(ga_des(=%g))>1d-4',ga_des/gmax/gdt);
end
if abs(g_off(1)-4d-6)<1d-10, warning('input changed; flip gdt and g_off'); end


%% calculate positive lobes; invert when necessary
sgn    = sign(ga_des);
if abs(sgn)<1d-10, sgn = 1; end
ga_des = abs(ga_des);
g_off  = sgn*g_off;

g_nom = sqrt(smax*abs(ga_des)+sum(0.5*g_off.^2)); % nominal gradient strength
% if g_nom<1.1*sum(g_off) % if offset larger, invert again for lobe pointing up
if verb>2, fprintf('g_nom1=%g\n',g_nom*1d3); end
if any(g_nom<1.1*g_off)
    if verb>0, fprintf('Inverting\n'); end
    sgn = -sgn;
    ga_des = -ga_des;
    g_off = -g_off;
    % g_nom = sqrt(-smax*abs(ga_des)+sum(0.5*g_off.^2));
    g_nom = sqrt(abs(-smax*abs(ga_des)+sum(0.5*g_off.^2)));
end
if ~isreal(g_nom), error('~isreal(g_nom)'); end
if verb>0, fprintf('g_nom=%g; g_off=[%g %g][mT/m]\n',g_nom*1d3,g_off*1d3); end


%% differentiate between trapezoids and triangles
if g_nom>gmax      % trapezoid: add gradient plateau in middle
    if verb>0, fprintf('Constructing trapezoid\n'); end
    
    % n_ramp = ceil((gmax-g_off)/smax/gdt);             % #ramp points    
    n_ramp = ceil(abs(gmax-g_off)/smax/gdt);             % #ramp points    
    ga_ramp = 0.5*(gmax^2-g_off.^2)/smax;             % ramp areas    
    n_plt = ceil((ga_des-sum(ga_ramp))/gmax/gdt);     % #plateau points
    if any(n_ramp<2)
        warning('n_ramp(=[%d %d])<2',n_ramp(1),n_ramp(2));
        n_ramp(n_ramp<2) = 2;
    end
else          % triangular gradient
    if verb>0, fprintf('Constructing triangle\n'); end
    % n_ramp = ceil((g_nom-g_off)/smax/gdt)+1;
    n_ramp = ceil(abs(g_nom-g_off)/smax/gdt)+1;
    n_plt = 0;
    if any(n_ramp<2)
        warning('n_ramp(=[%d %d])<2',n_ramp(1),n_ramp(2));
        n_ramp(n_ramp<2) = 2;
    end
end
nn = sum(n_ramp)+n_plt;
if rnd2even && (abs(nn/2-floor(nn/2))>eps)
    if verb>0, fprintf('Rounding up: n_plt+=1\n'); end
    n_plt = n_plt+1;
    nn = nn+1;
end


%% calculate exact g_nom_act to prevent rounding errors
lit = 0;
while true
    lit = lit+1;
    [smax1,g_nom_act] = sub_calc(n_ramp,n_plt,ga_des,g_off,gdt);
    
    if lit==100, warning('#iterations exceeded'); break; end
    if all(smax1<smax)
        break;
    else
        % increase n_ramp or n_plt if slewrate exceeded
        iis = smax1>smax;
        smax2 = zeros(4,2);
        smax2(1,:) = sub_calc(n_ramp,n_plt+2,ga_des,g_off,gdt);
        smax2(2,:) = sub_calc(n_ramp+[0 2],n_plt,ga_des,g_off,gdt);
        smax2(3,:) = sub_calc(n_ramp+[2 0],n_plt,ga_des,g_off,gdt);
        smax2(4,:) = sub_calc(n_ramp+[1 1],n_plt,ga_des,g_off,gdt);
        if verb>2
            fprintf('smax: 1 %.4g %.4g; 2 %.4g %.4g ;3 %.4g %.4g ;4 %.4g %.4g; 5 %.4g %.4g\n',...
                smax1,smax2(1,:),smax2(2,:),smax2(3,:),smax2(4,:));
            fprintf('n_ramp = [%d %d]; n_plt = %d\n',n_ramp,n_plt);
        end
        if all(smax2(1,iis)<smax1(iis))
            n_plt = n_plt+2;
        else
            if all(smax2(2,iis)<smax1(iis))
                n_ramp = n_ramp+[0 2];
            else
                if all(smax2(3,iis)<smax1(iis))
                    n_ramp = n_ramp+[2 0];
                else
                    if all(smax2(4,iis)<smax1(iis))
                        n_ramp = n_ramp+[1 1];
                    else
                        warning('nothing helps to reduce smax');
                        break;
                    end
                end
            end
        end
        nn = nn+2;
    end
end
if verb>0, fprintf('g_nom_act=%g [mT/m]\n',g_nom_act*1d3); end


%% construct actual gradient waveform
grad = sgn*[linspace(g_off(1),g_nom_act,n_ramp(1)),...
    g_nom_act*ones(1,n_plt),...
    linspace(g_nom_act,g_off(2),n_ramp(2))];


%% check waveform
ga_act = sum(grad)*gdt;                     % actual area
ga_err = (sgn*ga_des-ga_act)/(sgn*ga_des);  % area difference
if (abs(ga_err)>1d-6) && (abs(ga_act)>1d-10)
    warning('abs(ga_err(=%g%%))>1d-6: ga_des=%g; ga_act=%g',...
        ga_err*100,ga_des,ga_act); 
end
smax_act = max(abs(diff(grad)))/gdt;        % actual slewrate [T/m/s]
gmax_act = max(abs(grad));                  % actual grad amp [T/m]
if smax_act>smax, warning('smax_act(=%g)>smax(=%g)',smax_act,smax); end
if gmax_act>gmax, warning('gmax_act(=%g)>gmax(=%g)',gmax_act,gmax); end
npts = size(grad,2);                   % #grad points in time-dimension
if nn~=npts, warning('nn(=%d)~=npts(=%d)',nn,npts); end
if abs(grad(1,1)-sgn*g_off(1))>1d-10
    warning('grad(1,1)(=%g)-g_off(1)(=%g)',grad(1,1),sgn*g_off(1));
end
if abs(grad(1,end)-sgn*g_off(2))>1d-10
    warning('grad(1,end)(=%g)-g_off(2)(=%g)',grad(1,end),sgn*g_off(2));
end


%% print info
if verb>0
    fprintf('gmax=%g;\tgmax_act=%g [mT/m]\n',gmax*1d3,gmax_act*1d3);
    fprintf('smax=%g;\tsmax_act=%g [T/m/s]\n',smax,smax_act);    
    fprintf('n_ramp=[%d %d]; n_plt=%d\n',n_ramp,n_plt);
    fprintf('ga_des=%g; ga_act=%g; err=%g[%%]\n',...
        sgn*ga_des,ga_act,ga_err*100);
    fprintf('duration=%g [ms]\n',npts*gdt*1d3);
end
if rnd2even && (abs(size(grad,2)/2-floor(size(grad,2)/2))>eps)
    warning('waveform not even');
end

%% plotting
if verb>1
    clf
    subplot(2,1,1);
    plot((1:npts),grad*1d3,'b',...
        [1 npts],-gmax*1d3*[1 1],'r:',[1 npts],gmax*1d3*[1 1],'r:'); 
    grid on; axis([1 npts 1.1d3*gmax*[-1 1]]);
    ylabel('grad [mT/m]');
    
    subplot(2,1,2);
    plot((1:npts-1),diff(grad)/gdt,'b',...
        [1 npts-1],-smax*[1 1],'r:',[1 npts-1],smax*[1 1],'r:');
    grid on; axis([1 npts-1 1.1*smax*[-1 1]]);
    ylabel('slewrate [T/m/s]');
    xlabel('points');
end


end      % gradient_lobe.m


%% sub-function
function [smax_calc,g_nom_act]=sub_calc(n_ramp,n_plt,ga_des,g_off,gdt)
g_nom_act = (2*ga_des/gdt-sum(n_ramp.*g_off))/(sum(n_ramp)+2*n_plt);
smax_calc = abs((g_off-g_nom_act)./((n_ramp-1)*gdt));
end
