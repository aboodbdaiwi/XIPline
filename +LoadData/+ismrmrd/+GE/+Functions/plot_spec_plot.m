function [spec] = plot_spec_plot(handles)
% PLOT_SPEC_PLOT  plot spectra
%      plot_spec_plot(handles)
%
% 7/2006 Rolf Schulte
% See also PLOT_SPEC.
if nargin<1, help(mfilename); return; end


%% misc parameters
spec = handles.spec;
ppm  = handles.ppm;
pc   = read_editable(handles.phi_txt);  pc = pc(:);
phi1 = read_editable(handles.phi1_txt); phi1 = phi1(:);
turn = read_editable(handles.turn_txt); turn = turn(:);
zp   = read_editable(handles.zp_txt);   zp = zp(:);
lne  = read_editable(handles.line_txt);
si = size(spec);
figure(handles.figid);
switch get(handles.axis_unit,'Value')
 case 1, unit = 1;
 case 2, unit = handles.f0;
 otherwise, error('axis_unit not existing');
end
ax = (ppm+zp(1))*unit;


%% shift zp individually by linear phase in t-domain
if length(zp)>1,
    t = (0:si(2)-1)/handles.bw;
    for l1=2:si(1),
        fid = ifft(ifftshift(spec(l1,:),2),[],2);	% recon fid
        % move zp in t-domain
        zp_shift = (zp(l1)-zp(1))*handles.f0;
        spec(l1,:) = fftshift(fft(fid.*exp(1i*2*pi*t*zp_shift),[],2),2);
    end
end


%% retrieve zooming
zon = false;
x_lim = get(gca,'XLim');
y_lim = get(gca,'YLim');
apd = getappdata(handles.figid);
if isfield(apd,'ZoomOnState'),
  if strcmpi(apd.ZoomOnState,'on'),
    zon = true;
    apda = getappdata(gca);
    if isfield(apda,'matlab_graphics_resetplotview'),
      mgrpv = apda.matlab_graphics_resetplotview;
    end
  end
end
clf;


%% 0th order phase correction
if isempty(pc),
    spec = abs(spec);
else
    if any(pc~=0),
        if length(pc)~=si(1),
            spec = spec*exp(-1i*pc(1)/180*pi);
        else
            spec = spec.*(exp(-1i*pc/180*pi)*ones(1,si(2)));
        end
    end
end


%% 1st order phase correction
if any(phi1~=0),
    if (length(phi1)==1) && (length(turn)==1)
        spec = spec.*(ones(si(1),1)*exp(-1i*(ppm+zp(1)-turn)*phi1));
    else
        if length(phi1)~=si(1), phi1 = phi1*ones(si(1),1); end
        if length(turn)~=si(1), turn = turn*ones(si(1),1); end
        spec = spec.*exp(-1i*(ones(si(1),1)*ppm+zp(1)-turn*ones(1,si(2))).* ...
            (phi1*ones(1,si(2))));
    end
end

if get(handles.baseline,'Value')>1
    base = zeros(size(spec));
    navg = ceil(si(2)/512)*2;     % number of averages
    ind = [navg/2 si(2)/4 si(2)*3/4 si(2)-navg/2];
    for l1=1:si(1),
        ss = zeros(size(ind));
        for l2=1:length(ind)
            ss(l2) = mean(spec(l1,(ind(l2)-navg/2+1):(ind(l2)+navg/2)),2);
        end
        base(l1,:) = spline(ppm(ind),ss,ppm);
    end
    spec_cor = spec - base;
    if false
        % spec_cor = (linspace(0,40,si(1)).'*ones(1,si(2)))+spec_cor;
        spec_cor = (linspace(0,80,si(1)).'*ones(1,si(2)))+spec_cor;
    end
end


%% determine axis (without zooming)
dig  = floor(log10(max(max(abs(real(spec))))));
ymin = floor(min(min(real(spec)))/10^(dig-1)-1)*10^(dig-1);
ymax = ceil(max(max(real(spec)))/10^(dig-1)+1)*10^(dig-1);
xmin = min(ax);
xmax = max(ax);

linewidth = 1;
switch get(handles.baseline,'Value')
    case 1,
        plot(ax,real(spec),'',...
            [xmin xmax],[0 0],'r:','LineWidth',linewidth);
    case 2,
        plot(ax,real(spec_cor),'',...
            [xmin xmax],[0 0],'r:','LineWidth',linewidth);
    case 3,
        plot(ax,real(spec),'',...
            [xmin xmax],[0 0],'r:','LineWidth',linewidth);
        hold on; plot(ax,real(base),':');
        plot(ax,real(spec_cor),'--','LineWidth',linewidth); 
        hold off;
end
hold on; 
plot([turn(1) turn(1)],[ymin ymax],'g:');
plot([lne lne],[ymin ymax],'c:');
hold off;

if zon,
  % the 2d zoom of matlab saves the original axis in zlabel
  % 3d zoom works over camera positions
  set(gca,'XLim',x_lim);		% set the axes
  set(gca,'YLim',y_lim);		% 
  if exist('mgrpv','var'),
    setappdata(gca,'matlab_graphics_resetplotview',mgrpv);
  end
  zoom on;				% switch on zoom button
else
  axis([xmin xmax ymin ymax]);
end

set(gca,'XDir','reverse');
