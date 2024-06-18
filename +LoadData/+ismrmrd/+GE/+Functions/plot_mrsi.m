function varargout = plot_mrsi(varargin)
%PLOT_MRSI Plot GUI for MRSI data
% plot_mrsi(spec,par)
%  spec  MRSI spectra: matrix size = (spec,x,y,z)
%   par  Parameter structure (pc,verb,fax,fig_str)
%
% 7/2020 Rolf Schulte
if nargin<1, help(mfilename); return; end


% PLOT_MRSI MATLAB code for plot_mrsi.fig
%      PLOT_MRSI, by itself, creates a new PLOT_MRSI or raises the existing
%      singleton*.
%
%      H = PLOT_MRSI returns the handle to a new PLOT_MRSI or the handle to
%      the existing singleton*.
%
%      PLOT_MRSI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PLOT_MRSI.M with the given input arguments.
%
%      PLOT_MRSI('Property','Value',...) creates a new PLOT_MRSI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before plot_mrsi_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to plot_mrsi_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help plot_mrsi

% Last Modified by GUIDE v2.5 17-Aug-2020 12:02:47

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @plot_mrsi_OpeningFcn, ...
                   'gui_OutputFcn',  @plot_mrsi_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


%% --- Executes just before plot_mrsi is made visible.
function plot_mrsi_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to plot_mrsi (see VARARGIN)

na = length(varargin);
if na<2, warning('Insufficient #input arguments: na(=%d)<2',na); end
spec = varargin{1};
if na>1
    par = varargin{2};
else
    par = struct;
end

if length(size(spec))>4
    warning('spec>4D: displaying first timestep only');
    spec = spec(:,:,:,:,1);
end
if length(size(spec))<3, error('spec<3D'); end
[n1,n2,n3,n4] = size(spec);
spec = spec/max(abs(spec(:)));      % normalise max to 1

handles.spec = spec;
if isfield(par,'pc')
    handles.pc = par.pc;
else
    handles.pc = 0;
end
if isfield(par,'lb')
    if isempty(par.lb)
        lb = 0; 
    else
        lb = par.lb(1);
    end
else
    warning('par.lb not existing: setting lb=0');
    lb = 0;
end
handles.lb = lb;
if isfield(par,'verb')
    verb = par.verb;
else
    verb = false;
end
handles.verb = verb;
if isfield(par,'fax')
    if length(par.fax)==n1
        fax = par.fax(:);
    else
        warning('size(fax,1)(=%d)~=size(spec,1)(=%d) -> ignoring fax',...
            size(par.fax,1),n1);
    end
else
    warning('no field par.fax');
end
if ~exist('fax','var')
    warning('fax not existing; creating; assuming bw=5000[Hz]');
    fax = (-n1/2:n1/2-1)/n1*5000;
end
handles.fax = fax;
handles.df = abs(fax(2)-fax(1));
handles.show_line1 = false; % only show lines after moving slider
handles.show_line2 = false;
if isfield(par,'fig_str')
    fig_str = [par.fig_str ': '];
else
    fig_str = '';
end
handles.fig_str = fig_str;


% adjust sliders with Min,Max,Value
sub_adjust_sliders(handles,verb);

handles.output = hObject;   % Choose default command line output for plot_mrsi
guidata(hObject, handles);  % Update handles structure
% uiwait(handles.figure1);  % wait for user response (see UIRESUME)
sub_plot_mrsi(handles,verb); % actual plotting


%% --- Outputs from this function are returned to the command line.
function varargout = plot_mrsi_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



%% create functions
% --- Executes during object creation, after setting all properties.
% hObject    handle to slider_min3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: slider controls usually have a light gray background.
function slider_min1_CreateFcn(hObject, eventdata, handles)
set(hObject,'BackgroundColor',[.9 .9 .9]);

function slider_max1_CreateFcn(hObject, eventdata, handles)
set(hObject,'BackgroundColor',[.9 .9 .9]);

function slider_min2_CreateFcn(hObject, eventdata, handles)
set(hObject,'BackgroundColor',[.9 .9 .9]);

function slider_max2_CreateFcn(hObject, eventdata, handles)
set(hObject,'BackgroundColor',[.9 .9 .9]);

function slider_min3_CreateFcn(hObject, eventdata, handles)
set(hObject,'BackgroundColor',[.9 .9 .9]);

function slider_max3_CreateFcn(hObject, eventdata, handles)
set(hObject,'BackgroundColor',[.9 .9 .9]);

function slider_ind4_CreateFcn(hObject, eventdata, handles)
set(hObject,'BackgroundColor',[.9 .9 .9]);

function slider_scale_CreateFcn(hObject, eventdata, handles)
set(hObject,'BackgroundColor',[.9 .9 .9]);

function slider_line1_CreateFcn(hObject, eventdata, handles)
set(hObject,'BackgroundColor',[.9 .9 .9]);

function slider_line2_CreateFcn(hObject, eventdata, handles)
set(hObject,'BackgroundColor',[.9 .9 .9]);

function slider_pc_CreateFcn(hObject, eventdata, handles)
set(hObject,'BackgroundColor',[.9 .9 .9]);

function slider_lb_CreateFcn(hObject, eventdata, handles)
set(hObject,'BackgroundColor',[.9 .9 .9]);


%% callback functions
% --- Executes on slider movement.
% hObject    handle to slider_min3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
function slider_min1_Callback(hObject, eventdata, handles)
smax = get(handles.slider_max1,'Value');
sub_slider(hObject,1,1,smax);
sub_plot_mrsi(handles);

function slider_max1_Callback(hObject, eventdata, handles)
smin = get(handles.slider_min1,'Value');
sub_slider(hObject,1,smin);
sub_plot_mrsi(handles);

function slider_min2_Callback(hObject, eventdata, handles)
smax = get(handles.slider_max2,'Value');
sub_slider(hObject,1,1,smax);
sub_plot_mrsi(handles);

function slider_max2_Callback(hObject, eventdata, handles)
smin = get(handles.slider_min2,'Value');
sub_slider(hObject,1,smin);
sub_plot_mrsi(handles);

function slider_min3_Callback(hObject, eventdata, handles)
smax = get(handles.slider_max3,'Value');
sub_slider(hObject,1,1,smax);
sub_plot_mrsi(handles);

function slider_max3_Callback(hObject, eventdata, handles)
smin = get(handles.slider_min3,'Value');
sub_slider(hObject,1,smin);
sub_plot_mrsi(handles);

function slider_ind4_Callback(hObject, eventdata, handles)
sub_slider(hObject,1);
sub_plot_mrsi(handles);

function slider_scale_Callback(hObject, eventdata, handles)
sub_slider(hObject,0.1);
sub_plot_mrsi(handles);

function slider_line1_Callback(hObject, eventdata, handles)
sub_slider(hObject,0.01,'','',handles.verb);
handles.show_line1 = true;
sub_plot_mrsi(handles);
guidata(hObject, handles);

function slider_line2_Callback(hObject, eventdata, handles)
sub_slider(hObject,handles.df);
handles.show_line2 = true;
sub_plot_mrsi(handles);
guidata(hObject, handles);

function slider_pc_Callback(hObject, eventdata, handles)
sub_slider(hObject,0.1);
sub_plot_mrsi(handles);

function slider_lb_Callback(hObject, eventdata, handles)
sub_slider(hObject,0.1);
sub_plot_mrsi(handles);

function radiobutton_abs_Callback(hObject, eventdata, handles)
do_abs = get(handles.radiobutton_abs,'Value');
if do_abs
    set(handles.slider_pc,'Visible','off');
else
    set(handles.slider_pc,'Visible','on');    
end
sub_plot_mrsi(handles);

function closebutton_Callback(hObject, eventdata, handles)
close(handles.figure1)
if isdeployed, return; end


%% sub-function for actual plotting
function sub_plot_mrsi(handles,verb)
if ~exist('verb','var'), verb = []; end
if isempty(verb),        verb = false; end

% reading states of sliders
if verb, fprintf('reading states of sliders\n'); end
slice = round(get(handles.slider_ind4,'Value'));
scale = get(handles.slider_scale,'Value');
min1 = get(handles.slider_min1,'Value');
max1 = get(handles.slider_max1,'Value');
min2 = get(handles.slider_min2,'Value');
max2 = get(handles.slider_max2,'Value');
min3 = get(handles.slider_min3,'Value');
max3 = get(handles.slider_max3,'Value');
do_abs = get(handles.radiobutton_abs,'Value');
if min1>max1
    warning('min1(=%g)>max1(=%g): ignoring',min1,max1); 
    min1 = 1; max1 = size(handles.spec,1);
end
if min2>max2
    warning('min2(=%g)>max2(=%g): ignoring',min2,max2); 
    min2 = 2; max2 = size(handles.spec,2);
end
if min3>max3
    warning('min3(=%g)>max3(=%g): ignoring',min3,max3); 
    min3 = 1; max3 = size(handles.spec,3);
end

% adjust spectrum
if verb, fprintf('adjust spectrum\n'); end
spec = handles.spec(min1:max1,min2:max2,min3:max3,slice);
[n1,n2,n3] = size(spec);
% line-broadening
lb = get(handles.slider_lb,'Value')-handles.lb;
if abs(lb)>0.01
    if verb>1
        fprintf('lb = %g\n',lb); 
    end
    fid = ifft(fftshift(spec,1),[],1);
  
    fax = handles.fax(min1:max1);
    bw = max(fax)-min(fax);
    lbf = lb_fun(size(spec,1),bw,[0 lb]).';
    if lb<0, lbf = 1./lbf; end
    sca = sqrt((lb+handles.lb+2*pi)/(handles.lb+2*pi));
    
    spec = sca*ifftshift(fft(bsxfun(@times,fid,lbf),[],1),1);
    % spec = ifft(bsxfun(@times,fid,lbf),[],1);
end

% phase correction
if ~isreal(spec)
    if do_abs
        spec = abs(spec);
    else
        pc = get(handles.slider_pc,'Value');
        spec = real(spec*exp(1i*pc*pi/180));
    end
end
spec = scale*spec;
spec(spec>1) = 1;

% concatenate and plot spectra
if verb, fprintf('concatenate and plot spectra\n'); end
for l2=1:n2
    plot(linspace(-n3/2,n3/2,n1*n3),...
        reshape(spec(:,l2,:)-l2+n2/2,[n1*n3 1]),'b');
    if l2==1, hold on; end
end

% plot MRS grid
for l2=(-n2/2:n2/2), plot([-1 1]*n3/2,l2*[1 1],'k'); end
for l3=(-n3/2:n3/2), plot(l3*[1 1],[-1 1]*n2/2,'k'); end

% plot lines
if handles.show_line1
    fax = handles.fax;
    line1 = (get(handles.slider_line1,'Value')-fax(min1))/(fax(max1)-fax(min1));
    for l3=(-n3/2:n3/2-1)
        plot((line1+l3)*[1 1],[-1 1]*n2/2,'r:'); 
    end
end
if handles.show_line2
    fax = handles.fax;
    line2 = (get(handles.slider_line2,'Value')-fax(min1))/(fax(max1)-fax(min1));
    for l3=(-n3/2:n3/2-1)
        plot((line2+l3)*[1 1],[-1 1]*n2/2,'g:'); 
    end
end
hold off
bdr = 0.02;
axis([-n3/2-bdr n3/2+bdr -n2/2-bdr n2/2+bdr]);

if verb, fprintf('set XTick,YTick,Name\n'); end
set(gca,'XTick',[])
set(gca,'YTick',[])
fig_str = sprintf('plot_mrsi: %sspec(%d:%d,%d:%d,%d:%d,%d)',...
    handles.fig_str,min1,max1,min2,max2,min3,max3,slice);
set(gcf,'Name',fig_str);
if verb, fprintf('done sub_plot_mrsi\n'); end


%% read and round slider
function val = sub_slider(hdl,prec,vmin,vmax,verb)
if ~exist('prec','var'), prec = []; end
if ~exist('vmin','var'), vmin = []; end
if ~exist('vmax','var'), vmax = []; end
if ~exist('verb','var'), verb = []; end
if isempty(verb), verb = false; end

val = get(hdl,'Value');
if ~isempty(prec), val = round(val/prec)*prec; end
if ~isempty(vmin)
    if val<vmin
        if verb
            fprintf('val(=%g)<vmin(=%g): setting val=vmin\n',val,vmin); 
        end
        val = vmin; 
    end
end
if ~isempty(vmax)
    if val>vmax
        if verb
            fprintf('val(=%g)>vmax(=%g): setting val=vmax\n',val,vmax);
        end
        val = vmax; 
    end
end
if verb, fprintf('%s=%g\n',hdl.Tag,val); end
set(hdl,'Value',val);
set(hdl,'TooltipString',num2str(val));


%% adjust sliders
function sub_adjust_sliders(handles,verb)

[n1,n2,n3,n4] = size(handles.spec);
if verb, fprintf('Setting slider_min1 & max1\n'); end
set(handles.slider_min1,'Max',n1);
set(handles.slider_min1,'SliderStep',[1 10]/(n1-1));
if get(handles.slider_max1,'Max')~=n1
    set(handles.slider_max1,'Max',n1);
    set(handles.slider_max1,'Value',n1);
end
set(handles.slider_max1,'TooltipString',num2str(n1));
set(handles.slider_max1,'SliderStep',[1 10]/(n1-1));

if verb, fprintf('Setting slider_min2 & max2\n'); end
set(handles.slider_min2,'Max',n2);
set(handles.slider_min2,'SliderStep',[1 1]/(n2-1));
if get(handles.slider_max2,'Max')~=n2
    set(handles.slider_max2,'Max',n2);
    set(handles.slider_max2,'Value',n2);
end
set(handles.slider_max2,'TooltipString',num2str(n2));
set(handles.slider_max2,'SliderStep',[1 1]/(n2-1));

if verb, fprintf('Setting slider_min3 & max3\n'); end
set(handles.slider_min3,'Max',n3);
set(handles.slider_min3,'SliderStep',[1 1]/(n3-1));
if get(handles.slider_max3,'Max')~=n3
    set(handles.slider_max3,'Max',n3);
    set(handles.slider_max3,'Value',n3);
end
set(handles.slider_max3,'TooltipString',num2str(n3));
set(handles.slider_max3,'SliderStep',[1 1]/(n3-1));

if verb, fprintf('Setting slider_ind4'); end
if n4==1
    set(handles.slider_ind4,'Visible','off');
    set(handles.text_ind4,'String','2D');
else
    if get(handles.slider_ind4,'Max')~=n4
        set(handles.slider_ind4,'Max',n4);
        set(handles.slider_ind4,'Value',round(n4/2));
    end
    set(handles.slider_ind4,'TooltipString',num2str(round(n4/2)));
    set(handles.slider_ind4,'SliderStep',[1 1]/(n4-1));
end

if verb, fprintf('Setting slider_line1 & line2\n'); end
fmin = min(handles.fax);
fmax = max(handles.fax);
fmean = (fmin+fmax)/2;
set(handles.slider_line1,'Min',fmin);
set(handles.slider_line1,'Max',fmax);
set(handles.slider_line1,'Value',fmean);
set(handles.slider_line1,'SliderStep',[1 10]/(n1-1));

set(handles.slider_line2,'Min',fmin);
set(handles.slider_line2,'Max',fmax);
set(handles.slider_line2,'Value',fmean);
set(handles.slider_line2,'SliderStep',[1 10]/(n1-1));

if verb, fprintf('Setting slider_pc\n'); end
if isreal(handles.spec)
    set(handles.slider_pc,'Visible','off');
    set(handles.text_pc,'String','spec is real');
    set(handles.radiobutton_abs,'Visible','off');
else
    if get(handles.radiobutton_abs,'Value')
        set(handles.slider_pc,'Visible','off');
    else
        set(handles.slider_pc,'Visible','on');
    end
end
set(handles.slider_lb,'Value',handles.lb);
set(handles.slider_lb,'TooltipString',num2str(handles.lb));

if verb, fprintf('Done with sub_adjust_sliders\n'); end
