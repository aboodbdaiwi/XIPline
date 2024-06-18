function varargout = plot_spec(varargin)
% PLOT_SPEC
%    plot_spec(spec,par)
%      spec  spectrum (f-domain)
%       par  structure
%            f0 = synthesizer_frequency;
%            bw = sample_frequency;
%
%See also ECHO2SPEC.
% 1/2008 Rolf Schulte

% Last Modified by GUIDE v2.5 26-Jul-2006 11:03:05
if nargin<1, help(mfilename); return; end
% warning('f-axis reversed');

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @plot_spec_OpeningFcn, ...
                   'gui_OutputFcn',  @plot_spec_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
                   'gui_Callback',   []);
if nargin && isstr(varargin{1})
  gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
  [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
  gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% ----------------------------
% Executes just before plot_spec is made visible.
function plot_spec_OpeningFcn(hObject, eventdata, handles, varargin)
na = length(varargin);
handles.output = hObject;
spec = varargin{1};

if length(size(spec))>2, error('spec > 2D'); end
par = varargin{2};
f0 = par.synthesizer_frequency;
bw = par.sample_frequency;

si = size(spec);
handles.spec = spec;
if f0>1d3, f0 = f0*1d-6; end
if bw>100, bw = bw/f0; end
fprintf('f0 = %g [MHz]\nbw = %g [ppm]\n',f0,bw);
handles.f0 = f0;
handles.ppm = [-si(2)/2:si(2)/2-1]*bw/si(2);
% handles.ppm = -[-si(2)/2:si(2)/2-1]*bw/si(2);
% fprintf('Attention:\n');
% fprintf('Putting minus in front of ppm axis; if reversed->remove\n');
% same effect as fid=conj(fid) in echo2spec

switch round(f0)
    case 128, zp = 4.65;
    case 32,  zp = 183.2;
    otherwise, zp = 0;
end
fprintf('zp  = %g\n',zp);
set(handles.zp_txt,'String',num2str(zp));
set_slider(handles.zp,zp);
set(handles.line_txt,'String',num2str(zp));
set_slider(handles.line,zp);

handles.figid = gcf;
figure(gcf); clf;
zoom off;
set(handles.figure1,'Name','plot_spec');	% handle of gui
guidata(hObject, handles);
plot_spec_plot(handles);


% ----------------------------
% parse output
function varargout = plot_spec_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;


% ----------------------------
% create functions
% ----------------------------

function phi_CreateFcn(hObject, eventdata, handles)
set(hObject,'BackgroundColor',[.9 .9 .9]);

function phi_txt_CreateFcn(hObject, eventdata, handles)
set(hObject,'BackgroundColor','white');

function phi1_CreateFcn(hObject, eventdata, handles)
set(hObject,'BackgroundColor',[.9 .9 .9]);

function phi1_txt_CreateFcn(hObject, eventdata, handles)
set(hObject,'BackgroundColor','white');

function turn_CreateFcn(hObject, eventdata, handles)
set(hObject,'BackgroundColor',[.9 .9 .9]);

function turn_txt_CreateFcn(hObject, eventdata, handles)
set(hObject,'BackgroundColor','white');

function zp_CreateFcn(hObject, eventdata, handles)
set(hObject,'BackgroundColor',[.9 .9 .9]);

function zp_txt_CreateFcn(hObject, eventdata, handles)
set(hObject,'BackgroundColor','white');

function line_CreateFcn(hObject, eventdata, handles)
set(hObject,'BackgroundColor',[.9 .9 .9]);

function line_txt_CreateFcn(hObject, eventdata, handles)
set(hObject,'BackgroundColor','white');

function axis_unit_CreateFcn(hObject, eventdata, handles)
set(hObject,'BackgroundColor','white');

function baseline_CreateFcn(hObject, eventdata, handles)
set(hObject,'BackgroundColor','white');


% ----------------------------
% callback functions
% ----------------------------

function phi_Callback(hObject, eventdata, handles)
val = read_slider(hObject,0);
set(handles.phi_txt,'String',num2str(val));
plot_spec_plot(handles);

function phi_txt_Callback(hObject, eventdata, handles)
val = read_editable(hObject);
if isempty(val), val = 0; end
set_slider(handles.phi,val);
plot_spec_plot(handles);

function phi1_Callback(hObject, eventdata, handles)
val = read_slider(hObject,-2);
set(handles.phi1_txt,'String',num2str(val));
plot_spec_plot(handles);

function phi1_txt_Callback(hObject, eventdata, handles)
val = read_editable(hObject);
set_slider(handles.phi1,val);
plot_spec_plot(handles);

function turn_Callback(hObject, eventdata, handles)
val = read_slider(hObject,-2);
set(handles.turn_txt,'String',num2str(val));
plot_spec_plot(handles);

function turn_txt_Callback(hObject, eventdata, handles)
val = read_editable(hObject);
set_slider(handles.turn,val);
plot_spec_plot(handles);

function zp_Callback(hObject, eventdata, handles)
val = read_slider(hObject,-2);
set(handles.zp_txt,'String',num2str(val));
plot_spec_plot(handles);

function zp_txt_Callback(hObject, eventdata, handles)
val = read_editable(hObject);
set_slider(handles.zp,val);
plot_spec_plot(handles);

function line_Callback(hObject, eventdata, handles)
val = read_slider(hObject,-2);
set(handles.line_txt,'String',num2str(val));
plot_spec_plot(handles);

function line_txt_Callback(hObject, eventdata, handles)
val = read_editable(hObject);
set_slider(handles.line,val);
plot_spec_plot(handles);

function axis_unit_Callback(hObject, eventdata, handles)
apd = getappdata(handles.figid);
if isfield(apd,'ZoomOnState'),
    zoom off
end
plot_spec_plot(handles);

function baseline_Callback(hObject, eventdata, handles)
plot_spec_plot(handles);

function file_Callback(hObject, eventdata, handles)

function save_Callback(hObject, eventdata, handles)
par = handles.par;
spec = plot_spec_plot(handles);		% recon (phase) spectrum
fid = ifft(ifftshift(spec,2),[],2);	% recon fid
% move zp in t-domain
zp0 = 4.65;				% reference of LCModel
zp = read_editable(handles.zp_txt);
zp = (zp-zp0)*handles.f0;
t = [0:size(spec,2)-1]/par.sample_frequency;
fid = fid.*exp(i*2*pi*t*zp);
% figure(10); plot(t,real(fid));
% figure(11); plot(handles.ppm+zp0,real(fftshift(fft(fid,[],2),2)));
% set(gca,'XDir','reverse');
par.rows = size(fid,1);
par.samples = size(fid,2);
write_sdat([],fid,par);


function quit_Callback(hObject, eventdata, handles)
close(handles.figure1);
if isfield(handles,'figid'),
  figure(handles.figid);
  close(handles.figid);
end


% *******************************************************
% read in double from slider and round
function val = read_slider(hdl,prec)
val = get(hdl,'Value');
val = floor(val/10^prec+0.5)*10^prec;




