function [grad,rf,ssp,tt,data] = read_plotter(fname,gdt)
%READ_PLOTTER  Read xlm output (CaptureWaveform.xml) of plotter
% [grad,rf,tax,data] = read_plotter(fname)
%  fname  Filename of plotter output 'CaptureWaveform.xml.xxx'
%    gdt  Interpolate down from 1us to gdt            [us]      (4)
%   grad  Gradient axes                               [T/m]
%     rf  RF waveform
%     tt  Time axis                                   [s]
%   data  Corner points
%
%  9/2021 Rolf Schulte
if (nargin<1), help(mfilename); return; end

if ~exist('gdt','var'), gdt = []; end
if isempty(gdt),        gdt = 4; end
nmax = 5d5;       % maximum #points to read in
dt = 1;           % time discretisation


%% load xml file 
fprintf('Loading xml file ''%s''\n',fname);
if ~exist(fname,'file'), error('file %s not found',fname); end
fid = fopen(fname,'r');


%% parsing in input
lseq = 0;
ld = 0;
for l=1:nmax
    tline = fgetl(fid);
    
    if isempty(tline)
        disp(tline)
        fprintf('l=%g\n',l);
        warning('isempty(tline) -> break');
        break
    end
    if tline(1)==-1, break; end 
    if ~ischar(tline)
        disp(tline)
        fprintf('l=%g\n',l);
        warning('~ischar(tline) -> break');
        break
    end
    
    if ~isempty(regexp(tline,'</seq','match'))
        disp(tline)
        if ld==0, warning('data reading mode already off'); end
        ld = 0;
    end

    if ld>0
        fld = sscanf(tline,'%f',4);
        if ~isempty(fld)
            data{lseq}.fld(ld,:) =  fld;
            ld = ld+1;
        else
            disp(tline);
        end
    end
    if ~isempty(regexp(tline,'<seq','match'))
        disp(tline)
        lseq = lseq+1;
        data{lseq}.desc = tline;
        if ld>0, warning('already in data reading mode'); end
        ld = 1;
    end
    
end
fclose(fid);

if l==nmax
    warning('reaching nmax: l(=%g)==nmax(=%g); increase nmax',l,nmax); 
end


%% converting corner points into discretised vector
tmax = 0;
for l=1:length(data)
    if ~isempty(regexpi(data{l}.desc,'mgd::x')) || ...
            ~isempty(regexpi(data{l}.desc,'gradient x')), ix = l; end
    if ~isempty(regexpi(data{l}.desc,'mgd::y')) || ...
            ~isempty(regexpi(data{l}.desc,'gradient y')), iy = l; end
    if ~isempty(regexpi(data{l}.desc,'mgd::z')) || ...
            ~isempty(regexpi(data{l}.desc,'gradient z')), iz = l; end
    if ~isempty(regexpi(data{l}.desc,'rho')), ir = l; end
    if ~isempty(regexpi(data{l}.desc,'SSP\>')), is = l; end
    iend = size(data{l}.fld,1);
    while iend>1
        if abs(data{l}.fld(iend,2))>1d-10
            break;
        else
            iend = iend-1;
        end
    end
    if data{l}.fld(iend,1)>tmax
        tmax = data{l}.fld(iend,1); 
    end
end

tmax = ceil(tmax/dt)+1;
grad = zeros(tmax,3);
if exist('ix','var')
    grad(:,1) = sub_digitise_corner_points(data{ix}.fld,dt,tmax);
else
    warning('no x found'); 
end
if exist('iy','var')
    grad(:,2) = sub_digitise_corner_points(data{iy}.fld,dt,tmax);
else
    warning('no y found');
end
if exist('iz','var')
    grad(:,3) = sub_digitise_corner_points(data{iz}.fld,dt,tmax);
else
    warning('no z found'); 
end
% grad(:,4) = linspace(0,tmax/dt-1,tmax)*1d-6;
if exist('ir','var')
    rf = sub_digitise_corner_points(data{ir}.fld,dt,tmax);
else
    rf = [];
    warning('no rho found'); 
end

if exist('is','var')
    ssp = sub_digitise_corner_points(data{is}.fld,dt,tmax);
else
    ssp = [];
    warning('no ssp found'); 
end

grad = grad/100;   % convert to SI units: G/cm -> T/m 
% if plotter output stored in physical units;

if abs(gdt-dt)>1d-10
    grad = interp1((0:dt:tmax-1),grad,(0:gdt:tmax-1));
    rf = interp1((0:dt:tmax-1),rf,(0:gdt:tmax-1));
    ssp = interp1((0:dt:tmax-1),ssp,(0:gdt:tmax-1));
end

tt = (0:gdt:tmax-1)*1d-6;


end      % read_plotter.m


%% subfunction digitise_corner_points
function dig=sub_digitise_corner_points(cp,dt,np)

%% convert corner points to digitised gradient trajectory
% np = ceil(cp(1,end)/dt)+1;
dig = zeros(np,1);
if cp(1,1)~=0, warning('time not starting at zero: cp(1,1)=%g',cp(1,1)); end
for l=2:size(cp,1)
    n1 = round(cp(l-1,1)/dt)+1;
    n2 = cp(l,1)/dt+1;
    if n2<0, error('n2(=%g) < 0',n2); end
    if abs(n2-floor(n2)) > 1d-10
        warning('n2(=%g) not natural number; rounding',n2);
        n2 = round(n2);
    end
    dig(n1:n2,1) = linspace(cp(l-1,2),cp(l,2),n2-n1+1);
end

if size(dig,1)~=np
    fprintf('Attention: size(dig,1)(=%g)~=np(=%g); cropping\n',size(dig,1),np); 
    dig = dig(1:np,1);
end

end      % sub_digitise_corner_points
