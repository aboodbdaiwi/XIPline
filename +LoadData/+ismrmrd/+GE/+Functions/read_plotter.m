function [grad,rf,data] = read_plotter(fname)
%READ_PLOTTER  Read xlm output (CaptureWaveform.xml) of plotter
% [grad,rf,data] = read_plotter(fname)
%
% 7/2017  Rolf Schulte
if (nargin<1), help(mfilename); return; end;

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
    
    if ~ischar(tline),
        disp(tline)
        fprintf('l=%g\n',l);
        warning('~ischar(tline) -> break');
        break
    end
    if isempty(tline), 
        disp(tline)
        fprintf('l=%g\n',l);
        warning('isempty(tline) -> break');
        break
    end
    
    if ~isempty(regexp(tline,'</seq','match')),
        disp(tline)
        if ld==0, warning('data reading mode already off'); end
        ld = 0;
    end

    if ld>0
        fld = sscanf(tline,'%f',4);
        if ~isempty(fld), 
            data{lseq}.fld(ld,:) =  fld;
            ld = ld+1;
        else
            disp(tline);
        end
    end
    if ~isempty(regexp(tline,'<seq','match')),
        disp(tline)
        lseq = lseq+1;
        data{lseq}.desc = tline;
        if ld>0, warning('already in data reading mode'); end
        ld = 1;
    end
    
end
fclose(fid);

if l==nmax, 
    warning('reaching nmax: l(=%g)==nmax(=%g); increase nmax',l,nmax); 
end

%% converting corner points into discretised vector
tmax = 0;
for l=1:length(data)
    if regexpi(data{l}.desc,'mgd::x'), ix = l; end
    if regexpi(data{l}.desc,'mgd::y'), iy = l; end
    if regexpi(data{l}.desc,'mgd::z'), iz = l; end
    if regexpi(data{l}.desc,'mgd::rho1'), ir = l; end
    if data{l}.fld(end,1)>tmax, tmax = data{l}.fld(end,1); end
end
tmax = ceil(tmax/dt)+1;
tax = (0:dt:tmax-1);
grad = zeros(tmax,4);
grad(:,1) = sub_digitise_corner_points(data{ix}.fld,dt,tmax);
grad(:,2) = sub_digitise_corner_points(data{iy}.fld,dt,tmax);
grad(:,3) = sub_digitise_corner_points(data{iz}.fld,dt,tmax);
grad(:,4) = linspace(0,tmax/dt-1,tmax)*1d-6;
rf = sub_digitise_corner_points(data{ir}.fld,dt,tmax);


end   % end main function read_plotter.m


%% subfunction digitise_corner_points
function dig=sub_digitise_corner_points(cp,dt,np)

%% convert corner points to digitised gradient trajectory
% np = ceil(cp(1,end)/dt)+1;
dig = zeros(np,1);
if cp(1,1)~=0, warning('time not starting at zero: cp(1,1)=%g',cp(1,1)); end
for l=2:size(cp,1),
    n1 = round(cp(l-1,1)/dt)+1;
    n2 = cp(l,1)/dt+1;
    if n2<0, error('n2(=%g) < 0',n2); end
    if abs(n2-floor(n2)) > 1d-10, 
        warning('n2(=%g) not natural number; rounding',n2);
        n2 = round(n2);
    end
    dig(n1:n2,1) = linspace(cp(l-1,2),cp(l,2),n2-n1+1);
end

if size(dig,1)~=np, warning('size(dig,1)(=%g)~=np(=%g)',size(grad,1),np); end

end   % end subfunction digitise_corner_points

