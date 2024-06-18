function write_fdl(fname,x,what)
%WRITE_FDL  Write vector x into binary file readable for fidall
%  write_fdl(fname,x,what)
%  fname  output filename
%         use meaningful name and link it to 'vap_XXX<cv20>.fdl')
%      x  binary vector
%   what  'te','tr'                     [us]
%         'phase','flip','phi','theta'  [deg]
%         'freq'                        [Hz]
%         'dudri'        [0=I+Q on, 1=Ion+Qoff, 2=Ioff+Qon,3=I+Q off]
%         scaling factor
%
% 7/2019 Rolf Schulte
if (nargin<1), help(mfilename); return; end

if ~ischar(fname), error('fname not char'); end
if isempty(regexp(fname,'\.', 'once')), fname = [fname '.fdl']; end
if ~isnumeric(x),  error('x not numeric'); end
if ~isvector(x),   error('x not a vector'); end

if isnumeric(what)
    if length(what)~=1, warning('length(what)~=1'); end
    scaling = what(1);
    what = 'scaling';
end

min_val = 0;
switch lower(what(1:2))
    case 'te', scaling = 1;     max_val = 2^31-1;
    case 'tr', scaling = 1;     max_val = 2^31-1;
    case 'ph', scaling = 10000; max_val = 360;
    case 'th', scaling = 10000; max_val = 360;
    case 'fl', scaling = 10000; max_val = 180;
    case 'du', scaling = 1;     max_val = 3;
    case 'sc', min_val = -inf;  max_val = inf;
    case 'fr', scaling = 100;   min_val = -1d4; max_val = 1d4;
    otherwise, error('what (=%s) unknown',what);
end
if any(x>max_val), error('x>%g',max_val); end
if any(x<min_val), error('x<%g',min_val); end
if (abs(scaling-1)>eps), x = x*scaling; end
npars = length(x);
if npars>16382*16
    warning('npars(=%g) exceeding fidall limits',npars); 
end
if npars>2^31-1, error('npars(=%g) too long',npars); end
if any(abs(x))> 2^31-1
    error('max(abs(x)))(=%g) exceeding 32bits',max(abs(x))); 
end

%% actual writing
fid = fopen(fname,'wb');
fwrite(fid,int32(npars),   'int32');
fwrite(fid,double(scaling),'double');
fwrite(fid,int32(0),       'int32');   % spare parameter
fwrite(fid,int32(x),       'int32');
fclose(fid);

end      % main function write_fdl.m
