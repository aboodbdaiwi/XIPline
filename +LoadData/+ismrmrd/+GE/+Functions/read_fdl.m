function [x,npars,scaling]=read_fdl(fname)
%READ_FDL  Read vector x readable for fidall
%[x,npars,scaling]=read_fdl(fname)
%  fname  filename
%      x  binary vector
%  npars  length of binary vector
%scaling  scaling factor (already applied to x)
%
% 02/2018 Rolf Schulte
if (nargin<1), help(mfilename); return; end;

if ~ischar(fname), error('fname not char'); end
if isempty(regexp(fname,'\.', 'once')), fname = [fname '.fdl']; end


%% actual writing
if ~exist(fname,'file'), error('fname(=%s) not found',fname); end
fid =     fopen(fname,'rb');
npars =   fread(fid,1,'int32');
scaling = fread(fid,1,'double');
spare =   fread(fid,1,'int32');   % spare parameter
x =       fread(fid,inf,'int32');
fclose(fid);
x = x/scaling;


%% checks
if length(x)~=npars, warning('length(x)(=%g)~=npars(=%g)',length(x),npars); end
if ~isnumeric(x),  warning('x not numeric'); end
if ~isvector(x),   warning('x not a vector'); end
if any(x<0),       warning('x<0'); end

if npars>16382, warning('npars(=%g) exceeding fidall limits',npars); end
if npars>2^31-1, error('npars(=%g) too long',npars); end
if any(abs(x))> 2^31-1
    error('max(abs(x)))(=%g) exceeding 32bits',max(abs(x))); 
end

