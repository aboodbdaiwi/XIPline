function time = read_fidall_gating(fname)
%READ_FIDALL_GATING Read in gating times recorded by fidall -> rfs_recgat
%
% time = read_fidall_gating(fname)
%  fname  filename  [string]
%   time  recorded gating times [1,#times]
%
% 26.07.2012 Rolf Schulte
if (nargin<1), help(mfilename); return; end;

if ~exist(fname,'file'), error('File %s not found',fname); end 
fid = fopen(fname,'r'); 
tmp = fread(fid,'*char')'; 
i1 = find(tmp=='=');
i2 = find(tmp=='('); i2=[i2,length(tmp)];
for ii=1:length(i1)
    time(ii)=str2num(tmp(i1(ii)+2:i2(ii)-3));
end
fclose(fid);
time = time(2:end)-time(2);
