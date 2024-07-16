function pc = ssfp_phase_cycling(nexc,npc,sv,docen)
%SSFP_PHASE_CYCLING Generate phase cycling scheme for SSFP
%   pc = ssfp_phase_cycling(nexc,npc,sv,docen)
%                                                          (default)
% nexc   #views for one full ssfp encoding
%  npc   #ssfp phase cycles                                (4)
%   sv   save matrix to fdl file for fidall                (false)
%docen   start with centred band (0-180-0-...)             (true)
%
%   pc   phase cycling matrix
%  
% See also write_fdl.m
% 12/2020 Rolf Schulte
if (nargin<1), help(mfilename); return; end


%% default parameters + checks
if ~exist('npc','var'), npc = []; end
if isempty(npc),        npc = 4; end
if ~exist('sv','var'),  sv = []; end
if isempty(sv),         sv = false; end
if ~exist('docen','var'), docen = []; end
if isempty(docen),      docen = true; end
if abs((nexc/npc-floor(nexc/npc)))>1d-10
    error('nexc(=%g) must be multiple of npc(=%g)',nexc,npc);
end


%% calculate phase cycling scheme
dp = 360*(0:(npc-1))/npc;    % SSFP phase advance (normally 180°)
if docen, dp = mod(dp+180,360); end
xx = zeros(npc,npc);         % single phase cycling matrix
for l=1:npc
    xx(l,:) = mod((0:(npc-1))*dp(l),360);
end


%% repeat scheme for fidall
pc = []; 
for l=1:npc
    pc = [pc , repmat(xx(l,:),[1 nexc/npc])]; 
end


%% print info
fprintf('#phase cycles = %g\n',npc);
fprintf('phase advance = %g [deg]\n',dp);
fprintf('whole phase cycling matrix\n');
disp(xx);


%% save into fdl file
if sv
    fname = sprintf('vap_phase_ssfp_npc%g_nexc%g',npc,nexc);
    fprintf('saving file to ''%s''\n',fname);
    write_fdl(fname,pc,'phase');
end


end      % main function ssfp_phase_cycling.m
