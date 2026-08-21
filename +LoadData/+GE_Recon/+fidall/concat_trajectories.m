function concat_trajectories(wfns,cci_wfname,recons)
%CONCAT_TRAJECTORIES
% concat_trajectories(wfns,cci_wfname,recons)
%      wfns  Input waveform files (cell array)
%cci_wfname  Output with concatenated trajectories (.wav) and info (.mat)
%            (string)
%    recons  Specify which recons to use in recon_concat               ([])
%            (cell array; default from trajectory)
%
%  7/2024 Rolf Schulte
if (nargin<1), help(mfilename); return; end


%% input
if ~exist('recons','var'), recons = []; end
if ~isempty(recons) && ~iscell(recons), error('recons must be cell array'); end


%% checking waveform files/names
nwf = length(wfns);
suf = {'.mat','.wav'};
if ~iscell(wfns), error('wfns must be cell with multiple file names'); end
if nwf<2, error('length(wfns)(=%d)<2',nwf); end
for l1=1:nwf
    if ~isempty(regexpi(wfns{l1},'\.wav$', 'once')) || ...
            ~isempty(regexpi(wfns{l1},'\.mat$', 'once'))
        wfns{l1} = wfns{l1}(1:end-4);
    end
    if ~isempty(regexpi(wfns{l1},'\\', 'once')) || ...
            ~isempty(regexpi(wfns{l1},'/', 'once'))
        warning('Use file name only; omit dir');
    end
    for l2=1:2
        if ~exist([wfns{l1} suf{l2}],'file')
            error('file not found: wfns{%d}=%s',l1,[wfns{l1} suf{l2}]);
        end
    end
end
if  iscell(cci_wfname), error('wfn_out must not be cell'); end
if ~isempty(regexpi(cci_wfname,'\.wav$', 'once')) || ...
        ~isempty(regexpi(cci_wfname,'\.mat$', 'once'))
    cci_wfname = cci_wfname(1:end-4); 
end
if exist([cci_wfname '.wav'],'file') || exist([cci_wfname '.mat'],'file') 
    warning('file existing: wfn_out=%s; overwriting',cci_wfname); 
end


%% reading grad wav and check compatibility
for l1=1:nwf
    fprintf('Reading: ''%s.wav''\n',wfns{l1});
    [grad{l1},bw{l1},fov{l1},desc{l1},N{l1},params{l1}] = ...
        read_ak_wav([wfns{l1} '.wav']);
end
for l1=2:nwf
    if (bw{l1} ~= bw{1}),   error('bw{l1}(=%g) ~= bw{1}(=%g)',bw{l1},bw{1}); end
    if (fov{l1} ~= fov{1}), error('fov{l1}(=%g) ~= fov{1}(=%g)',fov{l1},fov{1}); end
end


%% concatenate and save grad
ngpt = 0; nexc = 0; ndim = 0; nkpt = 0;
for l1=1:nwf
    if size(grad{l1},1)>ngpt, ngpt = size(grad{l1},1); end
    if size(grad{l1},3)>ndim, ndim = size(grad{l1},3); end
    nexc = nexc + size(grad{l1},2);
    if params{l1}(7)>nkpt, nkpt = params{l1}(7); end
end
fprintf('ngpt = %d, nexc = %d, ndim = %d\n',ngpt,nexc,ndim);
desc = 'concatenated trajectories';
grad_out = zeros(ngpt,nexc,ndim);
lexc = 0;
for l1=1:nwf
    for l2=1:3, ii{l2} = 1:size(grad{l1},l2); end
    grad_out(ii{1},ii{2}+lexc,ii{3}) = grad{l1};
    ind_exc{l1} = false(1,nexc);
    ind_exc{l1}(1,ii{2}+lexc) = true;
    lexc = lexc + size(grad{l1},2);
end
fprintf('Writing: ''%s.wav''\n',cci_wfname);
write_ak_wav([cci_wfname '.wav'],grad_out,bw{1},fov{1},desc,nkpt,[],[],4);


%% write info for recon_concat.m
% save('tmp.mat','-struct','wf');
if isempty(recons)
    for l1=1:nwf
        reco = [];
        if ~isempty(regexpi(wfns{l1},'spiral', 'once')), reco = 'spiral'; end
        if ~isempty(regexpi(wfns{l1},'epi', 'once')), reco = 'epi'; end
        if isempty(reco)
            error('no reconstruction found for wfn=''%s''',wfns{l1});
        end
        recons{l1} = reco;
    end
end
if length(recons)~=nwf
    warning('length(recons)(=%d) ~= nwf(=%d)',length(recons),nwf);
end
save([cci_wfname '.mat'],'wfns','ind_exc','recons');


end      % concat_trajectories.m
