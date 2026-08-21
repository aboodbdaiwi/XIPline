function [d,h,spec] = psf_mrsi(wfn,verb)
%PSF_MRSI  Simulate Point-Spread-Function for arbitrarily sampled MRSI
%[d,h,spec] = psf_mrsi(wfn,verb)
%                                                                 (default)
%   wfn   Location of waveform .mat file
%  verb   Verbose mode                                              (false)
%     d   Raw p-file data [(nx*ny*nz),ns,1,nt,1,nc] (ns=spectral pts;nc=coils)
%     h   Header structure
%  spec   Simulated spectrum
%
% 8/2023 Rolf Schulte
if nargin<1, help(mfilename); return; end

timerVal = tic;               % record execution time


%% misc input
if ~exist('verb','var'), verb = []; end
if isempty(verb), verb = 0; end


%% reading waveform file
if isempty(wfn), error('wfn empty'); end
if isstruct(wfn)
    wf = wfn;
else
    if ~isempty(regexpi(wfn,'\.wav$', 'once')), wfn = wfn(1:end-4); end
    if isempty(regexpi(wfn,'\.mat$', 'once')),  wfn = [wfn '.mat']; end
    if ~exist(wfn,'file'), error('file not found: wfn=%s',wfn); end
    wf = load(wfn);          % load waveform
end


%% miscelaneous parameters
dim = wf.dim;
mtx = wf.mtx(:).'; mtx = [mtx ones(1,3-length(mtx))];
nx = mtx(1);
ny = mtx(2);
nz = mtx(3);
% ns = wf.npts;                               % number of sampling points
ns = 1; 
bw = wf.bw;   
fov = wf.fov;
nexc = wf.nexc;
if ns>1
    n_ind = wf.n_ind;
    ind = wf.ind;
else
    n_ind = 1;
    ind = true;
end



%% generate header
h.rdb_hdr.spectral_width = bw;    % sampling bandwidth [Hz]
h.rdb_hdr.nslices = 1;            % number of slices
h.image.numecho = 1;              % number of time steps
h.rdb_hdr.dab(1:2) = 0;           % number of Rx coil elements
h.rdb_hdr.user19 = 1;             % #reference frames (w/o WS)
h.rdb_hdr.user23 = 0;
h.rdb_hdr.user26 = 0;             % xloc
h.rdb_hdr.user27 = 0;             % yloc
h.rdb_hdr.user28 = 0;             % zloc
h.rdb_hdr.fov = fov;              % field-of-view
h.image.psd_iname = 'psf_mrsi.m'; 
h.rdb_hdr.user2 = wf.nucleus;
h.rdb_hdr.data_collect_type = 1;  % no chopping
h.data_acq_tab.rotate = 0;
h.data_acq_tab.transpose = 0;
h.image.rawrunnum = 0;


%% generate Dirac image: central voxel=1; rest=0
spec = complex(zeros([ns,mtx]));
spec(:,ceil(nx/2),ceil(ny/2),ceil(nz/2)) = 1;
% spec(:,ceil(nx*3/4),ceil(ny*3/4),ceil(nz*3/4)) = 1;


%% spectral apodisation
if ns>1
    lbf = lb_fun(ns,ns,[0 2]).';
    spec = bsxfun(@times,spec,lbf);
end

%%  inverse spatial reconstruction
if verb>0, fprintf('Inverse spatial recon\n'); end
dd = complex(zeros(ns,nexc));
ss = reshape(spec,[ns,prod(mtx)]);
kk = permute(-1i*2*pi*wf.k,[3 2 1]);
if isfield(wf,'xx')
    xx = permute(wf.xx,[3 2 1]);
    xx = xx+1;
else
    xx = sub_grid(mtx,dim);
end

% if verb>1, wbid = waitbar(0,'inverse gridding'); end
for l1=1:prod(mtx)
    if any(abs(ss(:,l1))>eps)
        dd = dd + bsxfun(@times,ss(:,l1),...
            exp(sum(bsxfun(@times,kk,xx(:,l1)),1)));
    end
    % if verb>1, waitbar(l1/prod(mtx),wbid); end
end
dd = dd/nexc;
if verb>0, fprintf('\n'); end
% if verb>1, close(wbid); end

%% reshaping data and indexing
d = zeros(nexc,n_ind);
d(:,ind) = dd.';


%%
if verb>0
    fprintf('psf_mrsi.m finished: ');
    fprintf('total simulation time = %.2f s\n',toc(timerVal));
end

end      % main function psf_mrsi.m


%% sub-functions
% Cartesian grid
function xx = sub_grid(mtx,dim)
for l=1:dim, x{l} = (-mtx(l)/2:mtx(l)/2-1); end

xx = zeros(1,prod(mtx),dim);
switch dim
    case 1, xx = x{l};
    case 2
        [x1,x2] = ndgrid(x{1},x{2});
        xx(1,:,1) = x1(:).';
        xx(1,:,2) = x2(:).';
    case 3
        [x1,x2,x3] = ndgrid(x{1},x{2},x{3});
        xx(1,:,1) = x1(:).';
        xx(1,:,2) = x2(:).';
        xx(1,:,3) = x3(:).';
    otherwise, error('dim(=%g)~=1,2 or 3');
end

% xx = permute(xx,[2 3 1]);
xx = permute(xx,[3 2 1]);

end      % sub_grid
