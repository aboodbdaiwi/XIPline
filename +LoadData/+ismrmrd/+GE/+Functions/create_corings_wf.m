function [wf,hh] = create_corings_wf(h,kfname,do_one_ring)
%CREATE_CORINGS_WF  Generate waveform structure for reconstruction of
% Concentric Rings MRSI acquired with Ralph's sLASER sequence
%    [wf,hh] = create_corings_wf(h,kfname,do_one_ring)
%          h   Header structure (from p-files)
%     kfname   File with trajectory information (/usr/g/bin/dwcrt.k)
%do_one_ring   index for data reconstruction of single ring 
%              (0=all; 1,2,...,nrings)
%         wf   waveform structure for recon_corings or recon_spiral
%         hh   adapt header for recon_spiral
%
% 10/2019 Rolf Schulte
if nargin<1, help(mfilename); return; end

if ~exist('do_one_ring','var'), do_one_ring = []; end
if isempty(do_one_ring),        do_one_ring = false; end


%% info from p-file header
n1 = h.rdb_hdr.nframes;                % #rings (=exc)
n2 = h.rdb_hdr.da_xres;                % #readout points/exc
mtx = h.rdb_hdr.xcsi;
n_spec = h.rdb_hdr.user1;
nk_circ = h.image.user4;
bw_spec = h.rdb_hdr.user0;
kdt = 1/(bw_spec*nk_circ);
t = (0:(n2-1))*kdt;
t = repmat(t,[1 n1]);

wf.t = t;
wf.npix = mtx;                         % matrix size
wf.mtx(1,1:2) = mtx;                   % matrix size
wf.npts = n2;                          % #sampling points
wf.n_spec = n_spec;                    % #spectral points
wf.nk_circ = nk_circ;                  % #sample points per 1 ring
wf.bw_spec = bw_spec;                  % spectral bandwidth [Hz]
wf.n_ring = n1;                        % #rings
wf.n_intlv = 1;                        % #spectral interleaves
wf.n_zpe = 1;                          % #z phase encodes
wf.dim = 2;                            % spatial dimensions
wf.nexc = n1;                          % #excitations


%% read k-space waveform file
if ~exist(kfname,'file')
    error('k-space file ''%s'' not found',kfname);
end
fp = fopen(kfname,'rt');
A = fscanf(fp, '%f', [4, inf]);
plot(A(3,:),A(4,:),'.');
fclose(fp);

k1 = permute(A([4,3],:),[3 2 1]);


%% calculate density-compensation function
dcf1 = voronoi_area(k1*mtx);
wf.dcf1 = dcf1;
wf.k1 = k1;
if do_one_ring>0
    wf.n_spec = 1;
    wf.k = k1;
    wf.dcf = dcf1;
    ind = false(n1,n2);
    if do_one_ring>1
        ii = (do_one_ring-1)*nk_circ+(1:nk_circ);
        ind(:,ii) = true;
    else
        ind(:,1:nk_circ) = true;
    end
else
    k = reshape(k1,[nk_circ n1 2]);
    k = repmat(k,[n_spec 1]);
    k = reshape(k,[1 nk_circ*n1*n_spec 2]);
    wf.k = k;
    
    dcf = reshape(dcf1,[nk_circ n1]);
    dcf = repmat(dcf,[n_spec 1]);
    dcf = dcf(:).';
    wf.dcf = dcf;
    
    ind = true(n1,n2);
end
wf.ind = ind;                          % index to select acquired data


%% adapt header
hh = h;
hh.rdb_hdr.user21 = 0;
hh.rdb_hdr.user22 = 0;


end      % main function create_corings_wf.m
