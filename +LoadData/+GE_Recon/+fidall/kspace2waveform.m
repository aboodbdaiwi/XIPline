function wf = kspace2waveform(fname,h,wfn,nozf)
%KSPACE2WAVEFORM  Convert spiral.e kspfile to recon_spiral.m waveform file
%
%      wf = kspace2waveform(fname,h,wfn)
%   fname   Input k-space filename                              ('kspfile')
%       h   Data header structure
%     wfn   Output waveform filename
%      wf   Waveform structure
%    nozf   Skip zero-filling; calculate npix from kmax
%
% 10/2023 Rolf Schulte
if (nargin<1) && (nargout<1), help(mfilename); return; end


%% iput args
if ~exist('fname','var'),  fname = ''; end
if isempty(fname),         fname = 'kspfile'; end
if ~exist('h','var'),      h = []; end
if ~exist('wfn','var'),    wfn = ''; end
if ~exist('nozf','var'),   nozf = []; end
if isempty(nozf),          nozf = true; end


%% reading data
[wf.k,wf.dcf] = read_kspace(fname);


%% fill other wf fields
if isempty(h)
    warning('isempty(h): wf fields incomplete');
    npix = 256;
else
    bw = h.rdb_hdr.bw*2d3;
    wf.ind = true(h.image.dim_Y,h.image.dim_X);
    wf.t = repmat((0:(h.image.dim_X-1))/bw,[h.image.dim_Y 1]);
    npix = h.rdb_hdr.im_size;
end

if nozf
    kmax = max(abs(wf.k));
    wf.npix = ceil(256 * kmax)*2;
else
    wf.npix = npix;
end


%% k-space scaling
if wf.npix ~= 256
    wf.k = wf.k * 256/wf.npix;
end


%% info
kmax = max(abs(wf.k));
fprintf('kmax = %g\n',kmax);
if (kmax > 0.501), warning('kmax(=%g) > 0.501\n',kmax); end


%% writing to file
if ~isempty(wfn), save(wfn,'-struct','wf'); end


end      % kspace2waveform.m
