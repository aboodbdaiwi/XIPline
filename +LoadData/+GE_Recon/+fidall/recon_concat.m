function out=recon_concat(dd,hh,cci_fname,fname,delay,mtx,lb)
%RECON_CONCAT Start reconstruction for concatenated trajectories
% recon_concat(dd,hh,cci_fname,fname,delay,mtx,lb)
%       dd   Raw data or filename
%       hh   Data header (or empty)
%cci_fname   File with concatenation info
%    fname   Output filename
%    delay   Gradient delay time [us]
%      mtx   Matrix size to reconstruct to
%       lb   Apodisation [Hz]
%      out   output data struct with bb and bbabs
%
%  7/2024 Rolf Schulte
if (nargin<1), help(mfilename); return; end


%% input
if ~exist('delay','var'),    delay = []; end
if isempty(delay),           delay = 0; end
if ~exist('fname','var'),    fname = []; end
if ~isempty(fname)
    if ~isempty(regexpi(fname,'\.7$', 'once')), fname = fname(1:end-2); end
    if ~isempty(regexpi(fname,'\.h5$','once')), fname = fname(1:end-3); end
end
if ~exist('mtx','var'),      mtx = []; end
if isempty(mtx),             mtx = 0; end
if ~exist('lb','var'),       lb = []; end
if isempty(lb),              lb = 0; end


%% reading in data, if pfile name is given
if ~isnumeric(dd) || isempty(dd)
    if ~exist(dd,'file') && ~isempty(dd)
        warning('strange input for dd (=''%s''); file not existing?',dd);
    end
    [dd,hh] = read_p(dd);
end


%% load concatenation info
if ~isempty(regexpi(cci_fname,'\.wav$', 'once')), cci_fname = cci_fname(1:end-4); end
if isempty(regexpi(cci_fname,'\.mat$', 'once')),  cci_fname = [cci_fname '.mat']; end
if ~exist(cci_fname,'file'), error('file not found: ccifn=%s',cci_fname); end
cci = load(cci_fname);


%% parameters
nwf = length(cci.wfns);
nexc = cci.ind_exc{nwf}(end);
nrep = size(dd,1)/nexc;
if abs(nrep-floor(nrep))>1d-10, warning('nrep(=%g) not int',nrep); end
wfpath = fileparts(cci_fname);
if ~isempty(wfpath)
    if isunix
        wfpath = [wfpath '/'];
    else
        wfpath = [wfpath '\'];
    end
end


%% loop through reconstructions
out = struct();
for l1 = 1:nwf
    fprintf('**********************************\n');
    fprintf('Reconstructing (%d/%d): %s\n',l1,nwf,cci.recons{l1});
    fprintf('\twaveform = %s\n',cci.wfns{l1});
    ii1 = repmat(cci.ind_exc{l1}(:),[ceil(nrep) 1]);
    ii1 = ii1(1:size(dd,1),1);
    fprintf('ii1 = %s\n',num2str(ii1));
    d = dd(ii1,:,:,:,:,:);
    wfn = [wfpath cci.wfns{l1}];
    out(l1).wfn=wfn;
    if length(mtx)==1, mtx = repmat(mtx,nwf); end
    if length(lb)==1,  lb  = repmat(lb,nwf); end
    switch cci.recons{l1}
        case 'spiral'
            h = hh;
            seadd = l1*100;   % add to series# for dicom export
            if h.rdb_hdr.user22>0
                warning('spsp: setting cv22(=%d)&cv40(=%d) to 0',...
                    h.rdb_hdr.user22,h.rdb_hdr.user40);
                h.rdb_hdr.user22 = 0;
                h.rdb_hdr.user40 = 0;
            end
            if isempty(fname), fn = [];
            else,              fn = [fname '_2D_cc' num2str(l1)]; end
            recon_spiral(d,h,wfn,mtx(l1),delay,lb(l1),[],fn,...
                [],[],[],[],[],[],seadd);

        case 'grid'
            if isempty(fname), fn = [];
            else,              fn = [fname '_grd_cc' num2str(l1)]; end
            [out(l1).bb,out(l1).bbabs] = recon_grid(d,hh,wfn,mtx(l1),...
                delay,lb(l1),fn);

        case 'grid3d'
            if isempty(fname), fn = [];
            else,              fn = [fname '_3D_cc' num2str(l1)]; end
            [out(l1).bb,out(l1).bbabs] = recon_grid3d(d,hh,wfn,mtx(l1),...
                delay,lb(l1),[],fn);

        case 'epi'
            if isempty(fname), fn = [];
            else,              fn = [fname '_epi_cc' num2str(l1)]; end
            recon_epi(d,h,wfn,mtx(l1),delay,fn);

        otherwise
            warning('cci.recon(=%s) not found/implemented',cci.recon);
    end
end      % for l1=1:nwf


end      % recon_concat.m
