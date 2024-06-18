function [d,h]=read_pfile(fname,do_single,verb)
%READ_PFILE Read GE P-file.7 using Orchestra
%     [d,h] = read_pfile(fname,do_single,verb)
%     fname   filename                   ('' -> last .7 file in directory)
% do_single   d in single precision      (false)
%      verb   verbose mode               (true)
%
%         d   raw data matrix (same dimensionality as p-files)
%             [nframes,frame_size,nphases,nechoes,nslices,nreceivers]
%         h   header structure
%
% 2/2020  Rolf Schulte
%See also GERecon, archive2header.
if (nargout<1), help(mfilename); return; end

use_read_mr = false;              % revert back to Matlab reading routines

if ~exist('verb','var'),  verb = []; end
if isempty(verb),         verb = true; end
if ~exist('fname','var'), fname = []; end
if isempty(fname)
    tmp=ls('*.7');
    fname = tmp(end,:);
    if verb>0, fprintf('fname = ''%s''\n',fname); end
end
if ~exist(fname,'file'),  error('file ''%s'' not found',fname); end
if ~exist('do_single','var'), do_single = []; end
if isempty(do_single),        do_single = false; end


%% open archive and load header

pfile = LoadData.ismrmrd.GE.GERecon('Pfile.Load', fname);
h = LoadData.ismrmrd.GE.GERecon('Pfile.Header',pfile);
if verb
    rdbm_rev = h.RawHeader.rdbm_rev;
    fprintf('header release = %g\n',rdbm_rev);
end
psd_iname =  h.ImageData.psd_iname;
if verb, fprintf('psd_iname=%s\n',psd_iname); end
if isfield(pfile,'numSpokesLowRes')
    warning('pfile.numSpokesLowRes existing: using read_MR_rawdata.m\n');
    use_read_mr = true;
end
if ~isempty(regexpi(psd_iname,'3drad')) || ...
          ~isempty(regexpi(psd_iname,'silen'))|| ...
          ~isempty(regexpi(psd_iname,'burzte'))
    warning('psd_iname suggest ZTE data: using read_MR_rawdata.m\n');
    use_read_mr = true;
end
nrcv = (h.RawHeader.dab(2)-h.RawHeader.dab(1)+1);
if pfile.channels~=nrcv
    warning('Bug in Orchestra; not reading multi-channel data correctly')
    fprintf('pfile.channels(=%d)~=nrcv(from dab)=%d\n',pfile.channels,nrcv);
    use_read_mr = true;
end


%% reverting back to read_MR_rawdata
if use_read_mr
    if verb, fprintf('Reverting back to read_MR_rawdata\n'); end
    [d,h] = read_MR_rawdata(fname);
    if do_single, d = single(d); end
    return
end


%% Scan parameters
nreceivers = pfile.channels;
nframes = pfile.yRes;
frame_size = pfile.xRes;
nphases = pfile.phases;
nechoes = pfile.echoes;
nslices = pfile.slices;
if verb
    fprintf('nframes=%g  frame_size=%g  nphases=%g  nechoes=%g  nslices=%g  nreceivers=%g\n',...
        nframes,frame_size,nphases,nechoes,nslices,nreceivers);
end


%% pre-allocate memory for data
if do_single
    d = complex(zeros(nframes,frame_size,nphases,nechoes,nslices,nreceivers,'single'));
else
    d = complex(zeros(nframes,frame_size,nphases,nechoes,nslices,nreceivers));
end


%% reading data
for l3=1:nphases
    for l4=1:nechoes
        for l5=1:nslices
            for l6=1:nreceivers
                kspace = GERecon('Pfile.KSpace',l5,l4,l6,l3,pfile);
                d(:,:,l3,l4,l5,l6) = kspace.';
            end
        end
    end
end


%% convert archive structure to p-file header structure
if nargout>1
    h = LoadData.ismrmrd.GE.Functions.pfile2header(h);
    for lslice=1:nslices
        h.orientation{lslice} = GERecon('Pfile.Orientation',lslice);
        h.corners{lslice} = GERecon('Pfile.Corners',lslice);
    end
end


end   % main function read_pfile
