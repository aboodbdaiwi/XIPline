function par = header2par(h)
% HEADER2PAR Convert GE MR header structure into par structure
%   par = header2par(header)
%
%8/2009 Rolf Schulte
if (nargin<1), help(mfilename); return; end;

if ~isfield(h,'rdb_hdr'),
    warning('header2par:rdb_hdr','rdb_hdr not a field: setting par=h');
    par = h;
    par.bw = h.sample_frequency;
    par.f0 = h.synthesizer_frequency;
    return;
end

% par section as for PMS convention
par.samples = h.rdb_hdr.frame_size;                 % sampling points
par.sample_frequency = h.rdb_hdr.spectral_width;    % sampling freq [Hz]
if par.sample_frequency==0,    % if converted MRI psd
    par.sample_frequency = h.image.vbw*2d3;         % 2*MRI BW = fullBW [Hz]
end
par.offset_frequency = h.rdb_hdr.default_offset;    % offset freq
f0 = h.rdb_hdr.ps_mps_freq/10;
par.synthesizer_frequency = f0;   % f0=B0*gamma/2/pi
par.echo_time = h.rdb_hdr.te/1d3;                   % echo time
par.averages  = h.rdb_hdr.navs;                     % averages (opnex)
par.receivers = (h.rdb_hdr.dab(2)-h.rdb_hdr.dab(1))+1; % receiver channels
par.phase_cycling = h.rdb_hdr.navs;                 % excitations per phase cycle
par.ref_rows  = h.rdb_hdr.user19;                   % First nx-rows: ref scan
par.rows      = h.rdb_hdr.nframes-par.ref_rows;     % actual frames
par.nslices   = h.rdb_hdr.nslices;
par.repetition_time = h.image.tr*1d-6;

% par section as required for gridding (p2grid function)
par.chp = ischop(h);
par.bw = par.sample_frequency;
% if (par.bw~=h.rdb_hdr.user0), 
%     warning('header2par:bw',...
%         'h.rdb_hdr.spectral_width (%g) ~= h.rdb_hdr.user0 (%g)',...
%         par.bw,h.rdb_hdr.user0);
% end
par.r1 = h.rdb_hdr.ps_mps_r1;
par.r2 = h.rdb_hdr.ps_mps_r2;
par.tg = h.rdb_hdr.ps_mps_tg;
par.f0 = f0;

if strcmpi(h.image.psd_iname,'fidak')
    par.tacq = h.rdb_hdr.user25*1d-3;  % acq delay (FID -> linear phase)
end

if strcmpi(h.image.psd_iname,'jpress') &&  h.image.user24==1,
    % 1 (no add) Store each single data acquisition in pfile.
    par.no_add = true;
    par.rows     = h.rdb_hdr.nframes/h.image.nex-par.ref_rows;  % actual frames
    if rem(par.rows,1)~=0, par.rows; end
end
