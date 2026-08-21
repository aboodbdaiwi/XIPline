function [grad,bw,fov,desc,N,params,n_kpts] = read_ak_wav(fname,verb)
%READ_AK_WAV  Reads waveforms stored to external file using Stanford format
% [grad,bw,fov,desc,N,params,n_kpts] = read_ak_wav(fname,verb)
% 
%      fname  Gradient waveform (wav) file name
%       verb  Print info                                             (true)
%       grad  Gradient waveforms with 4us time resolution    [T/m]
%             [#pts/interleave,#interleaves,#groups]
%             with: #groups = 2 for 2d-imaging, =3 for 3d-imaging
%         bw  Full bandwidth (opuser0 = 2d3*oprbw)           [Hz]
%        fov  Field-of-view (old=nucleus; new=1H equivalent) [m]
%       desc  Description string (256 chars)
%          N  N.gpts   # input gradient pts/interleave
%             N.groups # groups
%             N.intl   # interleaves
%             N.params # parameters
%     params  Header file parameters (scanner units)
%             [grad_type fov N.intl(1) gmax N.gpts gdt n_kpts kdt 0 ...
%              esp_time esp_dir wf4pg(1:3)]
%
%  8/2024 Rolf Schulte
if (nargin<1), help(mfilename); return; end


%% input parameters
if ~exist('verb','var'), verb = []; end
if isempty(verb),        verb = true; end
if ~exist(fname,'file')
    if exist([fname '.wav'],'file')
        fname = [fname '.wav'];
    else
        error('file not found: fname=''%s''',fname);
    end
end
if ~isempty(regexpi(fname,'\.mat$','once'))
    warning('fname (=''%s'') ending with .mat; wrong file?',fname)
end


%% Read header + waveforms
fid = fopen(fname,'r','b');
if (fid==-1), error(['Could not read file: ',fname]); end
desc = strtrim(fread(fid, 256, 'char=>char').');
N.gpts = fread(fid, 1, 'uint16');
N.groups = fread(fid, 1, 'uint16');
N.intl = fread(fid, N.groups, 'uint16');
N.params = fread(fid, 1, 'uint16');
params = fread(fid, N.params, 'float64');
wave = fread(fid, inf, 'int16');
fclose(fid);


%% waveform checks
if (abs(max(abs(wave))/32767 - 1)>0.01)
    warning('max(abs(wave))=%g not scaled to 32767',max(abs(wave)));
end
if any(diff(N.intl)), error('N.intl must be the same for each group'); end
intl = N.intl(1);


%% reshaping and scaling to SI unites [T/m]
grad = reshape(wave,[N.gpts intl N.groups]);
if any(grad(end,:,:)~=1), warning('waveform not ending with 1'); end
grad(end,:,:) = 0;                 % setting stopbit to zero
grad = grad/32767*params(4)/100;   % scaling to SI units [T/m]

bw = 1d6/params(8);                % full bandwidth [Hz]
fov = params(2)/100;               % field-of-view [m]
n_kpts = params(7);                % #acq points

%% print info
if verb
    pstr = {'grad_type','fov','N.intl','gmax','N.gpts', ...
        'gdt','n_kpts','kdt','0','esp_time','esp_dir','wf4pgX','wf4pgY','wf4pgZ'};

    fprintf('description = ''%s''\n',desc);
    fprintf('#gpts   = %g (grad points per interleave)\n',N.gpts);
    fprintf('#groups = %g (grad directions)\n',N.groups);
    fprintf('#intl   = %g (interleaves)\n',intl);
    fprintf('#params = %g\n',N.params);
    fprintf('params  = %s\n',num2str(params(:).','%.3g '));
    for ll=1:N.params
        fprintf('%s = %g\n',pstr{ll},params(ll));
    end
    if N.params~=length(pstr)
        warning('N.params(=%g)~=length(pstr)(=%g)',N.params,length(pstr));
    end

end      % read_ak_wav.m
