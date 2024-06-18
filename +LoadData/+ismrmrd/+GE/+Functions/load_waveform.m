function wf = load_waveform(wfn,do_single,verb)
%LOAD_WAVEFORM  Load 3D waveform and recreate full k+dcf+ind+t
%       wf = load_waveform(wfn,do_single,verb)
%      wfn   Waveform mat file corresponding to ak_grad
%do_single   Use single precision                               (true)
%     verb   Verbose mode                                       (true)
%       wf   Waveform structure
%
% 2/2021  Rolf Schulte
if (nargin<1), help(mfilename); return; end


%% default parameters
if ~exist('do_single','var'), do_single = []; end
if isempty(do_single),        do_single = true; end
if ~exist('verb','var'),      verb = []; end
if isempty(verb),             verb = true; end


%% reading waveform file
if isempty(wfn), error('wfn empty'); end
if isstruct(wfn)
    wf = wfn;
else
    if ~isempty(regexpi(wfn,'\.wav$')), wfn = wfn(1:end-4); end
    if isempty(regexpi(wfn,'\.mat$')),  wfn = [wfn '.mat']; end
    if ~exist(wfn,'file'), error('file not found: wfn=%s',wfn); end
    wf = load(wfn);                    % load waveform
end


%% check for fields required for reconstruction
if ~isfield(wf,'mtx')
    fprintf('wf.mtx not existing\n'); 
end
if isfield(wf,'nviews')
    % reduced mat size: only single view stored; regenerate info
    if ~isfield(wf,'k')
        if ~isfield(wf,'ks') && ~isfield(wf,'phi') && ~isfield(wf,'theta')
            error('wf.ks,phi,theta not existing');
        else
            ks = squeeze(wf.ks);
            nks = size(ks,1);
            if do_single, ks = single(ks); end
            k = zeros(nks,wf.nviews,3,class(ks));
            
            sinphi = sind(wf.phi);
            cosphi = cosd(wf.phi);
            sintheta = sind(wf.theta);
            costheta = cosd(wf.theta);
            
            k(:,:,1) =  ks(:,1)*(costheta.*cosphi) + ...
                ks(:,2)*(sinphi) - ...
                ks(:,3)*(sintheta.*cosphi);
            k(:,:,2) =  -ks(:,1)*(costheta.*sinphi) + ...
                ks(:,2)*(cosphi) + ...
                ks(:,3)*(sintheta.*sinphi);
            k(:,:,3) = ks(:,1)*(sintheta) + ...
                ks(:,3)*(costheta);
            wf.k = reshape(k,[1 nks*wf.nviews 3]);
        end
    end
    if ~isfield(wf,'dcf'), wf.dcf = []; end
    if isempty(wf.dcf)
        if ~isfield(wf,'dcfs')
            wfn_dcf = [wfn(1:end-4) '_dcf.mat'];
            if exist(wfn_dcf,'file')
                if verb>0
                    fprintf('Loading dcf from file ''%s''\n',wfn_dcf);
                end
                load(wfn_dcf);
            else
                if verb>0
                    fprintf('Calculating dcf via calc_sdc3\n');
                    fprintf('\tand storing to ''%s''\n',wfn_dcf);
                end
                dcf = calc_sdc3(wf.k,wf.mtx,5);
                save(wfn_dcf,'dcf');
            end
            wf.dcf = dcf;
        else
            dcft = [];
            if isfield(wf,'dcft'), dcft = wf.dcft; end
            dcfs = wf.dcfs(:).';
            if do_single
                dcfs = single(dcfs);
                dcft = single(dcft);
            end
            if ~isempty(dcft)
                wf.dcf = dcfs.'*dcft(:).';
                wf.dcf = wf.dcf(:).';
            else
                wf.dcf = repmat(dcfs,[1 wf.nviews]);
            end
        end
    else
        dcf = wf.dcf;
        if isa(dcf,'uint16')
            if do_single
                wf.dcf = single(dcf)/(2^16-1)*wf.dcf_range;
            else
                wf.dcf = double(dcf)/(2^16-1)*wf.dcf_range;
            end
        end
    end
    if ~isfield(wf,'ind')
        if ~isfield(wf,'inds')
            error('wf.inds not existing'); 
        else
            wf.ind = repmat(wf.inds,[wf.nviews 1]);
        end
    end
    if ~isfield(wf,'t')
        if ~isfield(wf,'ts')
            error('wf.ts not existing'); 
        else
            ts = wf.ts(:).';
            if do_single, ts = single(ts); end
            wf.t = repmat(ts,[1 wf.nviews]);
        end
    end
end

if isempty(wf.dcf),    error('dcf is empty'); end
if any(isnan(wf.dcf)), error('dcf contains NaN'); end
if any(isinf(wf.dcf)), error('dcf contains inf'); end


end      % load_waveform.m
