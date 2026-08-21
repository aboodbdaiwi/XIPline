function wf = load_waveform(wfn,do_single,verb)
%LOAD_WAVEFORM  Load 3D waveform and recreate full k+dcf+ind+t
%       wf = load_waveform(wfn,do_single,verb)
%      wfn   Waveform mat file corresponding to ak_grad
%do_single   Use single precision                                    (true)
%     verb   Verbose mode                                            (true)
%       wf   Waveform structure
%
%  9/2022  Rolf Schulte
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
    if ~isempty(regexpi(wfn,'\.wav$', 'once')), wfn = wfn(1:end-4); end
    if isempty(regexpi(wfn,'\.mat$', 'once')),  wfn = [wfn '.mat']; end
    if ~exist(wfn,'file'), error('file not found: wfn=%s',wfn); end
    if verb, fprintf('Loading waveform file=\n\t%s\n',wfn); end
    wf = load(wfn);                    % load waveform
end


%% check for fields required for reconstruction
if ~isfield(wf,'mtx')
    if isfield(wf,'npix')
        wf.mtx = wf.npix;
    else
        fprintf('wf.mtx and wf.npix not existing\n');
    end
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
            if ~isempty(wf.theta)
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
            else
                if isfield(wf,'kz')
                    ngz = length(wf.kz)/length(wf.phi);
                    if abs(ngz-round(ngz))>1d-10
                        warning('ngz(=%g) not even',ngz); 
                    end
                    ngz = round(ngz);
                    k(:,:,1) =  repmat(ks(:,1)*cosphi + ...
                        ks(:,2)*sinphi,[1 ngz]);
                    k(:,:,2) = repmat(-ks(:,1)*sinphi + ...
                        ks(:,2)*cosphi,[1 ngz]);
                    k(:,:,3) = repmat(wf.kz,[nks 1]);
                else
                    error('theta=[] & no wf.kz field');
                end
            end
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
                dcf = calc_sdc3(double(wf.k),wf.mtx,5);
                if do_single, dcf = single(dcf); end
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
else     % isfield(wf,'nviews')
    % wf.nviews = size(wf.ind,1);
    wf.nviews = sum(sum(wf.ind,2)>0);
end


%% reformat 2D complex k-space
if (size(wf.k,3)==1) && ~isreal(wf.k)
    kk = zeros(1,size(wf.k,2),2);
    kk(1,:,1) = real(wf.k);
    kk(1,:,2) = imag(wf.k);
    wf.k = kk;
    clear kk
end


%% checks
if isempty(wf.dcf),    error('dcf is empty'); end
if any(isnan(wf.dcf)), error('dcf contains NaN'); end
if any(isinf(wf.dcf)), error('dcf contains inf'); end

if size(wf.k,1)~=1,   warning('size(k,1)(=%d)~=1',size(wf.k,1)); end
if size(wf.dcf,1)~=1, warning('size(dcf,1)(=%d)~=1',size(wf.dcf,1)); end
if length(size(wf.k))~=3
    warning('length(size(k))(=%d)~=3',length(size(wf.k))); 
end
dim = size(wf.k,3);
switch dim
    case 2, if verb, fprintf('2D\n'); end
    case 3, if verb, fprintf('3D\n'); end
    otherwise
        warning('size(k))(=%d)~=2 or 3',size(wf.k)); 
end
if size(wf.k,2)~=size(wf.dcf,2)
    warning('size(k,2)(=%d)~=size(dcf,2)(=%d)',size(wf.k,2),size(wf.dcf,2)); 
end
if isfield(wf,'mtx')
    if length(wf.mtx)>dim
        warning('length(wf.mtx)(=%d)>dim(=%d); truncating');
        wf.mtx = wf.mtx(1:dim);
    end
    if ((length(wf.fov)==2) && (length(wf.mtx)==3))
        warning('((length(wf.fov)==2) && (length(wf.mtx)==3))');
        fprintf('replicating fov(1)\n');
        wf.fov = [wf.fov(1),wf.fov(1),wf.fov(2)];
    end
end
if size(wf.t,1)>1
    warning('size(wf.t,1)(=%d)>1; reshaping',size(wf.t,1));
    wf.t = wf.t.'; wf.t = wf.t(:).';
end


end      % load_waveform.m
