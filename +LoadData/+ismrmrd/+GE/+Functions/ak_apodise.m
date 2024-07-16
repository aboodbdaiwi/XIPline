function d = ak_apodise(d,k,lbk,t,lbt,plt)
%AK_APODISE Spatial and/or spectral apodisation of arbitrary k-space data
%    d = ak_apodise(d,k,lbk,t,lbt,plt)
%                                                               (default)
%    d   data (#coil,#samples); opt 3rd dim=#time pts
%    k   k-space positions (1,#samples) (max(abs(k))<0.5)       ([])
%        []->no spatial apo
%  lbk   Gaussian line-broadening k-space                       (0)
%    t   time (#arms,#sample/arm)                               ([])
%        []->no spectral apodisation
%  lbt   Gaussian line-broadening time [Hz]                     (5)
%  plt   Plot linebroadenings                                   (false)
%
% 10/2020 Rolf Schulte
if (nargin<1), help(mfilename); return; end


%% default input
if ~exist('k','var'),   k = []; end
if ~exist('lbk','var'), lbk = []; end
if isempty(lbk),        lbk = 0; end
if ~exist('t','var'),   t = []; end
if ~exist('lbt','var'), lbt = []; end
if isempty(lbt),        lbt = 0; end
if ~exist('plt','var'), plt = false; end


%% misc
[~,~,npts,~] = size(d);
if plt, figure; end
if ~isempty(k)
    if max(abs(k(:)))>0.5
        warning('max(abs(k))(=%g)>0.5',max(abs(k(:))));
    end
    if size(d,2)~=size(k,2)
        warning('size(d,2)(=%d)~=size(k,2)(=%d)',size(d,2),size(k,2));
    end
end



%% spatial apodisation
if ~isempty(k) && abs(lbk)>1d-10
    if isreal(k)
        abs_k = sqrt(sum(k.^2,3));
    else
        abs_k = abs(k);
    end
    fil = exp(-(pi*lbk*abs_k).^2/(4*log(2))); % RFS 10/2020: 2*log(2)->4*log(2)
    if plt
        if ~isempty(t); subplot(2,1,1); end
        plot(abs(k),fil,'bx');
        xlabel('abs(k)');
        ylabel('filter');
        axis([0 0.5 0 1]);
    end
    d = bsxfun(@times,d,fil);
end


%% spectral apodisation
if ~isempty(t) && abs(lbt)>1d-10
    nrep = npts/size(t,1);
    if abs(nrep-floor(nrep))>1d-10
        error('nrep(=%g): npts(=%g) not multiple of size(t,1)(=%g)',...
            nrep,npts,size(t,1));
    end
    if prod(size(t))<size(d,2)
        warning('prod(size(t))(=%d)<size(d,2)(=%d)',...
            prod(size(t)),size(d,2));
    end
    if size(t,2)>size(d,2)
        warning('size(t,2)(=%d)>size(d,2)(=%d): truncating t',...
            size(t,2),size(d,2));
        t = t(:,1:size(d,2));
    end
    
    tt = reshape(t.',[1 size(t.')]);
    fil = exp(-(pi*lbt*tt).^2/(4*log(2))); % RFS 10/2020: 2*log(2)->4*log(2)
    if plt
        if ~isempty(k); subplot(2,1,2); end
        plot(squeeze(tt),squeeze(fil));
        xlabel('t');
        ylabel('filter');
        axis([0 max(t(:)) 0 1]);
    end
    d = bsxfun(@times,d,fil);
end

end      % main function ak_apodise.m
