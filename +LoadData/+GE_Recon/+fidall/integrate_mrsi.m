function [area,iint] = integrate_mrsi(spec,fax,freqs,nint,verb)
%INTEGRATE_MRSI  Quatify MRSI by integrating peaks
%  area = integrate_mrsi(spec,fax,freqs,nint,verb)
%  spec   Reconstructed MRSI spectra (nspec,nx,ny,nz)
%   fax   Frequency axis
% freqs   Frequencies to integrate
%  nint   Integration boundary (points)
%  verb   Verbose mode                                          (true)
%
%  area   Quantified signal
%  iint   Index to integration boundary
%
% 12/2020 Rolf Schulte
if nargin<1, help(mfilename); return; end
if nargin<4, error('nargin(=%d)<4',nargin); end


%% input parameters
if ~exist('verb','var'), verb = []; end
if isempty(verb),        verb = true; end
fax = fax(:).';
if length(freqs)~=size(nint,1)
    warning('length(freqs)(=%d)~=size(nint,1)(=%d)',length(freqs),size(nint,1));
end
if size(nint,2)~=2, warning('size(nint,2)(=%d)~=2',size(nint,2)); end


%% reshape spectra + spatial-rms
nfreq = length(freqs);
[nspec,nx,ny,nz] = size(spec);
nn = nx*ny*nz;
ss = reshape(spec,[nspec nx*ny*nz]);
spec_rms = sqrt(mean(abs(ss.*conj(ss)),2));


%% find maximum of spec near freqs(1) by walking uphill
imax0 = sub_find_max(spec_rms,fax,freqs(1),verb);
imax = zeros(1,nn);
for lx=1:nn
    imax(1,lx) = sub_find_max(abs(ss(:,lx)),fax,fax(imax0),false);
end
idiff = zeros(nfreq,1);
for l=1:nfreq, [~,idiff(l,1)] = min(abs(fax-freqs(l))); end
idiff = idiff-idiff(1);
if verb, fprintf('idiff = '); disp(idiff.'); end


%% integrate around peak
do_warn1 = true; do_warn2 = true;
area = zeros(nfreq,nn);
for lf=1:nfreq
    for lx=1:nn
        ii = imax(1,lx)+(nint(lf,1):nint(lf,2))+idiff(lf);
        if any(ii>nspec)
            if do_warn1, warning('ii>nspec'); end
            ii = ii(ii<=nspec);
            do_warn1 = false;
        end
        if any(ii<1)
            if do_warn2, warning('ii<1'); end
            ii = ii(ii>0); 
            do_warn2 = false;
        end
        area(lf,lx) = sum(abs(ss(ii,lx)),1);
    end
end

%% index to integration boundary
if nargout>1
    iint = zeros(nfreq,2);
    for lf=1:nfreq
        iint(lf,:) = imax0(1,1)+nint(lf,:)+idiff(lf);
    end
    if any(iint(:)>nspec)
         warning('iint>nspec'); disp(iint)
         iint(iint>nspec) = nspec;
    end
    if any(iint(:)<1)
        warning('iint<1'); disp(iint)
        iint(iint<1) = 1;
    end
end


%% reshape to original spatial dimensions
area = reshape(area,[nfreq nx ny nz]);


%% plotting
if verb>1
    imagesc_row(permute(area,[2 3 4 1]))
end

end      % recon_xenon.m


%% sub-functions
function imax = sub_find_max(spec,fax,freq,verb)
% find maximum near freq by walking uphill
nspec = size(spec,1);
[~,imax] = min(abs(fax-freq));
if verb, fprintf('Initial imax = %d\n',imax); end
l = 0;
while true
    l = l+1;
    if spec(imax+1,1)>spec(imax,1)
        imax = imax+1;
    else
        if spec(imax-1,1)>spec(imax,1)
            imax = imax-1;
        else
            break;
        end
    end
    if imax==1, warning('imax==1; break'); break; end
    if imax==nspec, warning('imax==nspec(=%d); break',nspec); break; end
    if l>nspec, warning('l(=%d)>nspec(=%d); break',l,nspec); break; end
end
if verb
    fprintf('Final imax = %d\n',imax); 
    fprintf('#steps = %d\n',l);
end

end      % sub_find_max
