function d = truma(d,resmpl,new_mtx,lb,verb)
%TRUMA Centered truncation or zero-filling of n-dimensional matrix
%      including FFT -> re-sampling
%    [dd] = truma(d,resmpl,new_mtx,lb,verb)
%                                                          (default)
%       d   Input matrix data (max 10D)
%  resmpl   Truncate/zero-fill Fourier transformed data    (false)
% new_mtx   Matrix size to truncate/zero-fill to           (min(size(d)))
%      lb   Gaussian apodisation (if resmpl==true)         ([])
%    verb   Verbose mode                                   (false)
%
%      dd   Data to min(size(d)) + centred
%
% 11/2020 Rolf Schulte
if nargin<1, help(mfilename); return; end


%% default input parameters
if ~exist('resmpl','var'),   resmpl = []; end
if isempty(resmpl),          resmpl = false; end
if ~exist('new_mtx','var'),  new_mtx = []; end
if ~exist('lb','var'),       lb = []; end
if ~exist('verb','var'),     verb = []; end
if isempty(verb),            verb = false; end


%% misc parameters + checks
if resmpl
    if verb, fprintf('Downsampling data (fftn before truncation)\n'); end
    data_is_real = false;
    if isreal(d), data_is_real = true; end  % remember data type
end
mtx = size(d);                              % matrix size
dim = length(mtx);                          % dimensions of data
if dim<2, error('dim(=%g)<2',dim); end      % non-sense
if dim>10, error('dim(=%g)>10',dim); end    % limited to 10D
if verb, fprintf('Input matrix size = '); disp(mtx); end
if isempty(new_mtx)
    nn = min(mtx);                          % smallest matrix
    new_mtx = repmat(nn,[1 dim]);           % new matrix size
else
    new_mtx = new_mtx(:).';
    if length(mtx)~=length(new_mtx)
        if any(new_mtx(1,((length(mtx)+1):length(new_mtx)))~=1)
            warning('length(mtx)(=%g)~=length(new_mtx)(=%g); truncating new_mtx',...
                length(mtx),length(new_mtx));
        end
        if length(mtx)<length(new_mtx)
            new_mtx = new_mtx(1,1:length(mtx));
        else
            warning('length(mtx)(=%d)>length(new_mtx)(=%d): adding mtx to new_mtx',...
                length(mtx),length(new_mtx));
            ii = (length(new_mtx)+1):length(mtx);
            new_mtx = [new_mtx mtx(1,ii)];
        end
    end
    if any(new_mtx>mtx)
        % warning('any(new_mtx>mtx)');
        zf = new_mtx;
        new_mtx = min([mtx;new_mtx],[],1);
    end
end
if verb && (any(new_mtx<mtx))
    fprintf('Truncate data to'); disp(new_mtx); 
end
if verb && (min(new_mtx)<2), fprintf('min(new_mtx)(=%g)<2',min(new_mtx)); end

for l=1:10,  ind{l} = 1; end                % initialise index list
for l=1:dim, ind{l} = (1:new_mtx(l))+floor((mtx(l)-new_mtx(l))/2); end
if resmpl, d = ifftshift(ifftn(d)); end  % downsampling
% if resmpl, d = ifftshift(ifftn(fftshift(d))); end % downsampling


%% actual truncation
d = d(ind{1},ind{2},ind{3},ind{4},ind{5},ind{6},ind{7},ind{8},ind{9},ind{10});


%% if new_mtx>mtx: zerofilling
if exist('zf','var')
    if verb, fprintf('Zerofill data to'); disp(zf); end
    for l=1:10,  izf{l} = 1; end                % initialise index list
    for l=1:dim, izf{l} = (1:new_mtx(l))+floor((zf(l)-new_mtx(l))/2); end

    dtmp = zeros(zf);
    dtmp(izf{1},izf{2},izf{3},izf{4},izf{5},izf{6},izf{7},izf{8},izf{9},izf{10}) = d;
    d = dtmp;
    clear dtmp;
end


%% Gaussian apodisation
if resmpl
    if ~isempty(lb)
        if length(lb)~=1
            warning('length(lb)(=%g) =~ 1; taking only lb(1)',length(lb)); 
        end
        for l=1:dim
            nn = size(d,l);
            lbf = lb_fun(nn,nn,[0 lb(1)],'h').';
            d = bsxfun(@times,d,reshape(lbf,[ones(1,l-1) nn 1]));
        end
    end

    d = fftn(fftshift(d));                     % FFT back
    % d = ifftshift(fftn(fftshift(d)));            % FFT back
    if verb>1
        fprintf('RMS(real(d(:)))=%g; RMS(imag(d(:)))=%g\n',...
            sqrt(mean(real(d(:)).^2)),sqrt(mean(imag(d(:)).^2)));
    end
    if data_is_real, d = real(d); end       % discard imag (small error)
end


end           % end of main function (truma)
