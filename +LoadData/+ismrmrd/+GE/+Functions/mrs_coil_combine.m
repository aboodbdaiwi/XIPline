function dd = mrs_coil_combine(d,ind1,pts4svds,n_noise,verb)
%MRS_COIL_COMBINE  Coil combination for single voxel MRS
%   Determines phases and weights of individual coil elements from first 
%   points of FID/echo, phases and weights channels, and then averages them
%
% dd = mrs_coil_combine(d,ind1,ind2,verb)
%         d  Raw p-file data [nf,ns,n3,n4,n5,nc]    (unchopped)
%            nf=#frames; ns=#spectral pts; nc=#coils
%      ind1  Index list for averaging along frames         (default=nf)
%  pts4svds  #largest samples for svds                     (default=256)
%   n_noise  #samples along t2 for noise determination     (default=ns/2)
%      verb  Verbose mode                                  (default=true)
%
%        dd  coil-combined FID/echo [nf,ns,n3,n4,n5]
%
%  2/2020 Rolf Schulte
if nargin<1, help(mfilename); return; end


%% miscelaneous parameters
[nf,ns,n3,n4,n5,nc] = size(d);

if ~exist('ind1','var'),     ind1 = []; end
if isempty(ind1),            ind1 = (1:nf); end
if ~exist('pts4svds','var'), pts4svds = []; end
if isempty(pts4svds),        pts4svds = 16; end
if ~exist('n_noise','var'),  n_noise = []; end
if ~exist('verb','var'),     verb = []; end
if isempty(verb),            verb = true; end

% parameter checks
if any(ind1<1),   warning('ind1<1');  ind1 = ind1(ind1>0); end
if any(ind1>nf),  warning('ind1>nf'); ind1 = ind1(ind1<=nf); end
if (pts4svds<10), warning('pts4svds<10'); end
if (pts4svds>ns), warning('pts4svds>ns'); end


%% coil combination
if nc>1
    [~,ind] = sort(sqrt(sum(sum(sum(sum(sum(d.*conj(d),1),3),4),5),6)));
    if pts4svds<ns
        ind_svds = sort(ind((ns-pts4svds+1):ns));
    else
        ind_svds = (1:ns);
    end
    di = mean(d(ind1,:,:,:,:,:),1);
    ww = zeros(1,1,n3,n4,n5,nc);
    for l3=1:n3
        for l4=1:n4
            for l5=1:n5
                A = squeeze(di(1,ind_svds,l3,l4,l5,:));
                [U,S,V,flag] = svds(A,1);
                ww(1,1,l3,l4,l5,:) = V;
                if flag~=0
                    fprintf('l3=%g; l4=%g; l5=%g; flag=%g\n',l3,l4,l5,flag);
                end
            end
        end
    end
    ww = bsxfun(@rdivide,ww,sum(abs(ww),6));
    dd = sum(bsxfun(@times,d,ww),6);
    
else
    if verb, fprintf('Single channel data; skipping coil combination\n'); end
    dd = d;
end

if verb
    ddi = mean(dd,1);
    if isempty(n_noise), n_noise = ceil(ns/2); end
    for l3=1:n3
        for l4=1:n4
            for l5=1:n5
                tmp = sort(abs(ddi(:,:,l3,l4,l5)),2,'descend');
                signal = mean(tmp(:,1:10),2);
                noise = std(tmp(:,n_noise:end));
                fprintf('SNR = %g/%g = %g\n',signal,noise,signal/noise);
            end
        end
    end
end

end      % mrs_coil_combine.m
