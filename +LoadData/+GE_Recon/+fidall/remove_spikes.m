function dd=remove_spikes(d,fac_thresh,n0,nexcl,plt)
%REMOVE_SPIKES Remove spike artefacts from FIDs with gradient modulations.
%  Find spikes by looking for outliers below diff(abs(d)), with
%  Threshold = 2*std(d)*fac_thresh.
%  Detected spikes are widened by n0 and data is replaced by interpolated 
%  neighbouring data.
%
% dd = remove_spikes(d,fac_thresh,n0,nexcl,plt)
%         d  data (>=2D); detection on 2nd dimension only
%fac_thresh  Factor for threshold                                       (1)
%            higher->detect less; lower->detect more spikes
%        n0  #points to interpolate around spikes                       (2)
%     nexcl  #points to skip at start of FID (to avoid nulling data)   (20) 
%       plt  Plotting                                               (false)
%        dd  Cleaned up data
%
% 12/2023 Rolf Schulte
if (nargin<1), help(mfilename); return; end


%% default parameters
if ~exist('fac_thresh','var'), fac_thresh = []; end
if isempty(fac_thresh), fac_thresh = 1; end
if ~exist('n0','var'),  n0 = []; end
if isempty(n0),         n0 = 2; end
if ~exist('nexcl','var'), nexcl = []; end
if isempty(nexcl),      nexcl = 20; end

if ~exist('plt','var'), plt = []; end
if isempty(plt),        plt = false; end


%% checks
si = size(d);
[n1,n2,n3] = size(d);
dd = d;                             % initialise matrix
if length(si)<2, error('data must be at least 2d'); end
% if si(1)>si(2), warning('si(1)>si(2)'); end

std_data = 2*std(d(:));             % offset for plotting
x = (1:si(2));                      % x axis
thresh = std_data*fac_thresh;       % threshold
fprintf('remove_spikes: threshold=%g\n',thresh);
nspi = 0;


%% actual processing
for l3=1:n3
    for l1=1:n1
        tmp = diff(abs(d(l1,:,l3)));
        tmp(1,1:nexcl) = 0; % exclude first points; FID mistaken for spikes
        tmp = tmp-mean(tmp);
        ind = tmp<-thresh;        % index list to spikes
        sum_ind = sum(ind);
        nspi = nspi+sum_ind;
        
        % negated and widened index list (spikes=false)
        ind2 = true(1,si(2));
        for l0=-n0:n0
            ind2(1,(n0+1):(n2-n0-1)) = ((~ind(1,(n0+l0+1):(n2-n0+l0-1))) & ...
                (ind2(1,(n0+1):(n2-n0-1))));
        end
        if ind2(1,end)==false
            ind2(1,end) = true;
            d(l1,end) = 0;
        end
        if ind2(1,1)==false
            ind2(1,1) = true;
            d(l1,1) = 0;
        end
        
        % replace detected spikes by interpolated neighbouring data
        if sum_ind>0
            dd(l1,:,l3) = interp1(x(ind2),d(l1,x(ind2),l3),x);
        end
    end       % for l1=1:n1
end           % for l3=1:n3


%% plotting
if plt
    if n3>1, warning('n3(=%d)>1: plotting only first element',n3); end
    clf
    subplot(2,1,1);
    plot(abs(d(:,:,1)).');
    hold on; plot([1 size(d,2)],[1 1]*thresh,'r:'); hold off
    axis([1 n2 0 1.05*max(max(abs(d(:,:,1))))]); grid on
    
    subplot(2,1,2); 
    plot(abs(dd(:,:,1)).');
    hold on; plot([1 size(d,2)],[1 1]*thresh,'r:'); hold off
    axis([1 n2 0 1.05*max(max(abs(dd(:,:,1))))]); grid on
end


%% checks
if any(isnan(dd(:)))
    warning('dd contains NaN; zeroeing'); 
    dd(isnan(dd)) = 0;
end
if any(isinf(dd(:)))
    warning('dd contains inf; zeroeing'); 
    dd(isinf(dd)) = 0;
end


%% print info
fprintf('Detected and removed %g spikes\n',nspi);


end      % remove_spikes.m
