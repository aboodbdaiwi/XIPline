function dd=remove_spikes(d,fac_thresh,plt)
%REMOVE_SPIKES Remove spike artefacts from FIDs with gradient modulations.
%  Find spikes by looking for outliers below diff(abs(d)), with
%  Threshold = 2*std(d)*fac_thresh.
%  Detected spikes are widened by n0=2 and data is replaced by interpolated 
%  neighbouring data.
%
%dd=remove_spikes(d,fac_thresh,plt)
%           dd  Cleaned up data
%            d  data (>=2D); detection on 2nd dimension only
%   fac_thresh  Factor for threshold (default=1)
%               higher->detect less; lower->detect more spikes
%          plt  if true -> plot (default=true)
%
% 08/2009  Rolf Schulte
if (nargin<1), help(mfilename); return; end;

n0 = 2;       % #zeroes to add around spikes

if ~exist('fac_thresh','var'), fac_thresh = []; end
if isempty(fac_thresh), fac_thresh = 1; end
if ~exist('plt','var'), plt = true; end
si = size(d);
[n1,n2,n3] = size(d);
dd = d;                             % initialise matrix
if length(si)<2, error('data must be at least 2d'); end
% if si(1)>si(2), warning('remove_spikes:si','si(1)>si(2)'); end

std_data = 2*std(d(:));             % offset for plotting
x = (1:si(2));                      % x axis
thresh = std_data*fac_thresh;       % threshold
fprintf('Threshold=%g\n',thresh);
nspi = 0;
for l3=1:n3
    for l1=1:n1,
        tmp = diff(abs(d(l1,:,l3)));
        tmp(1,1:5) = 0; % exclude first points; FID mistaken for spikes
        tmp = tmp-mean(tmp);
        ind = tmp<-thresh;        % index list to spikes
        sum_ind = sum(ind);
        nspi = nspi+sum_ind;

        % negated and widened index list (spikes=false)
        ind2 = true(1,si(2));
%         for l2=(n0+1):(si(2)-n0),
%             if ind(l2),
%                 ind2(1,l2+(-n0:n0)) = false;
%             end
%         end
        for l0=-n0:n0,
            ind2(1,(n0+1):(n2-n0-1)) = ((~ind(1,(n0+l0+1):(n2-n0+l0-1))) & ...
                (ind2(1,(n0+1):(n2-n0-1))));
        end
        if ind2(1,end)==false,
            ind2(1,end) = true;
            d(l1,end) = 0;
        end
        if ind2(1,1)==false,
            ind2(1,1) = true;
            d(l1,1) = 0;
        end
        
        % replace detected spikes by interpolated neighbouring data
        if sum_ind>0,
            dd(l1,:,l3) = interp1(x(ind2),d(l1,x(ind2),l3),x);
        end
        
        if plt,  % plotting
            if ((n3>1)&&(l1==1)), figure; end
            plot(x,abs(d(l1,:,l3))+(l1-1)*std_data,'b',...
                x,abs(dd(l1,:,l3))+(l1-1)*std_data,'r',...
                x,ind2*pltoff*0.03+(l1-1)*std_data,'m');
            if l1==1, hold on; end
        end
    end
end
if plt, hold off; end
if any(isnan(dd(:))), 
    warning('rs:dd','dd contains NaN'); 
end
fprintf('Detected and removed %g spikes\n',nspi);
