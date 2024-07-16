function dcf = calc_sdc3(k,mtx,numIter,osf,verb)
%CALC_SDC3  Calculate 3D Density Compensation Function (DCF)
% Wrapper function to sdc3_MAT.c
% Compile under Windows via "mex -compatibleArrayDims sdc3_MAT.c"
%
%     dcf = calc_sdc3(k,mtx,numIter,osf,verb)
%       k   k-Space trajectory normalised to +-0.5
%           (1,#k-pts,dim=3)
%     mtx   Matrix size
% numIter   #iterations                                  (15)
%     osf   Oversampling factor                          (2.6)
%    verb   Verbose mode                                 (false)
%
%     dcf   Density Compensation Function (1,#k-pts)
%
% Literature: mrm67_701_2012_zwart.pdf
% Download sdc3 library from: https://www.ismrm.org/mri_unbound
%  11/2020 Rolf Schulte
if (nargin<1), help(mfilename); return; end

fufa_dcf = 1/1.3;           % fudge factor to correct dcf
thold = 5;                  % crop values to threshhold 

%% input paameters
if ~exist('numIter','var'), numIter = []; end
if isempty(numIter),        numIter = 15; end
if ~exist('osf','var'),     osf = []; end
if isempty(osf),            osf = 2.6; end
if ~exist('verb','var'),    verb = []; end
if isempty(verb),           verb = false; end


%% checks
if length(mtx)~=1, error('length(mtx)(=%g)~=1',length(mtx)); end
if ~isreal(k),     error('~isreal(k)'); end
if size(k,1)~=1,   error('size(k,1)(=%g)~=1',size(k,1)); end
if size(k,3)~=3,   error('size(k,3)(=%g)~=3',size(k,3)); end
if max(abs(k(:)))>0.501
    warning('max(abs(k(:)))(=%g)>0.5',max(abs(k(:)))); 
end
if ~isa(k,'double')
    warning('k not double; converting'); 
    k = double(k); 
end
mtx = double(mtx);
numIter = double(numIter);
osf = double(osf);
verb = double(verb); 


%% reshape k
k = permute(k,[3 2 1]);


%% actual function call, dcf calculation
dcf = sdc3_MAT(k,numIter,mtx,verb,osf);
dcf = dcf.'*osf.^3*fufa_dcf;                 % correct scaling


%% checks
if any(isinf(dcf(:))), warning('isinf(dcf)'); dcf(isinf(dcf)) = 0; end
if any(isnan(dcf(:))), warning('isnan(dcf)'); dcf(isnan(dcf)) = 0; end
if any(dcf(:)<0),      warning('dcf<0');      dcf(dcf<0) = 0; end
if any(dcf(:)>thold)
    warning('dcf>%g; cropping',thold);
    dcf(dcf>thold) = thold;
end

end      % main function calc_dcf.m
