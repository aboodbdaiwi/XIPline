function w = lambertw(k,x)
%LAMBERTW Lambert's W function.
%   W = LAMBERTW(Z) solves W*exp(W) = Z.
%   W = LAMBERTW(K,Z) is the K-th branch of this multi-valued function.
%
%   References:
%   [1] Robert M. Corless, G. H. Gonnet, D. E. G. Hare,
%   D. J. Jeffrey, and D. E. Knuth, "On the Lambert W Function",
%   Advances in Computational Mathematics, volume 5, 1996, pp. 329-359.
%
%   [2] Corless, Jeffrey, Knuth, "A Sequence of Series
%   for The Lambert W Function", ISSAC '97


%   More information available from:
%   http://www.apmaths.uwo.ca/~rcorless/frames/PAPERS/LambertW

%   Copyright 2009-2018 The MathWorks, Inc.
narginchk(1,2);

if nargin == 1
  x = k;
  k = 0;
end

if ~isscalar(k) && ~isscalar(x) && ~isequal(size(x),size(k))
    error(message('symbolic:sym:ScalarOrSameSize'));
end

if ~isnumeric(k) || ~isnumeric(x)
    error(message('symbolic:symbolic:NonNumericParameter'));
end

if any(~((isfinite(k) & k == fix(k)) | isnan(k)))
    error(message('symbolic:sym:lambertw:ExpectingInteger1'));
end

if isempty(k)
    w = k;
    return;
end
if isempty(x)
    w = x;
    return;
end

if isscalar(k)
    w0 = getStartingPoint(k,x);
else
    scale = ones(size(k)).*ones(size(x));
    k = scale .* k;
    x = scale .* x;
    w0 = arrayfun(@getStartingPoint,k,x);
end
w = w0;

% Haley's method
refine = true(size(w));
target = 1e-7;
if isa(x,'double')
    target = 1e-14;
end
loopCnt = 0;
while any(refine(:)) && (loopCnt < 25)
    e = exp(w);
    f = w.*e - x;  % Iterate to make this quantity zero
    w = w - f./((e.*(w+1) - (w+2).*f./(2*w+2)));
    refine = abs(f)./abs(x) > target;
    loopCnt = loopCnt + 1;
end

w(isnan(w)) = NaN; % do not return NaN+NaN*1i
w(isinf(x) & x > 0 & k==0) = inf;
w(x==0) = -Inf; % for consistency with sym.lambwertw
w(x==0 & k==0) = 0;
end

function w = getStartingPoint(k,x)
    switch(k)
    case 0
        w = getStartingPointBranch0(x);
    case -1
        w = getStartingPointBranchM1(x);
    otherwise
        w = getStartingPointBranchNonSpecialCase(k, x);
    end
end

function w = getStartingPointBranch0(z)
    w = nan(size(z),'like',z);
    % principal branch, series around zero:
    % w = z - z^2 + 3/2*z^3 - 8/3 * z^4 + O(z^5)
    startFound = abs(z) < 0.3;
    w(startFound) = z(startFound) - z(startFound).^2 + 3/2*z(startFound).^3;

    % close to -1/e, we need a special approximation
    closeToBranch = abs(z+0.3678794412) < 0.25;
    p = sqrt(2*(exp(1)*z(closeToBranch)+1));
    closeToBranch(closeToBranch) = abs(p) < sqrt(2);
    startFound = startFound | closeToBranch;
    w(closeToBranch) = -1 + p - 1/3*p.^2 + 11/72*p.^3;

    % close to 1, a constant is a good starting point:
    closeToBranch = abs(z-1) < 1e-3;
    startFound = startFound | closeToBranch;
    w(closeToBranch) = 0.56714329040;

    % the ln iteration in the generic branch easily drifts into the -1 branch
    % Use (28) from [1] instead, if possible
    closeToBranch = ~startFound & abs(log(z)-1) < sqrt(4+pi^2);
    startFound = startFound | closeToBranch;

    logZ1 = log(z(closeToBranch))-1;
    w(closeToBranch) = 1 + logZ1/2 + logZ1.^2/16 - logZ1.^3/192;

    % else, use generic approximation
    closeToBranch = ~startFound;
    if any(closeToBranch, 'all')
        w(closeToBranch) = getStartingPointBranchNonSpecialCase(0, z(closeToBranch));
    end
end

function w = getStartingPointBranchM1(z)
    w = nan(size(z),'like',z);
    % close to -1/e, we need a special approximation
    closeToBranch = imag(z) >= 0 & abs(z+0.3678794412) < 0.25;
    p = -sqrt(2*(exp(1)*z(closeToBranch)+1));
    closeToBranch(closeToBranch) = abs(p) < sqrt(2);
    p = p(abs(p) < sqrt(2));
    w(closeToBranch) = -1 + p - 1/3*p.^2 + 11/72*p.^3;

    w(~closeToBranch) = getStartingPointBranchNonSpecialCase(-1, z(~closeToBranch));
end

function w = getStartingPointBranchNonSpecialCase(k, z)
    v = log(z)+2i*k*pi;
    if v == 0 % can we get here for lambertw(0,1)?
        w = 0.5671;
        return
    end
    % [2], Section 4.4
    p = -log(v);
    w = v + v./(1+v).*p + 1/2.*v./(1+v).^3.*p.^2 - 1/6*v.*(2*v-1)./(1+v).^5.*p.^3;
end
