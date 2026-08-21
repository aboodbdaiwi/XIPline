function [T1,T2] = calc_T1T2_qti(t1val,t2val,n1,n2,cutfac,how,verb)
%CALC_T1T2_QTI Calculate T1 and T2 distribution for QTI
% [T1,T2] = calc_T1T2_qti(t1val,t2val,n1,n2,cutfac,how,verb)
%   t1val   T1 range [min,mid,max]                           ([0.1,1.0,5])
%   t2val   T2 range [min,mid,max]                          ([0.01,0.2,2])
%      n1   #T1s                                                  ([80 80])
%      n2   #T2s                                                  ([80 80])
%  cutfac   Keep: (T1>cutfac(1)*T2) & (T1<cutfac(2)*T2)          ([1.5,15])
%     how   Method to distribute points                                 (3)
%    verb   1=verbos, 2+=plotting                                       (1)
%
%      T1   Vector with T1s
%      T2   Vector with T2s
%
%  9/2022 Rolf Schulte
if ((nargin<1)&&(nargout<1)), help(mfilename); return; end


%% default values
if ~exist('t1val','var'),  t1val = []; end
if isempty(t1val),         t1val = [0.1,1,5]; end
if ~exist('t2val','var'),  t2val = []; end
if isempty(t2val),         t2val = [0.01,0.2,2]; end
if ~exist('n1','var'),     n1 = []; end
if isempty(n1),            n1 = [80 80]; end
if ~exist('n2','var'),     n2 = []; end
if isempty(n2),            n2 = [80 80]; end
if ~exist('cutfac','var'), cutfac = []; end
if isempty(cutfac),        cutfac = [1.5 15]; end
if ~exist('how','var'),    how = []; end
if isempty(how),           how = 3; end
if ~exist('verb','var'),   verb = []; end
if isempty(verb),          verb = 1; end


%% print info
if verb>0
    fprintf('T1 range: min=%g, mid=%g, max=%g\n',t1val(1),t1val(2),t1val(3));
    fprintf('T2 range: min=%g, mid=%g, max=%g\n',t2val(1),t2val(2),t2val(3));
    fprintf('n1=%g %g, n2=%g %g, cutfac=%g %g, how=%g\n',...
        n1(1),n1(2),n2(1),n2(2),cutfac(1),cutfac(2),how);
end


%% create T1/T2 vectors
switch how
    case 1
        tt1 = linspace(sqrt(t1val(1)),sqrt(t1val(2)),n1(1)).^2;
        tt2 = linspace(sqrt(t2val(1)),sqrt(t2val(2)),n2(1)).^2;
        T1 = repelem(tt1,1,length(tt2)); 
        T2 = repmat(tt2,1,length(tt1));
    case 2
        tt1 = sub_sampling(t1val,n1(1),0.05);
        tt2 = sub_sampling(t2val,n2(1),0.04);
        T1 = repelem(tt1,1,length(tt2)); 
        T2 = repmat(tt2,1,length(tt1));
    case 3
        tt1a = linspace(t1val(1),t1val(2),n1(1));
        tt2a = linspace(t2val(1),t2val(2),n2(1));
        T1a = repelem(tt1a,1,length(tt2a));
        T2a = repmat(tt2a,1,length(tt1a));
        tt1b = linspace(sqrt(t1val(1)),sqrt(t1val(3)),n1(2)).^2;
        tt2b = linspace(sqrt(t2val(1)),sqrt(t2val(3)),n2(2)).^2;
        T1b = repelem(tt1b,1,length(tt2b)); 
        T2b = repmat(tt2b,1,length(tt1b));
        ii = (T1b<t1val(2)) & (T2b<t2val(2));
        T1b = T1b(1,~ii);
        T2b = T2b(1,~ii);
        T1 = [T1a,T1b];
        T2 = [T2a,T2b];
    case 4
        tt1 = [linspace(t1val(1),t1val(2),n1(1)),...
            linspace(t1val(2),t1val(3),n1(2))];
        tt2 = [linspace(t2val(1),t2val(2),n1(1)),...
            linspace(t2val(2),t2val(3),n1(2))];
        T1 = repelem(tt1,1,length(tt2)); 
        T2 = repmat(tt2,1,length(tt1));
    otherwise, error('how (=%g) unkown');
end


%% keep only T1/T2 pairs in a certain range
if ~islogical(cutfac)
    ii = (T1>cutfac(1)*T2) & (T1<cutfac(2)*T2);
    T1 = T1(1,ii);
    T2 = T2(1,ii);
    if verb>0
        fprintf('Keeping %d/%d of values with (T1>%g*T2) & (T1<%g*T2)\n',...
            sum(ii),length(ii),cutfac(1),cutfac(2));
    end
end
if verb>1, plot(T1,T2,'x'); grid on; end


end      % calcT1T2_qti.m


%% sub-functions
function tt = sub_sampling(tval,nt,fufa)
tmin = tval(1);
tmain = tval(2);
tmax = tval(3);

g = (tmain-tmin)/(2*(tmax-tmain));
% n0 = round(1/(2*(1+g)) + sqrt(1/(4*(1+g)^2) + g/(1+g)*nt^2))
n0 = round((sqrt(g/(1+g))-fufa)*nt);

tmp1 = linspace(tmin,tmain,n0+1);
tmp2 = linspace(sqrt(tmain),sqrt(tmax),nt-n0).^2;
tt = [tmp1(1,1:n0),tmp2];

% subplot(2,1,1); plot(tt,'x-b')
% subplot(2,1,2); plot(diff(tt),'x-b')


end      % sub_sampling

