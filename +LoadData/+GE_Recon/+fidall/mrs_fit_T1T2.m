function mrs_fit_T1T2(d,h,what,tt,lb,int_width,int_roi)
%MRS_FIT_T1T2  Fit T1 and T2 relaxation times from inversion recovery and 
% spin echo spectra
%  mrs_fit_T1T2(d,h,what,tt,lb,int_width,int_roi)
%        d   Data (in cell structure)
%        h   Header (in cell structure)
%     what   1=fit T1, 2=fit T2                                         (2)
%       tt   T1->inversion times; T2->echo times      [ms]
%            T2 times stored in header
%       lb   Apodisation                              [Hz]         ([0 10])
%int_width   Integration width                        [Hz]            (500)
%  int_roi   Spectral region of interest              [Hz]     ([-700 700])
%            (for integration around max)
%
% 11/2023 Rolf Schulte
if (nargin<1), help(mfilename); return; end


%% default + checks
if ~iscell(d), error('~iscell(d)'); end
if ~exist('what','var'),     what = []; end
if isempty(what),            what = 2; end
if ~exist('tt','var'),       tt = []; end
if ~exist('lb','var'),       lb = []; end
if isempty(lb),              lb = [0 10]; end
if ~exist('int_width','var'), int_width = []; end
if isempty(int_width),       int_width = 500; end
if length(int_width)~=1,     error('length(int_width)~=1'); end
if ~exist('int_roi','var'),  int_roi = []; end
if isempty(int_roi),         int_roi = [-700 700]; end
if length(int_roi)~=2,       error('length(int_roi)~=2'); end


%% reading data
if ~isnumeric(d{1})
    fnames = d;
    clear d
    for ll=1:length(fnames)
        fprintf('Reading data %s\n',fnames{ll});
        [d{ll},h{ll}] = read_p(fnames{ll});
    end
end


%% reading times
if isempty(tt)
    switch what
        case 1, error('inversion time not stored in header');
        case 2
            for ll=1:length(d)
                tt(ll) = h{ll}.image.te*1d-3;
            end
        otherwise, error('what(=%d)~=1 or 2',what);
    end
end
if min(tt)<10, warning('min(tt)=%g<10: ms?',min(tt)); end
nt = length(d);


%% sorting data
if length(tt)~=nt, error('legnth(tt)(=%d)~=length(d)(=%d)',length(tt),nt); end
[tt,ii] = sort(tt(:));
for ll=1:nt
    tmpd{ll} = d{ii(ll)};
    tmph{ll} = h{ii(ll)};
end
d = tmpd; h = tmph;
clear tmpd tmph


%% checking header
R1 = zeros(1,nt); R2 = R1; TG = R1;
for ll=1:nt
    R1(1,ll) = h{ll}.rdb_hdr.ps_mps_r1;
    R2(1,ll) = h{ll}.rdb_hdr.ps_mps_r2;
    TG(1,ll) = h{ll}.rdb_hdr.ps_mps_tg;
end
if any(diff(R1)), warning('R1 = %s',num2str(R1)); end
if any(diff(R2)), warning('R2 = %s',num2str(R2)); end
if any(diff(TG)), warning('TG = %s',num2str(TG)); end



%% extracting signal b
for ll=1:nt
    [spec{ll},~,~,hz] = fid2spec(d{ll},h{ll},[],lb,[],[],false);
end
ns = size(spec{1},2);
ss = complex(zeros(nt,ns));
for ll=1:nt, ss(ll,:) = sum(spec{ll},1); end
if ~isempty(int_roi)
    ind_roi = (hz>int_roi(1)) & (hz<int_roi(2)); 
else
    ind_roi = true(1,ns);
end
[~,imax] = max(abs(ss(1,:).*double(ind_roi)));
fmax = hz(1,imax);
ind = ((hz-fmax)>-int_width/2) & ((hz-fmax)<int_width/2);
sig = zeros(1,nt);
for ll=1:nt
    sig(1,ll) = sum(ss(ll,ind),2);
end


%% normalisation of signal - real/abs
switch what
    case 1    % fit T1
        if abs(sig(1,1))> abs(sig(1,end))
            smax = -sig(1,1);
        else
            smax = sig(1,end);
        end
        sig = real(sig./smax);
    case 2    % fit T2
        sig = abs(sig./sig(1));
    otherwise, error('what(=%d)~=1 or 2',what);
end


sig = sig/max(abs(sig));


%% non-linear LSQ
x0 = [1,1000];        % scaling,T1/T2
lb = [0 0];
ub = [];
opt = optimoptions('lsqnonlin','Display','iter');
x = lsqnonlin(@(x)sub_decay_curve(x,what,tt.',sig),x0,lb,ub,opt);
fprintf('non-linear LSQ: T%d = %g [ms]\n',what,x(2));


%% fitted curve for plotting
t_fit = linspace(0,max(tt(~isinf(tt))),1000);
if what==1
    sig_fit = x(1)*(1-2*exp(-t_fit/x(2)));
else
    sig_fit = x(1)*exp(-t_fit/x(2));
end


%% plotting
clf; subplot(2,1,1);
[n1,n2] = size(ss);
xplt = repmat(hz,[n1 1]).';
zplt = repmat((1:n1).',[1 n2]).';
plot3(zplt,xplt,abs(ss),'r:',zplt(ind,:),xplt(ind,:),abs(ss(:,ind)),'g')
axis tight
view(-80,24);
ylabel('freq [Hz]');
grid on
title('Spectra');

subplot(2,1,2);
plot(tt,real(sig),'b-x',t_fit,sig_fit,'g');
grid on
if what==1
    xlabel('inversion time [ms]');
    title(sprintf('T1 = %g [ms]',round(x(2))));
    legend({'signal','fit'},'location','southeast')
    axis([0 max(t_fit) -1.1*x(1) 1.1*x(1)]);
else
    xlabel('TE [ms]');
    title(sprintf('T2 = %g [ms]',round(x(2))));
    legend({'signal','fit'})
    axis([0 max(t_fit) -0.1*x(1) 1.1*x(1)]);
end



end      % mrs_fit_T1T2.m


%% sub-function for non-linear LSQ
function ff = sub_decay_curve(x,what,tt,sig)
    switch what
        case 1, ff = x(1)*(1-2*exp(-tt/x(2))) - sig;
        case 2, ff = x(1)*exp(-tt/x(2)) - sig;
    end
end      % sub_decay_curve

