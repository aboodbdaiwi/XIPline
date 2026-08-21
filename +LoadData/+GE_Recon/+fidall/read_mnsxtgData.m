function [dd,h] = read_mnsxtgData(fname,nucleus)
%READ_MNSXTGDATA  Read raw MNS XTG data (prescan_config.cfg->debug->saved in logdir)
%     dd = read_mnsxtgData(fname)
%  fname   mnsxtgData.X or just input X                                 (1)
%nucleus   Nucleus (if h out)                                          (13)
%     dd   Raw FID signal (cell with mnsxtg repetitions)
%      h   Pseudo header for mns_prescan.m
%
%  9/2022 Rolf Schulte
if (nargin<1), help(mfilename); return; end

%% parameters
if ~exist('fname','var'),   fname = []; end
if isempty(fname),          fname = 1;  end
if isnumeric(fname),        fname = ['mnsxtgData.' num2str(fname)]; end
if ~exist('nucleus','var'), nucleus = []; end
if isempty(nucleus),        nucleus = 13; end
xres = 256;        % #FID points


%% read data
lcnt_old = -1;
fid = fopen(fname,'rb');

while true
    lcnt = fread(fid,1,'int');
    % fprintf('Reading mnsxtg_count=%d\n',lcnt); 
    if isempty(lcnt), fprintf('break\n'); break; end

    if lcnt~=lcnt_old
        d = complex(zeros(8,xres)); 
        lcnt_old = lcnt;
    end

    iframe = fread(fid,1,'int');
    xres = fread(fid,1,'int');
    num_rcv = fread(fid,1,'int');
    fprintf('mnsxtg_count=%d; iframe=%d; xres=%d; num_rcv=%d\n',...
        lcnt,iframe,xres,num_rcv);

    for lc=1:num_rcv
        ircv = fread(fid,1,'int');
        if (ircv+1)~=lc, warning('(ircv(=%d)+1)~=lc(=%d)',ircv,lc); end
        % fprintf('lc=%d ircv=%d\n',lc,ircv);
        tmp = fread(fid,2*xres,'double');
        if size(tmp,1)~=2*xres
            warning('size(tmp,1)(=%d)~=2*xres(=%d): filling with zeroes',...
                size(tmp,1),2*xres);
            tmp = [tmp ; zeros(2*xres-size(tmp,1),1)];
        end
        d(iframe+1,:,lc) = (tmp(1:2:end,1) + 1i*tmp(2:2:end,1)).';
    end
    dd{lcnt+1} = d;
    % fprintf('ftell(fid)=%g\n',ftell(fid));
    if feof(fid), break; end
end
tmp = fread(fid,'double');

fclose(fid);


%% fixing sign alterations
warning('bug in aps_getdata/ps_negate_data? every other sample negated; undoing');
fprintf('unchopping\n');
for lcnt=1:length(dd)
    dd{lcnt}(:,2:2:end,:,:) = -dd{lcnt}(:,2:2:end,:,:);
    dd{lcnt}(2:2:end,:,:,:) = -dd{lcnt}(2:2:end,:,:,:);
    fprintf('size(dd{%d}) = %s\n',lcnt,num2str(size(dd{lcnt})));
end


%% create pseudo header for mns_prescan.m
if nargout>1
    x = gyrogamma(1)/abs(gyrogamma(nucleus));
    h.rdb_hdr.ps_mps_freq = 1d4;
    h.image.slthick = 10;
    h.image.psd_iname = 'fidall';
    h.rdb_hdr.user32 = 0;              % #exc noise
    % B1max ref w/o pulse stretching [G]
    h.rdb_hdr.user34 = 12.8079*x/100;
    h.rdb_hdr.user0 = 2500;            % (full) BW [Hz]
    h.rdb_hdr.user2 = nucleus;
    h.image.user14 = 0;                % hardpulse excitation
    h.rdb_hdr.user19 = 2;              % # of blosi off-freq pulses/sign
    h.rdb_hdr.user30 = 2000;           % Bloch-Siegert frequency [Hz]
    h.rdb_hdr.user31 = 4000;           % pw blosi pulse [us] w/o stretching
    h.rdb_hdr.user33 = 0.0909*x;       % B1 of Bloch-Siegert pulse [G]
    h.rdb_hdr.ps_mps_tg = 0;           % current TG
    h.image.rawrunnum = 0; 
    h.exam.ex_no = str2double(fname(regexp(fname,'\d')));
    h.series.se_no = 0;
    h.rdb_hdr.data_collect_type = 1;   % chopping if odd
end



end      % read_mnsxtgData.m
