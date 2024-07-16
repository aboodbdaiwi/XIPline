function characterise_grad_noise(d,h,sv)
%CHARACTERISE_GRAD_NOISE  Quantify gradient noise measured with spirals
%   on different axes: none - Gx - Gy - Gz - all
%
% characterise_grad_noise(d,h,sv)
%   d  raw P-file data
%      optionally filename to P-file or image archive
%   h  header structure
%  sv  Save results; can be file name; if logical, create file name
%
% 11/2018 Rolf Schulte
if nargin<1, help(mfilename); return; end

if ~exist('sv','var'),  sv = []; end
if isempty(sv),         sv = false; end

%% misc pars
nax = 5;                                    % 5 different measurements
nrep = 4;                                   % each repeated 4 times
ax_str = {'G off','Gx','Gy','Gz','G all'};  % axes label


%% reading in data, if pfile name is given
if ~isnumeric(d)
    if exist(d,'file'),
        fname = d;
        [d,h] = read_p(d);
    else
        warning('strange input d');
    end
end

%% checks
bw = h.rdb_hdr.user0;                  % full acquisition bandwidth [Hz] 
if h.image.user4~=(nax*nrep)
    % #excitations
    warning('h.image.user4(=%g)~=(nax(=%g)*nrep(=%g))',...
        h.image.user4,nax,nrep);
end
if h.image.user3~=27
    % gradient trajectory
    warning('h.image.user3(=%g)~=27',h.image.user3);
end
if h.image.plane~=2
    % Scan Plane (different from opplane)
    % 2=axial, 4=sagittal, 8=coronal, 16=oblique
    warning('h.image.plane(=%g)~=1; ylabel wrong',h.image.plane);
end
if h.rdb_hdr.position~=1
    % Patient Position: 1=Supine, 2=Prone, 3=Left Decub, 4= Right Decub
    warning('h.rdb_hdr.position(=%g)~=1; ylabel wrong',h.rdb_hdr.position);    
end
if h.image.freq_dir~=1
    warning('h.image.freq_dir(=%g)~=1; ylabel wrong',h.image.freq_dir);
end


%% filenames
if ~exist('fname','var'), fname = sprintf('P%.5d.7',h.image.rawrunnum); end
tmp = regexpi(fname,'/');         % if directory, take name only
if ~isempty(tmp), fname = fname((tmp(end)+1):end); end
tmp = regexpi(fname,'\');         % same for win dir separators
if ~isempty(tmp), fname = fname((tmp(end)+1):end); end

if ischar(sv)
    % if filename given in sv field
    plt_fname = sv;
    if ~isempty(regexpi(plt_fname,'\.7$')),  plt_fname = plt_fname(1:end-2); end
    if ~isempty(regexpi(plt_fname,'\.h5$')), plt_fname = plt_fname(1:end-3); end
    sv = true;
else
    % remove suffix
    if sv
        tmp = regexpi(fname,'\.');
        if ~isempty(tmp)
            plt_fname = fname(1:(tmp(end)-1));
        else
            plt_fname = fname;
        end
    end
end

%% check + reshape data
[n1,n2,n3,n4,n5,n6] = size(d);
if n1~=(nax*nrep)
    warning('size(d,1)(=%g)~=(nax(=%g)*nrep(=%g))',n1,nax,nrep);
end
if n3~=1, warning('size(d,3)(=%g)~=1',n3); end
if n4~=1, warning('size(d,4)(=%g)~=1',n4); end
if n5~=1, warning('size(d,5)(=%g)~=1',n5); end
dd = zeros(nrep*n6,n2,nax);       % reshape: put nrep+ncoils into dim1
for lax=1:nax
    for lrep=1:nrep
        for lcoil=1:n6
            % fprintf('lax=%g; lrep=%g; lcoil=%g: dd ind1=%g; d ind1=%g\n',...
            %     lax,lrep,lcoil,lrep+nrep*(lcoil-1),lrep+nrep*(lax-1));
            dd(lrep+nrep*(lcoil-1),:,lax) = ...
                d(lrep+nrep*(lax-1),:,1,1,1,lcoil);
        end
    end
end


%% noise statistics
fd = ifftshift(ifft(dd,[],2),2);                 % freq domain
fd = abs(fd);
dd = abs(dd);
dd1 = reshape(dd,[size(dd,1)*size(dd,2),size(dd,3)]);
fd1 = reshape(fd,[size(dd,1)*size(dd,2),size(dd,3)]);

mean_td = mean(dd1,1);
std_td = std(dd1,[],1);
mean_fd = mean(fd1,1);
std_fd = std(fd1,[],1);
outp = sprintf('\t\ttime domain\t\t\tfrequency domain\n');
for lax=1:nax
    outp = sprintf('%s%5s:\t%.4g +- %.4g\t',...
        outp,ax_str{lax},mean_td(1,lax),std_td(1,lax));
    outp = sprintf('%s%.4g +- %.4g\n',outp,mean_fd(1,lax),std_fd(1,lax));
end
fprintf(outp);

%% plotting
tt = (0:(n2-1))/bw*1d3;                          % time axis
ax_td = [1 tt(end) 0 max(abs(dd(:)))];           % same td scaling
ff = (-n2/2:(n2/2-1))/n2*bw/1d3;                 % freq axis
ax_fd = [ff(1) ff(end) 0 max(abs(fd(:)))];       % same fd scaling

clf
posf = get(gcf,'position');
ss = get(0,'ScreenSize');
if (posf(2)+780)>ss(4), posf(2) = ss(4)-780; end
set(gcf,'position',[posf(1:2) 600 700]);

for lax=1:nax
    % time domain
    subplot(nax,2,lax*2-1);
    plot(tt,abs(dd(:,:,lax)).');
    grid on; axis(ax_td);
    ylabel(ax_str{lax});
    posa = [0.05 (0.8-0.9*(lax-1)/nax) 0.45 0.8/nax];
    set(gca,'Position',posa);
    if lax==1, title('time domain'); end
    if lax==nax
        xlabel('time [ms]');
    else
        set(gca,'XTickLabel','');
    end
    set(gca,'YTickLabel','');
    text(tt(ceil(0.51*n2)),0.92*ax_td(4),...
        sprintf('%.3g+-%.3g',mean_td(1,lax),std_td(1,lax)));
    drawnow

    % frequency domain
    subplot(nax,2,lax*2);
    plot(ff,abs(fd(:,:,lax)).');
    grid on; axis(ax_fd);
    posa(1) = 0.52;
    set(gca,'Position',posa);
    if lax==1, title('frequency domain'); end
    if lax==nax
        xlabel('frequency [kHz]');
    else
        set(gca,'XTickLabel','');
    end
    set(gca,'YTickLabel','');
    text(ff(ceil(0.51*n2)),0.92*ax_fd(4),...
        sprintf('%.3g+-%.3g',mean_fd(1,lax),std_fd(1,lax)));
    drawnow

    
end

figstr = sprintf('P%05d Exam%d Series%d',...
    h.image.rawrunnum,h.exam.ex_no,h.series.se_no);
set(gcf,'name',figstr);

if sv
    print('-dpng',[plt_fname '.png'],'-r600');
    saveas(gcf,[plt_fname '.fig'],'fig');
    fid = fopen([plt_fname '.txt'],'w');
    fprintf(fid,outp);
    fclose(fid);
end

end   % main function charactrise_grad_noise.m
