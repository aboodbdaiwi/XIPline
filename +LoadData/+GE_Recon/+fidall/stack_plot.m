function stack_plot(d,h,lb,how,sv)
%STACK_PLOT  Display spectra as stacked plot
%  stack_plot(d,h,lb,how,sv)
%   d  raw P-file data
%      optionally filename to P-file or image archive
%   h  header structure
%  lb  line broadening (exponential,Gaussian)   [Hz]  []
% how  Plotting style (1=plot, 2=plot3)               (2)
%  sv  Print pix and save data to file   [logical]    (false)
%      alternatively: filename for printing
%
% 12/2023 Rolf Schulte
if nargin<1, help(mfilename); return; end


%% input parameters
if ~exist('lb','var'),  lb = []; end
if ~exist('how','var'), how = []; end
if ~exist('sv','var'),  sv = []; end
if isempty(sv),         sv = false; end


%% reading in data, if pfile name is given
if ~isnumeric(d) || isempty(d)
    if ~exist(d,'file') && ~isempty(d)
        warning('strange input for d (=''%s''); file not existing?',d); 
    end
    [d,h] = read_p(d,true);
end


%% filenames
if isempty(h), error('isempty(h)'); end
if ~exist('fname','var'), fname = sprintf('P%.5d.7',h.image.rawrunnum); end
tmp = regexpi(fname,'/');         % if directory, take name only
if ~isempty(tmp), fname = fname((tmp(end)+1):end); end
tmp = regexpi(fname,'\');         % same for win dir separators
if ~isempty(tmp), fname = fname((tmp(end)+1):end); end

if ischar(sv)
    % if filename given in sv field
    plt_fname = sv;
    if ~isempty(regexpi(plt_fname,'\.7$', 'once')),  plt_fname = plt_fname(1:end-2); end
    if ~isempty(regexpi(plt_fname,'\.h5$', 'once')), plt_fname = plt_fname(1:end-3); end
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


%% reconstruction
[spec,par,fid,hz]=fid2spec(d,h,[],lb);
if sv, save([plt_fname '.mat'],'-mat','spec','par','fid','hz','h'); end


%% plotting
[n1,n2,n3]=size(spec);
if isempty(how)
    if n1==1
        how = 1;
    else
        how = 2;
    end
end
for l3=1:n3
    if n3>1
        figure(100+l3); clf
    else
        clf
    end
    abspec = abs(spec(:,:,l3)).';
    switch how
        case 1
            xx = repmat(linspace(0,max(abs(spec(:))),n1),[n2 1]);
            plot(hz,abspec+xx);
            axis tight;
            xlabel('freq [Hz]');
            set(gca,'XDir','reverse');
        case 2
            x = repmat(hz,[n1 1]).';
            z = repmat((1:n1).',[1 n2]).';
            plot3(z,x,abspec);
            axis tight
            view(-85,23);
            ylabel('freq [Hz]');
            % zlabel('signal [a.u.]')
        otherwise, error('how(=%g)~=1 or 2',how);
    end
    t_str = sprintf('%s  n1=%g  n2=%g',regexprep(fname,'_','\\_'),n1,n2);
    if n3>1, t_str = sprintf('%s   l3=%g',t_str,l3); end
    title(t_str);
    figstr = sprintf('P%05d Exam%d Series%d',...
        h.image.rawrunnum,h.exam.ex_no,h.series.se_no);
    if n3>1, figstr = sprintf('%s Slice%d',figstr,l3); end
    set(gcf,'name',figstr);
    grid on
    
    if sv
        if isdeployed
            renderer = '-painters';
        else
            renderer = '-opengl';
        end
        l3str = '';
        if n3>1, l3str = ['_l' num2str(l3)]; end
        print('-dpng',[plt_fname l3str '.png'],'-r300',renderer);
        saveas(gcf,[plt_fname l3str '.fig'],'fig');
        if true
            inp.InstanceNumber = l3;
            write_scdicom([plt_fname l3str '.dcm'],gcf,h,inp);
        end
    end
end


end      % stack_plot.m
