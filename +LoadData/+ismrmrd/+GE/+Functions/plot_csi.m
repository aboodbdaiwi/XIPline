function plot_csi(spec,pc,scale,figstr,fname,method,h)
%PLOT_CSI  Plots CSI data in matrix form
%  plot_csi(spec,pc,scale,figstr,fname,method)
%                                                               (default)
%      spec  Reconstructed spectra [#s,#x,#y,#z,#t]
%        pc  Phase correction [deg]                             ([])
%            if emtpy -> magnitude
%     scale  scaling of spectra                                 ('max')
%                       'all'  global maximum + minimum
%                       'ind'  individually
%             [minVal maxVal]  global to these value (method==2)
%                      maxVal  global (method==1)
%    figstr   Figure title string                               ('')
%     fname   Print to file (if string)                         ([])
%    method   1 = plot into one axes                            (1)
%             2 = create multiple axes
%             3 = same as 1 except orientation same as PlotSquareCuboid
%             +10 = invisible -> save only w/o displaying figure
%         h   p-file header:                                    ([])
%             if given -> export plots as secondary-capture dicom
%
% 8/2020 Rolf Schulte
%See also RECON_CSI, RECON_MRSI.
if nargin<1, help(mfilename); return; end


%% misc input defaults + checks
if ~exist('pc','var'),     pc = []; end
if ~exist('scale','var'),  scale = []; end
if isempty(scale),         scale = 'all'; end
if ~exist('figstr','var'), figstr = []; end
if ~exist('fname','var'),  fname = []; end
if ~exist('method','var'), method = []; end
if isempty(method),        method = 1; end
if method>10
    nofig = true;
    method = method-10;
else
    nofig = false;
end
if ~exist('h','var'), h = []; end
if isempty(h)
    export_dcm = false;
else
    export_dcm = true; 
end

[ns,nx,ny,nz,nt] = size(spec);
if ns<32, warning('ns(=%g)<32: spectra along 1st dim?',ns); end
if (length(size(spec))>5), error('length(si)>5'); end
if ((nx==1)&&(ny==1)&&(nz==1))
    warning('>=1D spatial required: nx=%g,ny=%g,nz=%g',nx,ny,nz);
end


%% phasing of spectra (global)
if isreal(spec)
    if ~isempty(pc), warning('spec real + pc=%g; ignoring pc',pc); end
    if isempty(figstr), figstr = 'CSI: spec w/o pc'; end
else
    if isempty(pc)
        spec = abs(spec);
        if isempty(figstr), figstr = 'CSI: abs(spec)'; end
    else
        spec = real(spec*exp(-1i*pc/180*pi));
        if isempty(figstr) 
            figstr = sprintf('CSI: real(spec with pc=%g)',pc);
        end
    end
end


%% plotting
for l5=1:nt
    % misc prep
    ss = squeeze(spec(:,:,:,:,l5));
    if method==3, ss = permute(ss,[1 3 2 4]); end
    [n1,n2,n3,n4] = size(ss);

    % actual plotting
    for l4=1:n4
        if nofig
           fid = figure('Visible','off');
        else
            if ((n4>1)||(nt>1))
                fid = figure;
            else
                fid = gcf;
                clf;
            end
        end
        
        switch method
            case {1,3}      % plotting into single axes
                if isnumeric(scale)
                    if scale(1)<max(ss(:))/3
                        warning('scale(1)(=%g) < max(ss(:))(=%g)',...
                            scale(1),max(ss(:)));
                    end
                    ss = ss/scale(1);
                else
                    ss = ss/max(ss(:));             % normalise to 1
                    if regexpi(scale,'ind')
                        [nn1,nn2]=size(ss);
                        for ll=1:nn2, ss(:,ll) = 0.8*ss(:,ll)/max(ss(:,ll)); end
                    end
                end
                axes('units','norm','pos',[0 0 1 1])
                hold on
                for l2=(-n2/2:n2/2), plot([-1 1]*n3/2,l2*[1 1],'k'); end   % MRS grid
                for l3=(-n3/2:n3/2), plot(l3*[1 1],[-1 1]*n2/2,'k'); end   % MRS grid
                for l2=1:n2
                    % plot(linspace(-n3/2,n3/2,n1*n3),...
                    %     reshape(ss(n1:-1:1,l2,:,l4)-l2+n2/2,[n1*n3 1]),'b');
                    plot(linspace(-n3/2,n3/2,n1*n3),...
                        reshape(ss(:,l2,:,l4)-l2+n2/2,[n1*n3 1]),'b');
                end
                hold off
                axis([-n3/2 n3/2 -n2/2-0.3 n2/2]);
                
            case 2          % plotting into multiple axes
                if isnumeric(scale)
                    minVal = scale(1);
                    maxVal = scale(2);
                else
                    if regexpi(scale,'all')
                        minVal = min(ss(:));
                        maxVal = max(ss(:));
                    end
                end
                for l2=1:n2
                    for l3=1:n3
                        axes('units','norm','pos',[(l2-1)/n2 (l3-1)/n3 1/n2 1/n3]);
                        plot((1:n1),ss(:,l3,l2,l4),'b',[1 n1],[0 0],'k:');
                        set(gca,'XtickLabel','','YTickLabel','');
                        if exist('maxVal','var')
                            axis([1 n1 minVal maxVal]);
                        else
                            axis tight;
                        end
                    end
                end
            otherwise, error('method(=%g) unknown',method);
        end
        str = figstr;
        if n4>1, str = sprintf('%s; slice=%g',str,l4); end
        if nt>1, str = sprintf('%s; time step=%g',str,l5); end
        set(fid,'Name',str);
        drawnow;
        if ~isempty(fname)
            print(sprintf('%s_sl%g_ts%g',fname,l4,l5),'-dpng','-r600','-painters');
            % print(sprintf('%s_sl%g_ts%g',fname,l4,l5),'-dsvg');
            % saveas(fid,sprintf('%s_sl%g_ts%g',fname,l4,l5),'fig');
            if export_dcm
                inp.InstanceNumber = (l5-1)*n4+l4;
                write_scdicom(sprintf('%s_sl%g_ts%g.dcm',fname,l4,l5),fid,h,inp);
            end
        end
        if nofig, close(fid); end
    end    % for l4=1:n4
end        % for l5=1:nt

end        % main function plot_csi.m
