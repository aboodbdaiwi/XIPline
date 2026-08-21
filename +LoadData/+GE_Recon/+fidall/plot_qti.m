function plot_qti(bb,dim,plt,fname,figstr,cmax,fidoff,rot,cmstyle)
%PLOT_QTI Plot MR Parameter maps
% plot_qti(bb,dim,plt,fname,figstr,cmax,fidoff,rot,cmstyle)
%       bb   Complex, coil-combined svd-compressed images  (#x,#y,#z,#svds)
%      dim   Spatial dimensions                                         (3)
%      plt   Visible plotting                                        (true)
%    fname   Save plots to fname                                       ('')
%   figstr   String for figure name                                    ('')
%     cmax   Max values for colormap (T1,T2,PD)      [ms] ([3000,300,1000])
%   fidoff   Figure ID offset                                         (200)
%      rot   Rotate 3D slices via rot90 [x,y,z]                        ([])
%  cmstyle   Colormap (T1,T2,PD)                                        (1)
%            0=orig='copper','hot','gray'
%            1=ismrm='lipari','navia','gray'
%
%  9/2024 Rolf Schulte
if nargin<1, help(mfilename); return; end


%% default input parameters
if ischar(bb), load(bb); end
if ~exist('dim','var'),     dim = []; end
if isempty(dim)
    if size(bb,3)<17, dim=2;
    else,             dim=3;
    end
end
if ((dim<2) || (dim>3)),    warning('dim(=%d) not 2D or 3D',dim); end
if ~exist('plt','var'),     plt = []; end
if isempty(plt),            plt = true; end
if length(plt)~=1,          warning('length(plt)~=1'); end
if ~exist('fname','var'),   fname = ''; end
if ~exist('figstr','var'),  figstr = ''; end
if ~exist('cmax','var'),    cmax = []; end
if isempty(cmax),           cmax = [3000,300,1000]; end
if length(cmax)~=3,         error('length(cmax)(=%d)~=3',length(cmax)); end
if ~exist('fidoff','var'),  fidoff = []; end
if isempty(fidoff),         fidoff = 200; end
if ~exist('rot','var'),     rot = []; end
if ~exist('cmstyle','var'), cmstyle = []; end
if isempty(cmstyle),        cmstyle = 1; end
if size(bb,4)~=3, warning('size(bb,4)(=%g)~=3',size(bb,4)); end


%% plotting
if plt || ~isempty(fname)
    fprintf('Plotting results\n');
    map_str = {'T1','T2','PD'};
%     switch cmstyle
%         case 0, 
            cmap = {'copper','hot','gray'};
%         case 1, cmap = {'lipari','navia','gray'};
%         otherwise, error('cmstype (=%g) ~=0 or 1',cmstyle);
%     end
    for ll=1:3
        if plt
            fid = figure(ll+fidoff); clf;
        else   
            fid = figure('Visible','off'); 
        end
        set(fid,'DefaultAxesFontSize',14);
        if dim==3
            imagesc_ind3d(1d3*abs(bb(:,:,:,ll)),[],[],[],[],[],[],rot,true);
        else
            imagesc_row(1d3*abs(bb(:,:,:,ll)),'','',size(bb,3)>5,true,[],[],true);
        end
        caxis([0 cmax(ll)]); colormap(cmap{ll}); 
        if (ll==1) && (cmstyle==1), colormap(lipari([],0,cmax(ll))); end
        if (ll==2) && (cmstyle==1), colormap(navia([],0,cmax(ll))); end

        % colorbar
        set(fid,'name',sprintf('%s%smap',figstr,map_str{ll}));
        if ~isempty(fname)
            ffn = [fname '_' map_str{ll} 'map.png'];
            if isdeployed
                print(ffn,'-dpng','-r300','-painters');
            else
                print(ffn,'-dpng','-r300');
            end
        end
        if ~plt, close(fid); end
    end
end


end      % plot_qti.m
