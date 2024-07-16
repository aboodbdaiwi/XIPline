function zte_b1shimming(fname,centre,radius,method,sv,verb,thresh_fac)
%ZTE_B1SHIMMING  Wrapper function for phase_b1map and b1optimiser
% zte_b1shimming(fname,centre,radius,method,sv,verb,thresh_fac)
%     fname  Filename for reconstructed data (*.mat)
%            or structure from "zte=load(fname);"
%            or directory name containing *.mat file
%    centre  Index to centre of circular mask     [pixel]
%            (default=(mtx/2+1)*[1,1,1])
%    radius  Radius of circular mask              [pixel]
%            (default=h.image.shim_fov) 
%    method  Method for b1optimiser           (default=2)
%        sv  Save results                     (default=false)
%      verb  Verbose (+plotting)              (default=1)
%thresh_fac  Factor for mask threshold        (default=0.5)
%
% See also B1MAP_PHASE, B1OPTIMISER.
% 10/2017 Rolf Schulte
if nargin<1, help(mfilename); return; end


%% defaults
if ~exist('centre','var'), centre = []; end
if ~exist('radius','var'), radius = []; end
if ~exist('method','var'), method = []; end
if isempty(method),        method = 2; end
if ~exist('sv','var'),     sv = []; end
if isempty(sv),            sv = false; end
if ~exist('verb','var'),   verb = []; end
if isempty(verb),          verb = 1; end
if ~exist('thresh_fac','var'),thresh_fac = []; end
if isempty(thresh_fac),       thresh_fac = 0.5; end


%% loading zte data files
if isstruct(fname)
    zte = fname; 
    clear fname
else
    if verb>0, fprintf('Loading zte data file %s\n',fname); end
    if exist(fname)~=2
        if isempty(fname), fname = '.'; end
        xx = dir([fname '/*.mat']);
        if isempty(xx), error('file %s not found',fname); end
        % fname = [xx(1).folder '/' xx(1).name];  % not existing in V8.1
        fname = [fname '/' xx(1).name];
        if length(xx)>1, warning('more than 1 file found; choosing %s',fname); end
        if ~exist(fname,'file'), error('file %s not found (bug)',fname); end
    end
    zte = load(fname);
end


%% filename for saving
if sv
    if ~exist('fname','var')
        % fname = [datestr(datetime('now'),'yyyymmdd') '_zte_b1shimming'];
        fname = [datestr(now,'yyyymmdd_hhMMss') '_zte_b1shimming'];
        warning('fname missing; generated=%s',fname);
    end
    if ~isempty(regexpi(fname,'\.mat$'))
        fname = fname(1:end-4);
    end
end


%% data checks
if ~isfield(zte,'h'), error('h missing'); end
if ~isfield(zte,'bb'), error('bb missing'); end
if zte.h.rdb_hdr.user10~=3,
    warning('zte.h.rdb_hdr.user10(=%g)~=3',zte.h.rdb_hdr.user10);
end
mtx = size(zte.bb);


%% generate cardiac mask: sphere around centre of shim volume
res = zte.h.image.slthick;                  % image resolution [mm]
if isempty(radius)
    fov = zte.h.image.shim_fov(1);          % shim FOV [mm]
    if fov>0
        radius = ceil(fov/res/2);
        if verb>0
            fprintf('res=%g [mm], fov=%g[mm],radius=%g[pixel]\n',...
                res,fov,radius);
        end
    else
        radius = 15;
        fprintf('fov=%g -> radius=%g\n',fov,radius);
    end
end
if isempty(centre)
    centre = ceil((mtx(1:3)/2+1));
end
% ctr = [hi.shim_ctr_R(1),hi.shim_ctr_A(1),hi.shim_ctr_S(1)]-...
%     [hi.ctr_R(1),hi.ctr_A(1),hi.ctr_S(1)];
% ind_cent = floor(ctr/res+mtx/2+0.5);
% ind_cent = floor(mtx/2);

mask = circular_mask(mtx(1:3),centre,repmat(radius,[1 3]));


%% if single Tx channel, do statistics
if size(zte.bb,4)==1
    fprintf('Single Tx data; showing data statistics only\n');
    bbabs = abs(sqrt(mean(conj(zte.bb).*zte.bb,5)));
    bbabs = bbabs/max(bbabs(:));
    bvec = bbabs(mask);
    fprintf('Signal over masked area: %.3g +- %.3g; STD=%.3g [%%]\n',...
        mean(bvec),std(bvec),std(bvec)/mean(bvec)*100);
    fprintf('Signal over whole image: %.3g +- %.3g; STD=%.3g [%%]\n',...
        mean(bbabs(:)),std(bbabs(:)),std(bbabs(:))/mean(bbabs(:))*100);
    figure(17); clf    
    imagesc_ind3d((bbabs.*mask + 0.25*bbabs.*~mask)*1.7,centre,[0 1]);
    
    title(sprintf('DD [%g, %g, %g]; STD=%.3g [%%]\n',...
        zte.h.rdb_hdr.dd_mode,...
        zte.h.ps.dd_q_ta_offset_from_qd,...
        zte.h.ps.dd_q_phase_delay_from_qd,...
        std(bvec)/mean(bvec)*100));
    
    if sv, print([fname '_b1statistics.png'],'-dpng'); end
    return;
end


%% relative b1mapping
[b1map,bbabs] = b1map_phase(zte.bb,thresh_fac,verb);


%% actual b1 shimming
if verb>1, figure(15); clf; end
[b1fit,txAtten,txPhase] = b1optimiser(b1map,mask,method,'',verb);


%% plotting
figure(16); clf;
bplt = zeros(mtx(1),mtx(2),mtx(3),3);
% bplt(:,:,:,1)   = (bbabs.*mask + 0.25*bbabs.*~mask)*1.7;
bplt(:,:,:,1)   = bbabs*1.7;
% bplt(:,:,:,2:3) = 2*abs(b1map);
bplt(:,:,:,2)   = (angle(b1map(:,:,:,2))+pi)/(2*pi);
bplt(:,:,:,3)   = abs(b1fit)/1.7;
imagesc_ind3d(bplt,centre,[0 1]);
plot_circle(centre(2),centre(1),radius,'w');
plot_circle(centre(3)+mtx(2),centre(1),radius,'w');
plot_circle(centre(3)+mtx(1)+mtx(2),centre(2),radius,'w');
plot_circle(centre(2),centre(1)+mtx(1),radius,'w');
plot_circle(centre(3)+mtx(2),centre(1)+mtx(1),radius,'w');
plot_circle(centre(3)+mtx(1)+mtx(2),centre(2)+mtx(1),radius,'w');


xlabel('xy                  xz                    yz');
% ylabel('b1fit  phase(Q-I)   b(Q)   b(I)   b(comb+mask)');
ylabel('b1fit      phase(Q-I)      bbabs');
title_str = sprintf('txAtten=%g[db/10]; txPhase=%g [deg]',...
    floor(txAtten+0.5),floor(txPhase+0.5));
title(title_str);


%% saving results
if sv
    print([fname '_b1shim.png'],'-dpng');
    b1fit = single(b1fit);
    save([fname '_b1shim.mat'],'-mat','-v7.3','txAtten','txPhase','method',...
        'b1map','b1fit','bbabs','centre','fov','mask','mtx','radius','res');
end



end      % main function (zte_b1shimming.m)
