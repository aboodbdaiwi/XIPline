function spec = csi_orientation(spec,rot90fac,trnsps,do_voxshift)
%CSI_ORIENTATION  Orient images (xy plane) from xyz to AP-RL-SI coordinates
%       spec = csi_orientation(spec,rot90fac,trnsps,do_voxshift)
%       spec   Spectrum        [#s,#x,#y,#z,#t,#c]
%   rot90fac   h.data_acq_tab.rotate
%     trnsps   h.data_acq_tab.transpose
%do_voxshift   Shift back by one voxel depending on orientation     (false)
%
%  5/2024 Rolf Schulte
if nargin<1, help(mfilename); return; end


%% misc parameters + checks
if ~exist('do_voxshift','var'), do_voxshift = []; end
if isempty(do_voxshift),        do_voxshift = false; end

if length(size(spec))>6, warning('spec(=%dD)>6D',length(size(spec))); end
% rot90fac = h.data_acq_tab.rotate;
% trnsps =   h.data_acq_tab.transpose;
if any(diff(rot90fac)),  warning('rot90fac'); disp(rot90fac); end
if any(diff(trnsps)),    warning('trnsps');   disp(tsnsps); end
rot90fac = rot90fac(1);
trnsps = trnsps(1);
fprintf('csi_orientation: rot90fac=%g trnsps=%g\n',rot90fac,trnsps);
if ~any(trnsps==[0,3]), warning('trnsps(=%g)~=[0,3]',trnsps); end
if ((rot90fac<0) || (rot90fac>3)), warning('rot90fac(=%g)<0 or >3',rot90fac); end 


%% orientation
spec = spec(:,end:-1:1,end:-1:1,end:-1:1,:,:);

spec = permute(spec,[2,3,4,5,6,1]);
mtx = size(spec);
mtx = [mtx , ones(1,6-length(mtx))];
if trnsps(1)==3
    if mod(rot90fac(1),2)<1
        mtx = mtx(1,[2,1,3,4,5,6]);
    end
else
    if mod(rot90fac(1),2)==1
        mtx = mtx(1,[2,1,3,4,5,6]);
    end
end
ss = complex(zeros(mtx,class(spec)));
for l3=1:size(spec,3)
    for l4=1:size(spec,4)
        for l5=1:size(spec,5)
            for l6=1:size(spec,6)
                % image orientation
                if trnsps(1)==3
                    ss(:,:,l3,l4,l5,l6) = rot90(fliplr(spec(:,:,l3,l4,l5,l6)).',rot90fac(1));
                else
                    ss(:,:,l3,l4,l5,l6) = rot90(fliplr(spec(:,:,l3,l4,l5,l6)),rot90fac(1));
                end
            end
        end
    end
end


%% correct one voxel shift
% test via:
% ss=zeros(1,4,4);ss(1,3,3)=1;
% ss2=csi_orientation(ss,0,0,true);squeeze(ss2)
if do_voxshift
    switch trnsps
        case 0
            switch rot90fac
                case 0, ss = ss([end,1:end-1],:,:,:,:,:);
                case 1, ss = ss([end,1:end-1],[end,1:end-1],:,:,:,:);
                case 2, ss = ss(:,[end,1:end-1],:,:,:,:);
                case 3
                otherwise, warning('rot90fac(=%g)~=0,1,2,3',rot90fac);
            end
        case 3
            switch rot90fac
                case 0, ss = ss(:,[end,1:end-1],:,:,:,:);
                case 1, ss = ss(:,:,:,:,:,:);
                case 2, ss = ss([end,1:end-1],:,:,:,:,:);
                case 3, ss = ss([end,1:end-1],[end,1:end-1],:,:,:,:);
                otherwise, warning('rot90fac(=%g)~=0,1,2,3',rot90fac);
            end
        otherwise, warning('trnsps(=%g)~=0,3',trnsps);
    end
    ss = ss(:,:,[end,1:end-1],:,:,:);
end


%% reorder to original shape
spec = permute(ss,[6,1,2,3,4,5]);


end      % csi_orientation.m
