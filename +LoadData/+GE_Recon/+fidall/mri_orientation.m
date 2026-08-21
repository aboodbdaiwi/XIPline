function bb = mri_orientation(bb,rot90fac,trnsps,do_voxshift)
%MRI_ORIENTATION  Orient images (xy plane) from xyz to AP-RL-SI coordinates
%         bb = mri_orientation(spec,rot90fac,trnsps,do_voxshift)
%         bb   Images        [#x,#y,#z,#timesteps,#coils]
%   rot90fac   h.data_acq_tab.rotate
%     trnsps   h.data_acq_tab.transpose
%do_voxshift   Shift back by one voxel depending on orientation     (false)
%
%  2/2025 Rolf Schulte
if nargin<1, help(mfilename); return; end


%% misc parameters + checks
if ~exist('do_voxshift','var'), do_voxshift = []; end
if isempty(do_voxshift),        do_voxshift = false; end
if length(size(bb))>5,  warning('bb(=%dD)>5D',length(size(bb))); end
if any(diff(rot90fac)), warning('rot90fac'); disp(rot90fac); end
if any(diff(trnsps)),   warning('trnsps');   disp(tsnsps); end
rot90fac = rot90fac(1);
trnsps = trnsps(1);
fprintf('mri_orientation: rot90fac=%g trnsps=%g\n',rot90fac,trnsps);
if ~any(trnsps==[0,3]), warning('trnsps(=%g)~=[0,3]',trnsps); end
if ((rot90fac<0) || (rot90fac>3))
    warning('rot90fac(=%g)<0 or >3',rot90fac);
end


%% orientation
% moved fliplr(bb) to pre-gradwarp
if get_mver>8
    if trnsps==3
        bb = rot90(permute(bb,[2 1 3 4 5 6]),rot90fac);
    else
        bb = rot90(bb,rot90fac);
    end
else
    mtx = size(bb);
    mtx = [mtx , ones(1,5-length(mtx))];
    if trnsps(1)==3
        if mod(rot90fac(1),2)<1
            mtx = mtx(1,[2,1,3,4,5]);
        end
    else
        if mod(rot90fac(1),2)==1
            mtx = mtx(1,[2,1,3,4,5]);
        end
    end
    bb1 = complex(zeros(mtx,class(bb)));
    for l3=1:size(bb,3)
        for l4=1:size(bb,4)
            for l5=1:size(bb,5)
                for l6=1:size(bb,6)
                    % image orientation
                    if trnsps(1)==3
                        bb1(:,:,l3,l4,l5,l6) = rot90(bb(:,:,l3,l4,l5,l6).',rot90fac(1));
                    else
                        bb1(:,:,l3,l4,l5,l6) = rot90(bb(:,:,l3,l4,l5,l6),rot90fac(1));
                    end
                end
            end
        end
    end
    bb = bb1;
    clear bb1
end


%% correct one voxel shift
% test via:
% ss=zeros(1,4,4);ss(1,3,3)=1;
% ss2=csi_orientation(ss,0,0,true);squeeze(ss2)
if do_voxshift
    warning('big mess: shift half a voxel before reco instead via shft!!!');
    input('');
    switch trnsps
        case 0
            switch rot90fac
                case 0, bb = bb([2:end,1],:,:,:,:,:);           % ok
                case 1, bb = bb([end,1:end-1],[end,1:end-1],:,:,:,:);
                case 2, bb = bb(:,[end,1:end-1],:,:,:,:);
                case 3, bb = bb(:,[end,1:end-1],:,:,:,:);       % ok
                otherwise, warning('rot90fac(=%g)~=0,1,2,3',rot90fac);
            end
        case 3
            switch rot90fac
                case 0, bb = bb(:,[2:end,1],:,:,:,:);           % ok
                case 1, bb = bb([end,1:end-1],:,:,:,:,:);       % ok
                case 2, bb = bb([end,1:end-1],:,:,:,:,:);
                case 3, bb = bb([end,1:end-1],[end,1:end-1],:,:,:,:);
                otherwise, warning('rot90fac(=%g)~=0,1,2,3',rot90fac);
            end
        otherwise, warning('trnsps(=%g)~=0,3',trnsps);
    end
    % bb = bb(:,:,[end,1:end-1],:,:,:);
end


end      % mri_orientation.m
