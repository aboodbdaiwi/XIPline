function bb2 = grad_warp(bb,h,gradcoil,threeD,gpath,zfov,verb,do_gwc)
%GRAD_WARP  Distort images to correct for gradient non-linearities
%     bb2 = grad_warp(bb,h,gradcoil,threeD,gpath,zfov,verb,do_gwc)
%      bb   n-dimensional real data: [x,y,z,rep]
%       h   p-file header structure
%gradcoil   Gradient coil                                              ('')
%           BRM, HRMBm HRMW, IRMW, ORM, VRMW, XRM, XRMW, RSNA, MAGNUS
%  threeD   Do 3D gradwarp                                          (false)
%   gpath   Path to gw_coils files
%           default: linux='/w/config', win='C:\BoxSync\data\gw_coils'
%    zfov   z-FOV for 3D-fidall to fill h.data_acq_tab [mm] (h.rdb_hdr.fov)
%    verb   Verbose mode                                             (true)
%  do_gwc   Load coefficients from gw_coils.dat                      (true) 
%
%     bb2   Warped images (real+single)
%
% See also GERecon('Gradwarp'), get_gradcoil.
%  2/2025 Rolf Schulte
if (nargin<1), help(mfilename); return; end


%% input parameters
if ~exist('gradcoil','var'), gradcoil = ''; end
if ~exist('threeD','var'),   threeD = []; end
if isempty(threeD),          threeD = false; end
if ~exist('gpath','var'),    gpath = ''; end
if ~exist('zfov','var'),     zfov = []; end
if ~exist('verb','var'),     verb = []; end
if isempty(verb),            verb = true; end
if ~exist('do_gwc','var'),   do_gwc = []; end
if isempty(do_gwc),          do_gwc = true; end


%% extract gradcoil from header.mrconfig
if isempty(gradcoil), gradcoil = get_gradcoil(h); end
if verb, fprintf('gradcoil = %s\n',gradcoil); end


%% sizes; reshape to 4D for >4D
nn = size(bb);
if length(nn)>4
    rshp = true;
else
    rshp = false;
end
nn = [nn ones(1,4-length(nn))];
nn4 = prod(nn(4:end));
if rshp
    bb = reshape(bb,[nn(1) nn(2) nn(3) nn4]);
end


%% abs+single; negative real values?
if ~isreal(bb)
    fprintf('bb complex: warping real+imag independently\n');
    data_positive = false;
    bb = cat(4,real(bb),imag(bb));     % convert to real matrix
    nn4 = 2*nn4;
    data_complex = true;
else
    if all(bb(:)>=-eps)
        data_positive = true;
    else
        data_positive = false;
    end
    data_complex = false;
end
if ~isa(bb,'single')
    warning('bb not single: converting to single');
    bb = single(bb);
end
if threeD && (abs(h.image.scanspacing)>0.01)
    warning('threeD && (abs(h.image.scanspacing(=%g))>0.01)',...
        h.image.scanspacing);
end


%% expand h.data_acq_tab to populate 3D slices acquired with fidall
if threeD && size(h.data_acq_tab.gw_point1,1)==1
    if verb, fprintf('Expanding h.data_acq_tab to %d slices\n',nn(3)); end
    sgn = 1;       % sign validated with fidspiral; not sure about fidall-3D
    h = sub_expand_dat(h,nn(3),zfov,sgn,verb);
end
if size(h.data_acq_tab.gw_point1,1)<nn(3)
    error('size(h.data_acq_tab.gw_point1,1)(=%d)<nn(3)(=%d)',...
        size(h.data_acq_tab.gw_point1,1),nn(3));
end

%% read in gw coefficients from file
if do_gwc
    gw = sub_read_gw_coils(gradcoil,gpath,verb);
end


%% start actual gradient warping
bb2 = zeros([nn(1) nn(2) nn(3) nn4],'single');
if threeD
    if verb, fprintf('3D gradwarp\n'); end
    corners = [sub_get_corners(h,1) sub_get_corners(h,nn(3))];
    for l4=1:nn4
        if do_gwc
            bb2(:,:,:,l4) = GERecon('Gradwarp',bb(:,:,:,l4),corners,...
                'SphericalHarmonicCoefficients',gw.coeff);
        else
            bb2(:,:,:,l4) = GERecon('Gradwarp',bb(:,:,:,l4),corners,gradcoil);
        end
    end
else
    if verb, fprintf('2D gradwarp\n'); end
    for l3=1:nn(3)
        corners = sub_get_corners(h,l3);
        for l4=1:nn4
            if do_gwc
                bb2(:,:,l3,l4) = GERecon('Gradwarp',bb(:,:,l3,l4),corners,...
                    'SphericalHarmonicCoefficients',gw.coeff);
            else
                bb2(:,:,l3,l4) = GERecon('Gradwarp',bb(:,:,l3,l4),corners,...
                    gradcoil);
            end
        end
    end
end
if ~isreal(bb2), warning('~isreal(bb2)'); end
if data_positive
    sneg = abs(sum(bb2(bb2<0)));
    spos = sum(bb2(bb2>0));
    if verb && any(bb2(:)<0)
        fprintf('Zeroing negative bb2: sum(neg)/sum(pos)=%g\n',sneg/spos);
    end
    if sneg/spos>0.01
        warning('sum(negative(bb2))/sum(positive(bb2))(=%g)>0.01',sneg/spos);
    end
    bb2(bb2<0) = 0;
end


%% split again into complex numbers
if data_complex
    bb2 = complex(bb2(:,:,:,1:(nn4/2)),bb2(:,:,:,(1:nn4/2)+nn4/2));
end


%% reshape to original data size
if rshp
    bb2 = reshape(bb2,nn);
end


%% final checks
if any(isnan(bb2(:))), warning('isnan(bb2)'); end
if any(isinf(bb2(:))), warning('isinf(bb2)'); end

end      % grad_warp.m


%% sub-functions
function corners = sub_get_corners(h,slice)
%SUB_GET_CORNERS Convert p-file header to corners (gradwarp input) structure

corners.UpperLeft  = single(squeeze(h.data_acq_tab.gw_point1(slice,1,:))).';
corners.UpperRight = single(squeeze(h.data_acq_tab.gw_point3(slice,1,:))).';
corners.LowerLeft  = single(squeeze(h.data_acq_tab.gw_point2(slice,1,:))).';
corners.Type = 'SliceCorners';
corners.Others = single([...
    h.data_acq_tab.freq_loc_shift(slice),...
    h.data_acq_tab.phase_loc_shift(slice),...
    h.data_acq_tab.slice_loc_shift(slice),...
    h.data_acq_tab.fov_freq_scale(slice),...
    h.data_acq_tab.fov_phase_scale(slice),...
    h.data_acq_tab.slthick_scale(slice)]);

end      % sub_get_corners


%% sub-function to read in grad warp coefficients
function gw  = sub_read_gw_coils(gradcoil,gpath,verb)
%READ_GW_COILS Reads gw_coils.dat file /w/config/gw_coils/
%    gw = read_gw_coils(gradcoil,verb)
%gradcoil   Gradien coil
%    verb   Verbose mode
%      gw   Output structure with grad warp info
%           coeff: Nx3 matrix with all of the coefficents
%           type:  gradwarptype
%           delta: 
%           line1: grad coil info
%
% Author: Jason A. Polzin GE Medical Systems
if (nargin<1), help(mfilename); return; end

if exist(gradcoil,'file')
    fname = gradcoil;    % fname_gwc=input
else                         % assemble fname_gwc
    if isempty(gpath)
        if isdeployed
            gpath = '/w/config';
        else
            gpath = 'C:\BoxSync\data\gw_coils';
        end
    end
    if ~exist(gpath,'dir'), error('dir gpath not found: ''%s''',gpath); end
    fname = [gpath filesep 'gw_coils.dat'];      % default on scanner
    if ~isempty(gradcoil) && ~isdeployed, fname = [fname '.' upper(gradcoil)]; end
    if ~exist(fname,'file')
        error('fname (=''%s'') not found',fname);
    end
end

if ~exist(fname,'file')
    error('fname (=''%s'') not found',fname);
end


%% read file
if verb, fprintf('Reading gradwarp file ''%s''\n',fname); end
fp = fopen(fname, 'rt');
if (fp == -1), error('Cannot open ''%s''',fname); end

ax = 'XYZ';             % axis labels
gw.line1 = fgetl(fp);   % read first line
if verb, fprintf('%s\n',gw.line1); end
% Use fgetl to grap entire line and save as string.
gw.type = sscanf(fgetl(fp), 'GRADWARPTYPE%d');
gw.coeff = zeros(10,3);
for j = 1:3
    for i = 1:10
        % This reads in a line of the form SCALE<XYZ>C where C is a
        % coefficient for GradWarp.
        gw.coeff(i,j) = sscanf(fgetl(fp), sprintf('SCALE%c%d%%g\n',ax(j), i));
    end
end
gw.delta = sscanf(fgetl(fp), 'DELTA%d');

fclose(fp);


%% final checks
if (gw.type ~= 1), warning('GRADWARPTYPE(=%g)~=1',gw.type); end
if abs(gw.delta)>1d-5, warning('DELTA(=%g)~=0',gw.delta); end


end      % sub_read_gw_coils.m


%% expand h.data_acq_tab.gw_point for 3D with fidall
function  h = sub_expand_dat(h,nz,zfov,sgn,verb)

warning('fidall-3D: not sure about direction of perpendicular: sign=%d',sgn);
if isempty(zfov)
    fprintf('Attention: zfov=[]; setting=h.rdb_hdr.fov(=%g)\n',h.rdb_hdr.fov);
    zfov = h.rdb_hdr.fov;
end

% extract gw corner points (z-direction=centre)
for ll=1:3
    fld = ['gw_point' num2str(ll)];
    gwp{ll} = h.data_acq_tab.(fld);
    if size(gwp{ll},2)==3
        warning('size(gwp{%d},2)==3: reshaping',ll);
        gwp{ll} = reshape(gwp{ll},[1 1 3]);
    end
end

% calculate perpendicular line (z-direction)
perpen = cross(gwp{3}-gwp{1},gwp{2}-gwp{1});     % through-slice direction
perpen = perpen/sqrt(sum(perpen.^2));            % normalise
if verb
    fprintf('perpendicular = %s\n',num2str(squeeze(perpen).'));
end
zloc_obl = bsxfun(@times,(-nz/2:nz/2-1).'/nz*zfov,perpen);

for ll=1:3
    fld = ['gw_point' num2str(ll)];
    h.data_acq_tab.(fld) = bsxfun(@plus,gwp{ll},sgn*zloc_obl);
    if verb>1
        fprintf('%s\n',fld);
        disp(num2str(squeeze(h.data_acq_tab.(fld))));
    end
end

end      % sub_expand_dat
