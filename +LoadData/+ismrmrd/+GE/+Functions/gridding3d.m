function [img, gdat] = gridding3d(kvec,kdata,dcf,imgN,a,kwidth,padN,vdim,crop)

% ========================================================================
% [img, gdat] = grecon3d(kvec, kdata, dcf, imgN, a, kwidth, padN, vdim)
% Uses mex function gr3dKB_mex()
%     =INPUT=
%     kvec   -- Npts x 3, k-trajectory k[x y z], scaled -0.5 to 0.5
%     kdata  -- Npts x 1, k-space data, 
%     dcf    -- Npts x 1, k-space weighting
%     imgN   -- 3x1, image size in each dim, [n1 n2 n3] for [x y z]
%                    (grid will be a*n1 x a*n2 x a*n3) 
%     a      -- 1x1, oversampling ratio 
%     kwidth -- 1x1, kernel width 
%     =optional=
%     padN   -- 3x1, zero pad output image matrix to pad1 x pad2 x pad3
%     vdim   -- 3x1, relative voxel dimensions in [x y z], e.g.,
%                        [1 1 1] : isotropic resolution
%                        [1 1 4] : z has 4x coarser resolution than xy
%                    it's assumed that the dim w/finest res has k-traj
%                        scaled to max abs value of 0.5
%     crop   -- 1x1, whether to crop to desired imaging FOV
%
%     =OUTPUT=
%     img    -- (n1 x n2 x n3) or (pad1 x pad2 x pad3), reconstructed image  
%     gdat   -- a*(n1 x n2 x n3) or a*(pad1 x pad2 x pad3), gridded k-space data 
%
% NOTE: So far, have only tested for cases where n1==n2
% ========================================================================
% Modified:
%           2010/08/27 new input param "crop"
%           2009/11/20 better input checking and handling for padN
%           2009/11/10 better support for anisotropic spatial resolution
%           2009/11/08 now takes padN instead of pad
%                      (assimilates function of grecon3dPad.m)
% Holden H. Wu; holdenhwu@stanford.edu
%           2/2017 RFS bug fixes
% ========================================================================

if( nargin<5 )
    a = 1.25; % 1.25X grid
end
if( nargin<6 )
    kwidth = 4;
end
if( nargin<7 )
    padN = imgN;
end
if( nargin<8 )
    vdim = [1 1 1];
end
vdim = vdim / min(vdim);
if( nargin<9 )
    crop = 1;
end    

if( numel(imgN)==1 )
    imgN = imgN*ones(3,1); % Assume same size in each dim
end
if( numel(padN)==1 )
    padN = padN*ones(3,1); % Assume same size in each dim
end 
% NOTE: Should do some more checking.

kvec = squeeze(kvec);
dcf = dcf.';
kdata = kdata.';

a = ceil(a*imgN/2)*2/imgN;  % RFS: round up to multiple of 2; recalc a

% Pre-compute KB kernel, cf. Beatty et al.
b = pi*sqrt( (kwidth/a)^2*(a-0.5)^2 - 0.8  ); % b=11.441 for a=2&kwidth=5
G = 1;
% sample the kernel and store as a lookup table
max_eps1 = 1e-4;
S = sqrt(0.37/max_eps1)/a; % linear interp
S = 10*ceil( S/10 ); % S samples per grid unit
kx_samp = 0:1/(G*S):kwidth/(2*G); % only need one side
C_samp = (G/kwidth) * besseli( 0, b*sqrt(1-(2*G*kx_samp/kwidth).^2) );

% 3D Gridding
% convert to single column
% kdata = kdata(:);
% kdata = kdata(:) + 1i*eps; % make sure it's complex
kdata = kdata(:);                                    % rfs
if isreal(kdata), kdata = complex(kdata); end        % rfs

dcf   = dcf(:);
% initialize k-space grid
% grsize = ceil(a*imgN); %grcenter = (a*imgN/2+1);
grsize = floor(a*imgN+0.5); % RFS
gdat   = zeros(grsize(1), grsize(2), grsize(3));
% call mex function
gdat = gr3dKB_mex(kvec(:,1), kvec(:,2), kvec(:,3), kdata, dcf, grsize(:), kwidth/2, C_samp, vdim);

% Zero-pad to larger image matrix
if( padN(1)==imgN(1) && padN(2)==imgN(2) && padN(3)==imgN(3) )
    changeN = 0;
else
    changeN = 1;
end    
if( changeN && padN(1)>=imgN(1) && padN(2)>=imgN(2) && padN(3)>=imgN(3) )
    grsizePad = ceil(a*padN);
    gdat2 = zeros( grsizePad(1), grsizePad(2), grsizePad(3) );
    x0 = floor( (grsizePad-grsize)/2 );
    gdat2(x0(1)+[1:grsize(1)], x0(2)+[1:grsize(2)], x0(3)+[1:grsize(3)]) = gdat;
    imgN = padN;
    grsize = grsizePad;
    gdat = gdat2;
    clear gdat2;
end

% check data and exclude nan and inf
if any(isnan(gdat(:))) 
    warning('isnan(gdat); nulling NaNs'); 
    gdat(isnan(gdat)) = 0; 
end
if any(isinf(gdat(:)))
    warning('isinf(gdat); nulling Infs'); 
    gdat(isinf(gdat)) = 0; 
end


% 3D-IFFT
img = fftshift( ifftn( ifftshift(gdat) ) );

% Deapodization
% using analytical expression of FT{kernel}, cf. Beatty et al.
argvecX = (pi*kwidth*[-grsize(1)/2:grsize(1)/2-1]/grsize(1)).^2 - b^2;
argvecY = (pi*kwidth*[-grsize(2)/2:grsize(2)/2-1]/grsize(2)).^2 - b^2;
argvecZ = (pi*kwidth*[-grsize(3)/2:grsize(3)/2-1]/grsize(3)).^2 - b^2;

% create deap matrix
[argX, argY, argZ] = ndgrid(argvecX, argvecY, argvecZ);
c_3D = ones( size(img) );    
c_3D = c_3D .* (sin(sqrt( argX )) ./ sqrt( argX ));
c_3D = c_3D .* (sin(sqrt( argY )) ./ sqrt( argY ));
c_3D = c_3D .* (sin(sqrt( argZ )) ./ sqrt( argZ ));
c_3D = c_3D + 10000; % baseline offset added to deap function
% now deapodize
img = img ./ c_3D;

% Extract FOV
if( crop )
    x0  = floor( (a-1)*imgN/2 );
    img = img( x0(1)+[1:imgN(1)], x0(2)+[1:imgN(2)], x0(3)+[1:imgN(3)] );
end 

return;
