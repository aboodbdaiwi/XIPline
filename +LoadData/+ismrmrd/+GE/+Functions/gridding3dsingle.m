function [img, gdat] = gridding3dsingle(kvec,kdata,dcf,imgN,a,kwidth,padN,vdim,crop)
%GRIDDING3DSINGLE 3D Gridding reconstruction
% [img, gdat] = grecon3dsingle(kvec, kdata, dcf, imgN, a, kwidth, padN, vdim)
% Uses mex function gr3dKB_float_mex()
%     =INPUT=
%     kvec   -- Npts x 3, k-trajectory k[x y z], scaled -0.5 to 0.5
%     kdata  -- Nrep x Npts, k-space data
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
% Literature: itomi24_799_2005_beatty.pdf
% Holden H. Wu; holdenhwu@stanford.edu
%
% 1/2020 RFS bug fixes + converted to single
% ========================================================================
if (nargin<1), help(mfilename); return; end

% a = 1.25  for kwidth=4
% a = 1.125 for kwidth=3
if nargin<5, a = 1.125; end   % 1.25X grid 
if nargin<6, kwidth = 3; end
if kwidth<2.001, warning('kwidth(=%g)<2.001',kwidth); end
if nargin<7, padN = imgN; end
if nargin<8, vdim = [1 1 1]; end
vdim = vdim/min(vdim);
if nargin<9, crop = true; end
if numel(imgN)==1, imgN = imgN*ones(3,1); end % Assume same size in each dim
if numel(padN)==1, padN = padN*ones(3,1); end % Assume same size in each dim
% NOTE: Should do some more checking.
fufa_deapo = 0;   % RFS: fudge factor for deapodisation; originally 10000
do_debug = false;    % nan+inf -> saving variable for debugging 

nrep = size(kdata,1);             % #repetitions (eg coils)
if ~isa(kvec, 'single'), kvec =  single(kvec); end
if ~isa(dcf,  'single'), dcf =   single(dcf); end
if ~isa(kdata,'single'), kdata = single(kdata); end
if isreal(kdata),        kdata = complex(kdata); end

kvec = squeeze(kvec);
dcf = dcf.';

a = ceil(a*imgN/2)*2/imgN;        % RFS: round up to multiple of 2; recalc a


%% Pre-compute KB kernel (Beatty Eq 5)
b = pi*sqrt( (kwidth/a)^2*(a-0.5)^2 - 0.8  ); % b=11.441 for a=2&kwidth=5
G = 1;
% sample the kernel and store as a lookup table
max_eps1 = 1e-4;
S = sqrt(0.37/max_eps1)/a;        % linear interp
S = 10*ceil( S/10 );              % S samples per grid unit
kx_samp = 0:1/(G*S):kwidth/(2*G); % only need one side
C_samp = (G/kwidth) * besseli( 0, b*sqrt(1-(2*G*kx_samp/kwidth).^2) );


%% 3D Gridding
dcf = dcf(:);
% initialize k-space grid
% grsize = ceil(a*imgN); %grcenter = (a*imgN/2+1);
grsize = floor(a*imgN+0.5); % RFS
if false
    % c-loop faster than parfor
    gdat = complex(zeros(grsize(1), grsize(2), grsize(3),nrep,'single'));
    parfor l=1:nrep
        gdat(:,:,:,l) = gr3dKB_float_mex(kvec(:,1), kvec(:,2), kvec(:,3), ...
            kdata(l,:), dcf, grsize(:), kwidth/2, single(C_samp), single(vdim));
    end
else
    % call mex function
    gdat = gr3dKB_float_mex(kvec(:,1), kvec(:,2), kvec(:,3), ...
        kdata, dcf, grsize(:), kwidth/2, single(C_samp), single(vdim));
end

%% Zero-pad to larger image matrix
% if padN(1)==imgN(1) && padN(2)==imgN(2) && padN(3)==imgN(3)
if all(abs(padN-imgN)<eps)
    changeN = false;
else
    changeN = true;
end
if changeN && padN(1)>=imgN(1) && padN(2)>=imgN(2) && padN(3)>=imgN(3)
    grsizePad = ceil(a*padN);
    gdat2 = complex(zeros(grsizePad(1),grsizePad(2),grsizePad(3),nrep,'single'));
    x0 = floor( (grsizePad-grsize)/2 );
    gdat2(x0(1)+(1:grsize(1)), x0(2)+(1:grsize(2)), x0(3)+(1:grsize(3)),:) = gdat;
    imgN = padN;
    grsize = grsizePad;
    gdat = gdat2;
    clear gdat2;
end


%% check data and exclude nan and inf
if any(any(any(any(isnan(gdat))))) 
    warning('isnan(gdat); nulling NaNs');
    if do_debug
        save('gridding3dsingle_nan_debug.mat','-v7.3')
    end
    gdat(isnan(gdat)) = 0; 
end
if any(any(any(any(isinf(gdat)))))
    warning('isinf(gdat); nulling Infs');
    if do_debug
        save('gridding3dsingle_inf_debug.mat','-v7.3')
    end
    gdat(isinf(gdat)) = 0; 
end


%% 3D-IFFT
%img = fftshift( ifftn( ifftshift(gdat) ) );
img = fftshift(fftshift(fftshift(ifft(ifft(ifft(...
    ifftshift(ifftshift(ifftshift(gdat,3),2),1),[],3),[],2),[],1),3),2),1);
if nargout<2, clear gdat; end

%% Deapodization
% using analytical expression of FT{kernel}, cf. Beatty et al.
argvecX = (pi*kwidth*(-grsize(1)/2:grsize(1)/2-1)/grsize(1)).^2 - b^2;
argvecY = (pi*kwidth*(-grsize(2)/2:grsize(2)/2-1)/grsize(2)).^2 - b^2;
argvecZ = (pi*kwidth*(-grsize(3)/2:grsize(3)/2-1)/grsize(3)).^2 - b^2;

% create deap matrix
[argX, argY, argZ] = ndgrid(sqrt(argvecX), sqrt(argvecY), sqrt(argvecZ));
c_3D = (sin(argX)./argX).*(sin(argY)./argY).*(sin(argZ)./argZ)+fufa_deapo;
% now deapodize
% img = img ./ c_3D;
img = bsxfun(@rdivide,img,c_3D/(a^3));


%% Extract FOV
if crop
    x0  = floor( (a-1)*imgN/2 );
    img = img( x0(1)+(1:imgN(1)), x0(2)+(1:imgN(2)), x0(3)+(1:imgN(3)), :);
end

end  % main function gridding3dsingle.m

