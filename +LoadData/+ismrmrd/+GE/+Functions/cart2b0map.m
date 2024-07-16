function [df0,bb1,bb2,mask,df0fit] = cart2b0map(d1,d2,dte,res,fac_thresh,fftshiftdir,offset)
%CART2B0MAP  Reconstruct B0 phase map of Cartesian sampled data
%[df0,b1,b2,mask] = cart2b0map(d1,d2,dte,res,fac_thresh,fftshiftdir,offset)
%
%         d1  data for TE1
%         d2  data for TE2
%        dte  Delta TE = TE2 - TE1                   [s]
%        res  anticipated resolution of B0 map
% thresh_fac  factor on threshold for mask
%fftshiftdir  1=1;2=2;3=all;rest=none
%     offset  phase offset (to move phase jumps out) [-pi pi]
%             optional; if omitted, calc offset with minimal jumps
%
%        df0  B0 map                                 [Hz]
%         b1  Reconstructed image (d1)
%         b2  Reconstructed image (d2)
%       mask  mask (fill_holes(df0~=0))
%
%  10/2015 Rolf Schulte
if (nargin<3), help(mfilename); return; end;

if ~exist('res','var'),         res = []; end
if ~exist('offset','var'),      offset = []; end
if ~exist('fac_thresh','var'),  fac_thresh = []; end
if isempty(fac_thresh),         fac_thresh = 1; end
if ~exist('fftshiftdir','var'), fftshiftdir = []; end
if isempty(fftshiftdir),        fftshiftdir = 1; end


% rotation factor: 
%   0=no rotation
%   1=90 anti-clock
rot90fac = 0;

if dte>0.1, warning('cart2b0map:dte','dte>0.1: unit in s?'); end

if any(size(d1)~=size(d2)), 
    warning('cart2b0map:sizes','size(d1)~=size(d2): setting equal'); 
    tmpres = max([size(d1) ; size(d2)]);
else
    tmpres = [];
end

%% reconstruction of images
b1 = recon2d(d1,tmpres,fftshiftdir);
b2 = recon2d(d2,tmpres,fftshiftdir);
b1 = conj(b1);  % required to have right convention for phase (mfi_recon)
b2 = conj(b2);

si = size(b1);
n_coils = size(d1,6);
if n_coils~=size(b1,3), 
    warning('n_coils(=%g) = size(b1,3)(=%g)',n_coils,size(b1,3)); 
end
switch length(res)
    case 0, res = [si(1) si(2)];
    case 1, res = [res res];
    case 2,
    otherwise, error('length(res) ~=1 or 2');
end

%% calculate mask
bb1 = sqrt(sum(b1.*conj(b1),3));
bb2 = sqrt(sum(b2.*conj(b2),3));
tmp = [bb1(1,:) bb1(end,:) bb1(:,1).' bb1(:,end).' ...
    bb2(1,:) bb2(end,:) bb2(:,1).' bb2(:,end).'];
thresh2 = (1.5*mean(tmp)+5*std(tmp))*fac_thresh;
fprintf('Mask threshold = %g\n',thresh2);

mask = ((bb1+bb2)/2)>thresh2;
% clf; subplot(2,2,1); imagesc(mask);

if n_coils>1, 
    warning('Coil combination of phase not SNR optimal'); 
    if isempty(offset), 
        warning('Noise will lead to imprecise B0 maps'); 
        fprintf('Switching off pump calculation (offset=0)\n');
        offset = 0;
    end
end



%% determine offset with least phase jumps
if isempty(offset),
    min_jumps = si(1)*si(2);
    larr = (-180:10:180); larr = larr(randperm(length(larr)));
    for l=larr                                                      % [deg]
        dphi = mean(angle(conj(b1).*b2*exp(1i*pi/180*l)),3);        % [rad]
        jumps1 = mask.*diff([zeros(1,si(2)) ; dphi],1,1);
        jumps2 = mask.*diff([zeros(si(1),1) , dphi],1,2);
        jumps = sum(sum(((abs(jumps1)>1.9*pi)+(abs(jumps2)>1.9*pi)),1),2)/2;
        
        if min_jumps>jumps
            offset = pi/180*l;                                      % [rad]
            min_jumps = jumps;
        end
        if jumps==0, break; end
    end
    fprintf('min_jumps = %g \tphase offset = %g [rad]\n',min_jumps,offset);
end

%% generation of phase maps
dphi = mean(angle(conj(b1).*b2*exp(1i*offset))-offset,3);           % [rad]
df0 = (dphi/(2*pi*dte)).*mask; 
df0off = mean([ mean(df0(mask~=0)) ...
    median(df0(mask~=0)) ...
    df0(floor(si(1)/2),floor(si(2)/2))]);
fprintf('df0off = %g (subtracting from df0)\n',df0off);
half_phaseHz = 0.5/dte;
if abs(df0off)/half_phaseHz>0.8, 
    warning('cart2b0map:ssfp','df0off=%g; ssfp_flag=1? subtracting %g',...
        df0off,half_phaseHz); 
    df0 = (df0-half_phaseHz).*mask;
end
df0 = (df0-df0off).*mask;

% resample (down)
if any(res~=[si(1) si(2)]),
    switch 2
        case 1,
            x = (-si(2)/2:si(2)/2-1).'/si(2);
            y = (-si(1)/2:si(1)/2-1)/si(1);
            xi = (-res(2)/2:res(2)/2-1).'/res(2);
            yi = (-res(1)/2:res(1)/2-1)/res(1);
            df0 = interp2(x,y,df0,xi,yi);
            mask = abs(df0)>1d-10;
        case 2,
            % mask interpolation: interp2
            x = (-si(2)/2:si(2)/2-1).'/si(2);
            y = (-si(1)/2:si(1)/2-1)/si(1);
            xi = (-res(2)/2:res(2)/2-1).'/res(2);
            yi = (-res(1)/2:res(1)/2-1)/res(1);
            mask = interp2(x,y,double(mask),xi,yi);
            mask(isnan(mask)) = 0;
            mask = logical(abs(mask)>1d-10);

            % df0: (filtered) Fourier interpolation 
            i11 = ceil(res(1)/2)-1;
            i12 = ceil(si(1)-(res(1)/2));
            i21 = ceil(res(2)/2)-1;
            i22 = ceil(si(2)-(res(2)/2));
            tmp = ifft2(fftshift(df0));
            tmp = tmp([1:i11 i12:end],[1:i21 i22:end]);
            if false
                filt = lb_fun(res(1),1000,[0 10],res(1)/2/1000).'* ...
                    lb_fun(res(2),1000,[0 10],res(2)/2/1000);
                tmp = tmp.*fftshift(filt);
                clf; imagesc(filt); drawnow; input('');
            end
            df0 = fftshift(fft2(tmp));
            fprintf('mean(real(df0(:))/mean(imag(df0(:))=%g\n',...
                mean(real(df0(:)))/mean(imag(df0(:))));
            df0 = real(df0).*mask;
    end
end

if any(isnan(df0(:))),
    warning('cart2b0map:nan','df0 contains NaN: setting to zero'); 
    df0(isnan(df0)) = 0;
end

if rot90fac~=0,
    df0  = rot90(df0,rot90fac);
    bb1 = rot90(bb1,rot90fac);
    bb2 = rot90(bb2,rot90fac);
    dphi = rot90(dphi,rot90fac);
end

mask = logical(fill_holes(mask));
df0fit = fiex(df0);

% plotting
clf;
set(gcf,'name','cart2b0map');
b1m = sort(bb1(:)); b1m = 1.2*b1m(floor(0.97*si(1)*si(2)));
b2m = sort(bb2(:)); b2m = 1.2*b2m(floor(0.97*si(1)*si(2)));
subplot(3,2,1); imagesc(bb1,[0 b1m]);
axis square off; title('b1'); colorbar;
subplot(3,2,2); imagesc(bb2,[0 b2m]);
axis square off; title('b2'); colorbar;
subplot(3,2,3); imagesc(dphi,[-pi pi]); 
axis square off; title('dphi'); colorbar;
subplot(3,2,4); imagesc(df0);
axis square off; title('B0map [Hz]'); caxis(0.7*[-1 1]*max(abs(df0(:)))); colorbar;
subplot(3,2,5); imagesc(mask); 
axis square off; title('mask'); colorbar;
subplot(3,2,6); imagesc(df0);
axis square off; title('fitted B0map [Hz]'); caxis(0.7*[-1 1]*max(abs(df0(:)))); colorbar;
drawnow;


