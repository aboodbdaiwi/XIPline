function [b1fit,txAtten,txPhase] = b1optimiser(b1map,opt_mask,method,inp,verb)
%B1OPTIMISER  Optimise dual-drive settings with given B1+ map
%[b1fit,txAtten,txPhase] = b1optimiser(b1map,mask,method,inp,verb)
%    b1map  B1map (complex; arbitrary units)
%           size(b1map) -> [mtx1,mtx2,mtx3,ntx]; with ntx=2
%     mask  Mask for optimisation (optional; zeros in b1map equivalent)
%   method  Optimisation method
%           1=linear LSQ (requires inp = target b1)
%           2=lsqnonlin
%           3=exhaustive search
%           4=quadrature mode (x=[0,1]
%           5=input [txAtten txPhase]
%           cost function (2-5): squared deviation from mean of abs of 
%              combined map:
%              w = x(1)*exp(1i*x(2)); f = abs(BI + w*BQ); f = f - mean(f);
%      inp  input for methods 2 & 5
%     verb  verbose mode: 0=off; 1=verbose; 
%           2=verb+imagesc_row; 3=verb+imagesc_ind3d (defaults:2D=2,3D=3)
%
%    b1fit  Resulting b1map
%  txAtten  Amplitude difference between Q and I channel Tx  [dB/10]
%           txAtten = -200*log10(x(1));
%           Prescan UI: Amplitude Attenuation
%  txPhase  Phase difference between Q and I channel Tx      [deg]
%           txPhase = x(2)*180/pi; 
%           Prescan UI: Phase Delay
%
% See also BLOSI_B1MAP
% 6/2017 Rolf Schulte
if (nargin<1), help(mfilename); return; end;

if ~exist('opt_mask','var'), opt_mask = []; end
if ~exist('method','var'), method = []; end
if isempty(method),        method = 2; end
if ~exist('inp','var'),    inp = []; end
if ~exist('verb','var'),   verb = []; end
if ~isa(b1map,'double') && (method==2), 
    warning('lsqnonlin requires double; converting b1map to double'); 
    b1map = double(b1map);
end

%% boundaries given by scanner
maximum_attenuation = 45;         % maximum attenuation = 45 [dB/10]
minimum_attenuation = 0;          % minimum attenuation = 0 [dB/10]
maximum_phase_delay = 0;          % maximum phase = 0 [deg]
minimum_phase_delay = -45;        % minimum phase = 45 [deg]

% convert to complex polar numbers
r_min = 10^(-maximum_attenuation/10/20);
r_max = 10^(-minimum_attenuation/10/20);
phi_min = minimum_phase_delay*pi/180;
phi_max = maximum_phase_delay*pi/180;


%% misc parameters and checks
if isreal(b1map), warning('isreal(b1map)'); end
[mtx1,mtx2,mtx3,ntx] = size(b1map);
if ntx~=2, error('ntx(=%g)~=2',ntx); end
if mtx3>6    % only display 3D for z>6 
    if isempty(verb), verb = 3; end
else
    if isempty(verb), verb = 2; end    
end

if ~isempty(opt_mask),
    if ((size(opt_mask,1)~=mtx1) || (size(opt_mask,2)~=mtx2)) 
        warning('size mismatch mask & b1map');
    end
    if (size(opt_mask,3)~=mtx3)
        fprintf('Attention: size mismatch mask & b1map in dim=3\n');
        fprintf('\tassuming multi-slice, averaging over all\n');
    end
    b1map = bsxfun(@times,b1map,opt_mask);
end

% mask = (b1map~=0);
mask = abs(b1map)>eps;
mask = mask(:,:,:,1) & mask(:,:,:,2);

selvox = sum(mask(:))/numel(mask);
fprintf('Voxels selected by mask for b1optimiser: %g/%g = %g[%%]\n',...
    sum(mask(:)),numel(mask), selvox*100);
if (selvox<0.01), warning('<1%% voxels selected'); end
if (selvox==0), error('no voxels selected'); end


%% optimisation
BI = b1map(:,:,:,1);              % b1map I-channel
BI = BI(mask);                    % vectorise
BQ = b1map(:,:,:,2);              % b1map Q-channel
BQ = BQ(mask);                    % vectorise
b1fit = zeros(mtx1,mtx2,mtx3);    % resulting combined map


%% choose different optimisation/combination methods
switch method
    case 1,
        method_str = 'Linear LSQ';
        b1target = inp;
        if isempty(b1target), b1target = (mean(abs(BI)) + mean(abs(BQ)))/2; end
        A = BI./BQ;
        b = ones(mtx1,mtx2,mtx3)*b1target;
        b = b(mask);
        x1 = pinv(A)*b;
        x = [abs(x1),angle(x1)];
    case 2,
        method_str = 'Non-linear LSQ';
        x0 = [1,0];
        lb = [r_min,phi_min]; 
        ub = [r_max,phi_max];
        x = lsqnonlin(@(x)b1optim_costfun(x,BI,BQ),x0,lb,ub);
    case 3,
        method_str = 'Exhaustive search';
        nr = 100; np = 100;
        r_arr =   linspace(r_min,r_max,nr);
        phi_arr = linspace(phi_min,phi_max,np);
        f_arr = zeros(nr,np);
        for lr=1:nr,
            for lp=1:np,
                f = b1optim_costfun([r_arr(lr),phi_arr(lp)],BI,BQ);
                f_arr(lr,lp) = sum(f.^2);
            end
        end
        [tmp,ind] = min(f_arr(:)); 
        % fprintf('min(f) = %g\n',tmp);
        [x1,x2] = ind2sub(size(f_arr),ind); x = [r_arr(x1),phi_arr(x2)];
        if verb>1,
            % figure(20); 
            clf;
            % imagesc(phi_arr,r_arr,f_arr,[0 max(f_arr(:))]);
            % hold on; plot(phi_arr(x2),r_arr(x1),'xw'); hold off
            
            imagesc(phi_arr*180/pi,-200*log10(r_arr),f_arr,[0 max(f_arr(:))]);
            hold on; plot(180/pi*phi_arr(x2),-200*log10(r_arr(x1)),'xw'); hold off

            ylabel('r'); xlabel('phi [rad]'); title('Exhaustive search');
            colorbar
            figure;
        end
    case 4,
        method_str = 'Quadrature mode';
        x = [1,0];
    case 5,        
        if length(inp)~=2, warning('length(inp)(=%g)~=2',length(inp)); end
        method_str = sprintf('Combination with given input (%g %g)',inp);
        % x = inp;
        x = [10^(-inp(1)/200) inp(2)*pi/180];
        
    otherwise, error('method(=%g) not existing');
end


%% calculate and check resulting values
w = x(1)*exp(1i*x(2));            % complex weights (normal complex numb)
b1fit(mask) = BI+w*BQ;            % combined image
txAtten = -200*log10(x(1));       % Q-channel transmit attenuation [dB/10]
txPhase = x(2)*180/pi;            % Q-channel phase [deg]
f = b1optim_costfun(x,BI,BQ);
f = sum(f.^2);                    % cost function value
if txAtten < minimum_attenuation,
    warning('txAtten(=%g) < minimum_attenuation(=%g)',...
        txAtten,minimum_attenuation);
end
if txAtten > maximum_attenuation,
    warning('txAtten(=%g) > maximum_attenuation(=%g)',...
        txAtten,maximum_attenuation);
end
if txPhase < minimum_phase_delay,
    warning('txPhase(=%g) < minimum_phase_delay(=%g)',...
        txPhase,minimum_phase_delay);
end
if txPhase > maximum_phase_delay,
    warning('txPhase(=%g) > maximum_phase_delay(=%g)',...
        txPhase,maximum_phase_delay);
end


%% plotting
if verb>1,
    cmax = 1.5*max(abs([BI ; BQ]));    % maximum for plotting
    % figure(20+method);
    i1 = (sum(sum(b1fit,2),3)>0);
    i2 = (sum(sum(b1fit,1),3)>0);
    i3 = (sum(sum(b1fit,1),2)>0);
    
    clf
    if verb==2
        imagesc_row(abs(b1fit(i1,i2,i3)),[0 cmax],'',true,false); 
    else
        imagesc_ind3d(abs(b1fit(i1,i2,i3)),'',[0 cmax],'',false,false); 
    end
    colorbar
    title(sprintf('%s: f=%g; attn=%g; phase=%g',...
        method_str,f,floor(txAtten+0.5),floor(txPhase+0.5)));
    drawnow;
end


%% summarise results
if verb>0,
    fprintf('%s\n',method_str);
    fprintf('x(1) = r = %g; x(2) = phi = %g\n',x);
    fprintf('w = x(1)*exp(1i*x(2)) = %g + %gi\n',real(w),imag(w));
    fprintf('txAtten = %g [dB/10]\n',floor(txAtten+0.5));
    fprintf('txPhase = %g [deg]\n',floor(txPhase+0.5));
    fprintf('f = %g\n',f);
end


end   % main function b1optimiser.m



%% subfunctions
function f = b1optim_costfun(x,BI,BQ)

w = x(1)*exp(1i*x(2));            % complex weight
f = abs(BI + w*BQ);               % combine I and Q
f = f - mean(f);                  % minimise deviation from mean

% if ~isempty(mask),
%     b1_fit = zeros(size(mask));
%     b1_fit(mask) = f;
%     imagesc(abs(b1_fit),[0 3]); axis image off;
%     drawnow;
% end

end
