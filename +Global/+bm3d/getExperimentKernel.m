function kernel = getExperimentKernel(noiseType, noiseVar, sz)

% Case gw / g0
kernel = ones(1);

noiseTypes = {'gw', 'g0', 'g1', 'g2', 'g3', 'g4', 'g1w', 'g2w', 'g3w', 'g4w'};

found = false;
for nt = noiseTypes
    if strcmp(noiseType, nt)
        found = true;
    end
end
if ~found
    disp('Error: Unknown noise type!')
    return;
end

if (~strcmp(noiseType, 'g4') && ~strcmp(noiseType, 'g4w')) || ~exist('sz', 'var')
    % Crop this size of kernel when generating,
    % unless pink noise, in which
    % Case we want to use the full image size
    sz = [101, 101];
   
end

% Sizes for meshgrids
sz2 = -(1 - mod(sz, 2)) * 1 + floor(sz/2);
sz1 = floor(sz/2);
[uu, vv] = meshgrid(-sz1(1):sz2(1), -sz1(2):sz2(2)); 
alpha = 0.8;

switch (noiseType(1:2))
    case 'g1'
        % Horizontal line
        kernel = 16 - abs((1:31)-16);
        
    case 'g2'
        % Circular repeating pattern
        scale = 1;
        dist = (uu).^2 + (vv).^2;
        kernel = cos(sqrt(dist) / scale) .* fspecial('gaussian', [sz(1), sz(2)], 10);
        
    case 'g3' 
        % Diagonal line pattern kernel
        scale = 1;
        kernel = cos((uu + vv) / scale) .* fspecial('gaussian', [sz(1), sz(2)], 10);    
        
    case 'g4'
        % Pink noise
        dist = (uu).^2 + (vv).^2;
        n = sz(1)*sz(2);
        spec3 = sqrt((sqrt(n)*1e-2)./(sqrt(dist) +  sqrt(n)*1e-2));
        kernel = fftshift(ifft2(ifftshift(spec3)));
end    

% -- Noise with additional white component --
if numel(noiseType) == 3 && noiseType(3) == 'w'
    kernel = kernel / norm(kernel(:));
    kalpha = sqrt((1 - alpha) + (alpha) * abs(fft2(kernel, sz(1), sz(2))).^2);
    kernel = fftshift(ifft2(kalpha));
end

% Correct variance
kernel = kernel / norm(kernel(:)) * sqrt(noiseVar);
