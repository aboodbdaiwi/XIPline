function [noise, PSD, kernel] = getExperimentNoise(noise_type, noise_var, realization, sz)

randn('seed',realization);

% Get pre-specified kernel
kernel = getExperimentKernel(noise_type, noise_var, sz);

% Create noisy image
half_kernel = ceil(size(kernel) ./ 2);
if(numel(sz) == 3 && numel(half_kernel) == 2)
    half_kernel(3) = 0;
end

% Crop edges
noise = convn(randn(sz + 2 * half_kernel), kernel(end:-1:1,end:-1:1, :), 'same');
noise = noise(1+half_kernel(1):end-half_kernel(1), 1+half_kernel(2):end-half_kernel(2), :);

PSD = abs(fft2(kernel, sz(1), sz(2))).^2 * sz(1) * sz(2);

end