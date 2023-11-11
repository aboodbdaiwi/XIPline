function psnr = getCroppedPSNR(y, y_est, half_kernel)
    psnr = getPSNR(y(half_kernel(1)+1:end-half_kernel(1), half_kernel(2)+1:end-half_kernel(2), :), ...
        y_est(half_kernel(1)+1:end-half_kernel(1), half_kernel(2)+1:end-half_kernel(2), :));
end