function pdf_ref_trunc = ADC_RefGamma(age)
% ADC_RefGamma  Returns truncated Gamma reference PDF for a given age
% 
% INPUT:
%    age  = subject age in years (scalar)
%
% OUTPUT:
%    pdf_ref_trunc = truncated Gamma PDF evaluated on x = 0–0.14 (1×400 vector)
%
% NOTES:
% Uses linear-age equations for mean and SD:
%    Mean(age) = 0.0274639 + 0.000233964 * age
%    SD(age)   = 0.0100438 + 3.25466e-05 * age
%
% Gamma parameterization:
%    shape k = (mu/sigma)^2
%    scale θ = sigma^2 / mu
%
% Truncated to ADC range [0, 0.14]
%

%% ---------------- Linear model coefficients -------------------------
bMean = [0.0274639; 0.000233964];     % mean intercept, slope
bSD   = [0.0100438; 3.25466e-05];     % SD intercept, slope

mu_age_raw    = @(a) bMean(1) + bMean(2).*a;
sigma_age_raw = @(a) bSD(1);%   + bSD(2).*a;

sigma_floor = 1e-6;
sigma_age   = @(a) max(sigma_age_raw(a), sigma_floor);
mu_age      = @(a) mu_age_raw(a);

% Gamma parameters
shape_age = @(a) (mu_age(a)./sigma_age(a)).^2;
scale_age = @(a) (sigma_age(a).^2)./mu_age(a);

adc_min = 0;
adc_max = 0.14;

%% ---------------- Build x-axis for output PDF -----------------------
x_adc = linspace(adc_min, adc_max, 400);

k     = shape_age(age);
theta = scale_age(age);

pdf_ref_trunc = truncatedGammaPDF(x_adc, k, theta, adc_min, adc_max);

end  % END OF MAIN FUNCTION




function pdf_t = truncatedGammaPDF(x, k, theta, L, U)
    x = max(min(x,U),L);
    base_pdf = gampdf(x, k, theta);

    cdf_L = gamcdf(L, k, theta);
    cdf_U = gamcdf(U, k, theta);
    Z = max(cdf_U - cdf_L, eps);  % normalizing constant

    pdf_t = base_pdf ./ Z;
    pdf_t(x < L | x > U) = 0;
end
