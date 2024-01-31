
function ADC_Model_1D=ADC_Model_Human_1D(ky,num_bavlue)
%   Inputs:
%      
%   Outputs:
%                   
%   Package: https://github.com/aboodbdaiwi/HP129Xe_Analysis_App
%
%   Author: Abdullah S. Bdaiwi
%   Work email: abdullah.bdaiwi@cchmc.org
%   Personal email: abdaiwi89@gmail.com
%   Website: https://www.cincinnatichildrens.org/research/divisions/c/cpir
%
%   Please add updates at the end. Ex: 3/10/24 - ASB: updated flip angle caluclation for spiral diffusion

    formatSpec_y = '    for (i in 1:%d){';
    formatSpec_b = '                for (l in 1:%d){';
    str_y = sprintf(formatSpec_y,ky);
    str_b = sprintf(formatSpec_b,num_bavlue);
    ADC_Model_1D = fopen('ADC_Model_Human_1D.txt','wt');
    fprintf(ADC_Model_1D, '%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s', 'model <- function()'...
                             ,'{'...
                             ,    str_y...
                             ,                str_b...
                             ,'                 lambda[i, l] <- ( pow(nu0[i],2)* exp(-2*alpha[i]*b[l]) )/(2*pow(sigma[l],2)) '...
                             ,'                  p[i, l] ~ dpois(lambda[i, l])'...
                             ,'                  kk[i, l] <- 2 * p[i, l] + 2'...
                             ,'                  M[i, l] ~ dchisqr(kk[i, l])'...
                             ,'                }'...
                             ,'                nu0[i] ~ dunif(0.0001, 200)'...
                             ,'                alpha[i] ~ dunif(1.00000E-10, 0.14)'...
                             ,'    }'...
                             ,'}');
    fclose(ADC_Model_1D);
end












