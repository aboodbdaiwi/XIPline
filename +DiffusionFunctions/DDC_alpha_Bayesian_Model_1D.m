
function DDC_Model_1D=DDC_alpha_Bayesian_Model_1D(ky,num_bvalues)
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
    str_b = sprintf(formatSpec_b,num_bvalues);
    DDC_Model_1D = fopen('DDC_alpha_Bayesian_Model_1D.txt','wt');
    fprintf(DDC_Model_1D, '%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s', 'model <- function()'...
                             ,'{'...
                             ,    str_y...
                             ,                str_b...
                             ,'                 lambda[i, l] <- ( pow(So[i],2)* exp(-2*pow((DDC[i]*b[l]),alpha[i])))/(2*pow(sigma[l],2)) '...
                             ,'                  p[i, l] ~ dpois(lambda[i, l])'...
                             ,'                  kk[i, l] <- 2 * p[i, l] + 2'...
                             ,'                  M[i, l] ~ dchisqr(kk[i, l])'...
                             ,'                }'...
                             ,'                So[i] ~ dunif(0.0001, 200)'...
                             ,'                DDC[i] ~ dunif(1.00000E-10, 0.14)'...
                             ,'                alpha[i] ~ dunif(1.00000E-10, 1)'...
                             ,'    }'...
                             ,'}');
    fclose(DDC_Model_1D);
end

% Please add updates here
% 
% 
% 
% 










