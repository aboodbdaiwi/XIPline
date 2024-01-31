
function fermi2 = fermi_filter_2D_AB(ky,kx,Wf)
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
%   Please add updates at the end. Ex: 3/10/24 - ASB: update .... 

    %GE default for Wf is 10
%      Centering corrected to window_size/2 on 10/16/17 by ZIC.
     Fermi_xx = LoadData.fermi_filter_1D(kx,Wf);
     Fermi_xx=repmat(Fermi_xx, ky,1);
     Fermi_xx=rot90(reshape(Fermi_xx, [kx, ky]));
    % figure; mesh(Fermi_xx)

     Fermi_yy = LoadData.fermi_filter_1D(ky,Wf);
     Fermi_yy=repmat(Fermi_yy, kx,1);
     Fermi_yy=(reshape(Fermi_yy, [ky, kx]));
    % figure; mesh(Fermi_yy)
    fermi2=Fermi_xx.*Fermi_yy; 
    % figure; mesh(fermi2)
end

  