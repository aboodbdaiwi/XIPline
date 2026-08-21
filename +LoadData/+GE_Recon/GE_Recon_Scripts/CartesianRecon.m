function CartesianRecon(pfilePath)
%% CartesianRecon - Reconstruct 2D or 3D Cartesian K-Space
%
% Copyright (c) 2019, 2025, GE HealthCare
% All Rights Reserved
% 
% This material is proprietary to GE HealthCare. The methods and
% techniques described herein are considered trade secrets and/or 
% confidential. Reproduction or distribution, in whole or in part, is
% forbidden except by express written permission of GE HealthCare.
% GE is a trademark of General Electric Company used under trademark license.
% Resulting outputs are not for diagnostic purposes.
%
% CartesianRecon(pfilePath)
% will reconstruct the Cartesian K-Space in the given pfile. Except for
% 3D scans with ARC, this includes Pfiles from both 2D and 3D
% acquisitions.
%
% Limitations: Parallel imaging, intensity correction

    % Load Pfile
    pfile = GERecon('Pfile.Load', pfilePath);
    zEncoded = pfile.isZEncoded;
    clear pfile;
    
    if zEncoded == 1
        Cartesian3DRecon(pfilePath);
    else
        Cartesian2DRecon(pfilePath);
    end
