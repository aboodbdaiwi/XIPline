% -GE CONFIDENTIAL-
% Type: Source Code
%
% Copyright (c) 2022, GE Healthcare
% All Rights Reserved
%
% This unpublished material is proprietary to GE Healthcare. The methods and
% techniques described herein are considered trade secrets and/or
% confidential. Reproduction or distribution, in whole or in part, is
% forbidden except by express written permission of GE Healthcare.

function  a = read_grad_header(my_file, rdbm_rev)% RDB revision 24.000
% RDBM revision 24.000
if rdbm_rev == 24.000 
  for id1 = 1 : 512
    for id2 = 1 : 3
        a.diffusion_grad_amp(id1,id2) = fread(my_file, 1, 'float32');
    end
  end
  a.hoec_bases.hoec_coef = fread(my_file, [3, 56], 'float32')';
  a.hoec_bases.hoec_xorder = fread(my_file, [1,56], 'int32')';
  a.hoec_bases.hoec_yorder = fread(my_file, [1,56], 'int32')';
  a.hoec_bases.hoec_zorder = fread(my_file, [1,56], 'int32')';


end

% RDBM revision 25.001
if rdbm_rev == 25.001 
  for id1 = 1 : 512
    for id2 = 1 : 3
        a.diffusion_grad_amp(id1,id2) = fread(my_file, 1, 'float32');
    end
  end
  a.hoec_bases.hoec_coef = fread(my_file, [3, 56], 'float32')';
  a.hoec_bases.hoec_xorder = fread(my_file, [1,56], 'int32')';
  a.hoec_bases.hoec_yorder = fread(my_file, [1,56], 'int32')';
  a.hoec_bases.hoec_zorder = fread(my_file, [1,56], 'int32')';


end

% RDBM revision 25.002
if rdbm_rev == 25.002 
  for id1 = 1 : 512
    for id2 = 1 : 3
        a.diffusion_grad_amp(id1,id2) = fread(my_file, 1, 'float32');
    end
  end
  a.hoec_bases.hoec_coef = fread(my_file, [3, 56], 'float32')';
  a.hoec_bases.hoec_xorder = fread(my_file, [1,56], 'int32')';
  a.hoec_bases.hoec_yorder = fread(my_file, [1,56], 'int32')';
  a.hoec_bases.hoec_zorder = fread(my_file, [1,56], 'int32')';


end

% RDBM revision 25.003
if rdbm_rev == 25.003 
  for id1 = 1 : 512
    for id2 = 1 : 3
        a.diffusion_grad_amp(id1,id2) = fread(my_file, 1, 'float32');
    end
  end
  a.hoec_bases.hoec_coef = fread(my_file, [3, 56], 'float32')';
  a.hoec_bases.hoec_xorder = fread(my_file, [1,56], 'int32')';
  a.hoec_bases.hoec_yorder = fread(my_file, [1,56], 'int32')';
  a.hoec_bases.hoec_zorder = fread(my_file, [1,56], 'int32')';


end

% RDBM revision 25.004
if rdbm_rev == 25.004 
  for id1 = 1 : 512
    for id2 = 1 : 3
        a.diffusion_grad_amp(id1,id2) = fread(my_file, 1, 'float32');
    end
  end
  a.hoec_bases.hoec_coef = fread(my_file, [3, 56], 'float32')';
  a.hoec_bases.hoec_xorder = fread(my_file, [1,56], 'int32')';
  a.hoec_bases.hoec_yorder = fread(my_file, [1,56], 'int32')';
  a.hoec_bases.hoec_zorder = fread(my_file, [1,56], 'int32')';


end

% RDBM revision 26.000
if rdbm_rev == 26.000 
  for id1 = 1 : 512
    for id2 = 1 : 3
        a.diffusion_grad_amp(id1,id2) = fread(my_file, 1, 'float32');
    end
  end
  a.hoec_bases.hoec_coef = fread(my_file, [3, 56], 'float32')';
  a.hoec_bases.hoec_xorder = fread(my_file, [1,56], 'int32')';
  a.hoec_bases.hoec_yorder = fread(my_file, [1,56], 'int32')';
  a.hoec_bases.hoec_zorder = fread(my_file, [1,56], 'int32')';


end

% RDBM revision 26.001
if rdbm_rev == 26.001 
  for id1 = 1 : 512
    for id2 = 1 : 3
        a.diffusion_grad_amp(id1,id2) = fread(my_file, 1, 'float32');
    end
  end
  a.hoec_bases.hoec_coef = fread(my_file, [3, 56], 'float32')';
  a.hoec_bases.hoec_xorder = fread(my_file, [1,56], 'int32')';
  a.hoec_bases.hoec_yorder = fread(my_file, [1,56], 'int32')';
  a.hoec_bases.hoec_zorder = fread(my_file, [1,56], 'int32')';


end

% RDBM revision 26.002
if rdbm_rev == 26.002 
  for id1 = 1 : 512
    for id2 = 1 : 3
        a.diffusion_grad_amp(id1,id2) = fread(my_file, 1, 'float32');
    end
  end
  a.hoec_bases.hoec_coef = fread(my_file, [3, 56], 'float32')';
  a.hoec_bases.hoec_xorder = fread(my_file, [1,56], 'int32')';
  a.hoec_bases.hoec_yorder = fread(my_file, [1,56], 'int32')';
  a.hoec_bases.hoec_zorder = fread(my_file, [1,56], 'int32')';


end

% RDBM revision 27.000
if rdbm_rev == 27.000 
  for id1 = 1 : 1024
    for id2 = 1 : 3
        a.diffusion_grad_amp(id1,id2) = fread(my_file, 1, 'float32');
    end
  end
  a.hoec_bases.hoec_coef = fread(my_file, [3, 56], 'float32')';
  a.hoec_bases.hoec_xorder = fread(my_file, [1,56], 'int32')';
  a.hoec_bases.hoec_yorder = fread(my_file, [1,56], 'int32')';
  a.hoec_bases.hoec_zorder = fread(my_file, [1,56], 'int32')';


end

% RDBM revision 27.001
if rdbm_rev == 27.001 
  for id1 = 1 : 1024
    for id2 = 1 : 3
        a.diffusion_grad_amp(id1,id2) = fread(my_file, 1, 'float32');
    end
  end
  a.hoec_bases.hoec_coef = fread(my_file, [3, 56], 'float32')';
  a.hoec_bases.hoec_xorder = fread(my_file, [1,56], 'int32')';
  a.hoec_bases.hoec_yorder = fread(my_file, [1,56], 'int32')';
  a.hoec_bases.hoec_zorder = fread(my_file, [1,56], 'int32')';


end

% RDBM revision 28.000
if rdbm_rev == 28.000 
  for id1 = 1 : 1024
    for id2 = 1 : 3
        a.diffusion_grad_amp(id1,id2) = fread(my_file, 1, 'float32');
    end
  end
  a.hoec_bases.hoec_coef = fread(my_file, [3, 56], 'float32')';
  a.hoec_bases.hoec_xorder = fread(my_file, [1,56], 'int32')';
  a.hoec_bases.hoec_yorder = fread(my_file, [1,56], 'int32')';
  a.hoec_bases.hoec_zorder = fread(my_file, [1,56], 'int32')';


end

% RDBM revision 28.002
if rdbm_rev == 28.002 
  for id1 = 1 : 1024
    for id2 = 1 : 3
        a.diffusion_grad_amp(id1,id2) = fread(my_file, 1, 'float32');
    end
  end
  a.hoec_bases.hoec_coef = fread(my_file, [3, 56], 'float32')';
  a.hoec_bases.hoec_xorder = fread(my_file, [1,56], 'int32')';
  a.hoec_bases.hoec_yorder = fread(my_file, [1,56], 'int32')';
  a.hoec_bases.hoec_zorder = fread(my_file, [1,56], 'int32')';


end

% RDBM revision 28.003
if rdbm_rev == 28.003 
  for id1 = 1 : 1024
    for id2 = 1 : 3
        a.diffusion_grad_amp(id1,id2) = fread(my_file, 1, 'float32');
    end
  end
  a.hoec_bases.hoec_coef = fread(my_file, [3, 56], 'float32')';
  a.hoec_bases.hoec_xorder = fread(my_file, [1,56], 'int32')';
  a.hoec_bases.hoec_yorder = fread(my_file, [1,56], 'int32')';
  a.hoec_bases.hoec_zorder = fread(my_file, [1,56], 'int32')';


end

% RDBM revision 30.000
if rdbm_rev == 30.000 
  for id1 = 1 : 1024
    for id2 = 1 : 3
        a.diffusion_grad_amp(id1,id2) = fread(my_file, 1, 'float32');
    end
  end
  a.hoec_bases.hoec_coef = fread(my_file, [3, 56], 'float32')';
  a.hoec_bases.hoec_xorder = fread(my_file, [1,56], 'int32')';
  a.hoec_bases.hoec_yorder = fread(my_file, [1,56], 'int32')';
  a.hoec_bases.hoec_zorder = fread(my_file, [1,56], 'int32')';


end

