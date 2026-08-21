fov = [480 480];
mtx = [64 32];

bw = 1920;
kdt = 1/1920*1e6;           % [us]
fname = true;
gmax = 40;                  % [mT/m]
smax = 130;                 % [T/m/s]
nucleus = '129Xe';

sampfac = 1;
dda = 0;
pe_sort = 0;    % sequential

[grad,ind] = design_cart(fov,mtx,bw,gmax,smax,true,false,...
    nucleus,sampfac,dda,pe_sort,[]);