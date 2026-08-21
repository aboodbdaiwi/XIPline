
fov = [400 400];
mtx = [64 64];
bw = 1.666667e+04;
gmax = 40;
smax = 150;

fname = true;
blncd = false;
nuc = '129Xe';
sampfac = 1;
dda = 0;
pe_sort = 1;    % center out
gdt2 = 4;       % gradient update time
bvalues = [0 12];
tau = 300e-6;   % rise time
del = 3.2e-3;   % single lobe duration
DEL = del;  % time to second lobe

[grad,ind] = design_diffusion_weighted_cart(fov,mtx,bw,gmax,smax,fname,...
                    blncd,nuc,sampfac,dda,pe_sort,gdt2,[0 12],tau,del,DEL);