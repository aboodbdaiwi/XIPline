% Diffusion, CPIR, 2D, Cartesian, 4bs
% Set cartesian params
fovs = [400, 350, 300, 250, 200]; %use min fov expected since shrinking will likely lead to grad/slew/pns limits; alternatively (here) scale gmax/smax 
voxel_size = 5; %using ~5 to match philips
mtxs = ceil(fovs./voxel_size); % determine matrix for desired voxel size 
bw = 20d3; %20kHz
gmax = 70; % T1 limits from config
smax = 100; % 150 T1 limits from config; using 100 so PNS<=80 for normal mode
fname = true; % create name
blncd = false;
nucleus = '129Xe';
sampfac = 1; %1x read oversampling
dda = 0; % no dummys
pe_sort = 1; % centric
gdt2 = 4;
do_vap = false; % save all trajs
inv_grad = false;

%Diffusion parameters
b_vals = [0, 10, 20, 30];
diff_time = 5; 
diff_gap = 0;
    
for acq = 1:length(fovs)
    fov = fovs(acq);
    mtx = mtxs(acq);
    [grad,ind] = design_cart_diffusion  ([fov, fov],[mtx, mtx],bw,gmax,smax,fname,blncd,nucleus,...
       sampfac,dda,pe_sort,gdt2,do_vap,inv_grad,...
       b_vals, diff_time, diff_gap);
end

%DESIGN_CART  Design Cartesian readout (1D,2D or 3D)
% [grad,ind] = design_cart(fov,mtx,bw,gmax,smax,fname,blncd,nucleus,...
%    sampfac,dda,pe_sort,gdt2,do_vap,inv_grad)
%                                                       [unit]    (default)
%        fov  Field of view                             [mm]
%        mtx  Matrix size                               [1]
%         bw  Full sampling bandwidth (+-)              [Hz]        ('max')
%       gmax  Maximum gradient strength                 [mT/m]         (33)
%       smax  Maximum slewrate                          [T/m/s]       (120)
%****  fov,mtx,gmax,smax: 1D=[freq]; 2D=[freq phase]; 3D=[freq,phase,phase]
%             dimensions=max(length(fov,mtx,gmax,smax))
%      fname  Filename (string, or logical)                         (false)
%      blncd  Balance gradient moments (for SSFP)                    (true)
%    nucleus  Nucleus                                                ('1H')
%    sampfac  Sampling factor (>1 -> oversampling along freq)           (2)
%        dda  Dummy scans (acquired)                                    (0)
%    pe_sort  Sorting of phase encoding:                                (0)
%             0=sequential, 1=centre-out
%       gdt2  Gradient update time for storing          [us]            (4)
%             must be multiple of 4
%     do_vap  Write agy+z scalings into vap file; use static grad wf (true)
%   inv_grad  Invert gradients                                      (false)
%
%       grad  Gradient trajectory                       [T/m/s]
%             (points,1,direction) 
%             direction: 1=x=freq, 2=y=phase1, 3=z=phase2
%        ind  Index for acquired data
%
%  9/2024 Rolf Schulte