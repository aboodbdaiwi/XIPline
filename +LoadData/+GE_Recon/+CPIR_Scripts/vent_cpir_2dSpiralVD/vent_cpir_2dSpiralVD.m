% Ventilation, CPIR, 2D, Variable Density Spiral
% Set params
fovs = [400, 350, 300, 250, 200]; %use min fov expected since shrinking will likely lead to grad/slew/pns limits; alternatively (here) scale gmax/smax 
voxel_size = 3; %using ~5 to match philips
mtxs = ceil(fovs./voxel_size); % determine matrix for desired voxel size 
tacq_xe = 7;
tacq_h = 2;
vardens_xe = [2 0.5 1 2]; % hanning var dens with undersmapling starting at 0.5kmax and going to 2 fold undersampling at edge
vardens_h = [0, 0, 1, 1]; % no var dens
fname = true; % create name
gmax = 70; % T1 limits from config
smax_xe = 120; % 150 T1 limits from config; using 120 so PNS<=80 for normal mode
smax_h = 150; % 150 T1 limits from config
samp_xe = [1.2 2]; % slight oversample in read, 2 fold oversample for center of kspace
samp_h = [1.2 1]; % slight oversample in read
balanced = false;
calc_dcf = true;
plt = true;
gdt = 4;
order = [2 4]; % skip every 4
do_rot_file = false; % save all trajs
bonus_spec_xe = 2;
bonus_spec_h = 0;
    
for acq = 1:length(fovs)
    fov = fovs(acq);
    mtx = mtxs(acq);
    [k,dcf,t,grad,ind,out] = design_spiral_bni(fov,mtx,tacq_xe,vardens_xe,fname,gmax,smax_xe,...
                samp_xe,'129xe',balanced,calc_dcf,plt,gdt,order,do_rot_file,bonus_spec_xe);
    [k,dcf,t,grad,ind,out] = design_spiral_bni(fov,mtx,tacq_h,vardens_h,fname,gmax,smax_h,...
                samp_h,'1h',balanced,calc_dcf,plt,gdt,order,do_rot_file,bonus_spec_h);
end

%DESIGN_SPIRAL_BNI  Design a 2D spiral trajectory from BNI code
%[k,dcf,t,grad,ind,out] = design_spiral_bni(fov,mtx,tacq,vardens,fname,gmax,smax,...
%                 samp,nucleus,balanced,calc_dcf,plt,gdt,order,do_rot_file)
%
%        fov  Field of view                                  [mm]
%        mtx  #pixels (Cartesian resolution after gridding)
%       tacq  Duration of single spiral arm                  [ms]       (5)
%    vardens  Variable density (ustype,us_0,us_1,us_r)            (0,1,1,1)
%      fname  Filename: if given->write_ak_wav 
%             if logical true: generate name                        (false)
%       gmax  Max gradient amp                               [mT/m]    (33)
%       smax  Max slew rate                                  [T/m/s]  (120)
%       samp  Sampling factor: >1 over, <1 under (freq,arms)      ([1.2,1])
%    nucleus  Nucleus                                                ('1H')
%   balanced  Balancing gradient area                               (false)
%   calc_dcf  Calculate density compensation function                (true)
%        plt  Plotting                                               (true)
%        gdt  Gradient update time                           [us]       (4)
%      order  Trajectory order:                                         (1)
%             0=consecutive,1=random,2=ll:order(2):end
%do_rot_file  Write single static spiral into gradient 
%             waveform file & rotate via vap_phiXX.fdl (false)
%
%          k  k-space trajectory [1,#samples*#intlv,2] - [-0.5..0.5]
%        dcf  Density compensation function [1,#samples*#intlv]
%          t  Time               [#intlv,#samples]
%       grad  Gradient waveform [T/m] with dt=4us
%             [#pts/intlv,#intlv,3]
%        ind  Index (2nd dim) for k-space points on spiral (excl ramps,etc)
%        out  Output structure of wrt_wavs