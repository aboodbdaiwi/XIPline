% Ventilation, CPIR, FLORET
% Set constant params
fov = 375; %use min fov expected since shrinking will likely lead to grad/slew/pns limits; using 375 to match philips
voxel_size = 3.9062; %using ~4 to match philips
mtx = ceil(fov/voxel_size); % determine matrix for desired voxel size 
hub = [3 0.4]; % 3 hubs with 36deg alpha
fname = true; % create name
gmax = 21; % 70 = T1 limits from config- reduced to limit gradient heating limitations
smax = 140; % 150 = T1 limits from config, use lower to enable normal mode scanning
samp = [1.2, 1]; % slight increased sampling along arms, fully sample arms
balanced = false;
calc_dcf = true;
plt = true;
gdt = 4;
do_rot_file = false; % save all trajs
bonus_spec_xe = 10;
bonus_spec_h = 0;

% Design xenon acquisition
nucleus = 'xe129';
order = [2, 21]; % use skip order to minimize change in eddy currents
tacq = 7.6; % philips was 897 arms at 7.42 ms, match arms for similar scan duration
[k_xe,dcf_xe,t_xe,grad_xe,ind_xe,out_xe] = design_floret(fov,mtx,tacq,hub,fname,gmax,smax,...
                 samp,nucleus,balanced,calc_dcf,plt,gdt,order,do_rot_file,bonus_spec_xe);
% Determine file
files = dir("floret_xe129_fov*");
files = struct2table(files);
files = sortrows(files,'date');
% Analyse acquisition
[snr_xe,lw_xe,sz_xe,dur_xe,mdcf_xe,psf_vol_xe] = analyse_trajectory(files.name{1},18,...
            [],[],[],[],fov);

% Design proton acquisition
nucleus = 'h1';
gmax = 70; % 70 = T1 limits from config
smax = 150; % T1 limits from config
order = [2, 8]; % use skip order to minimize change in eddy currents
tacq = 1.49; % philips was 2340 arms at 1.25 ms, match arms for similar scan duration
[k_h,dcf_h,t_h,grad_h,ind_h,out_h] = design_floret(fov,mtx,tacq,hub,fname,gmax,smax,...
                 samp,nucleus,balanced,calc_dcf,plt,gdt,order,do_rot_file,bonus_spec_h);
% Determine file
files = dir("floret_h1_fov*");
files = struct2table(files);
files = sortrows(files,'date');
% Analyse acquisition
[snr_h,lw_h,sz_h,dur_h,mdcf_h,psf_vol_h] = analyse_trajectory(files.name{1},325,...
            [],[],[],[],fov);



%DESIGN_FLORET  Design a 3D floret trajectory
%[k,dcf,t,grad,ind,out] = design_floret(fov,mtx,tacq,hub,fname,gmax,smax,...
%                 samp,nucleus,balanced,calc_dcf,plt,gdt,order,do_rot_file)
%
%        fov  Field of view                                  [mm]
%        mtx  #pixels (Cartesian resolution after gridding)
%       tacq  Duration of single spiral arm                  [ms]       (5)
%        hub  [#hubs alpha0]                                        ([1 1])
%             #hubs=1-3  -> orthogonal rotations of base hub
%             alpha0=0-1 -> max angle of base hub: 1=90°=full sampling
%      fname  Filename: if given->write_ak_wav 
%             if logical true: generate name                        (false)
%       gmax  Max gradient amp                               [mT/m]    (33)
%       smax  Max slew rate                                  [T/m/s]  (120)
%       samp  Samplig factor: >1 over, <1 under (freq,arms)       ([1.2,1])
%    nucleus  Nucleus                                                ('1H')
%   balanced  Balancing gradient area                               (false)
%   calc_dcf  Calculate density compensation function                (true)
%        plt  Plotting                                               (true)
%        gdt  Gradient update time                           [us]       (4)
%      order  Trajectory order:                                         (1)
%             0=consecutive,1=random,2=ll:order(2):end
%do_rot_file  Save hub rotations into phi and theta files           (false)
%
%          k  k-space trajectory [1,#samples*#intlv,3] - [-0.5..0.5]
%        dcf  Density compensation function [1,#samples*#intlv]
%          t  Time               [#intlv,#samples]
%       grad  Gradient waveform [T/m] with dt=4us
%             [#pts/intlv,#intlv,3]
%        ind  Index (2nd dim) for k-space points on spiral (excl ramps,etc)
%        out  Output structure of wrt_wavs
%
% Litrature: mrm66_39_2011_pipe.pdf
% 12/2022  Rolf Schulte