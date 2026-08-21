% Gas Exchange, CPIR, 1-point Dixon
% Set params
fov = 325; %match philips
mtx = 56; %match philips
kdt = 18; %match philips (~1.1ms readout) 
fname = true;
gmax = 70; % T1 limits from config
smax = 150; % 150 = T1 limits from config
nucleus = '129xe';
dim = 3;
typ = 0;
sampfac = [2, 1000]; %2fold oversample, 1000 spokes/intlv
rewinder = [2, 1.25]; %spoiler pushing out to 2.25 kmax, worse case scenario of 1.25
order = [3, 5]; %golden angle - tiny(5)
do_rot_file = false; % bonus fids require saving all
dda = 0;
flips = [20, 0.5]; %dissolved first, gas second; not match cpir's 3 intlve approach as not necessary
freqs = [7143, 0];  %dissolved first, gas second; match cpir but may be bad for contam
bonus_fids = [1, 1]; %bonus fid for both

[k,dcf,t,grad,ind,fname] = design_interleaved_radial(fov,mtx,kdt,fname,...
    gmax,smax,nucleus,dim,typ,sampfac,rewinder,order,do_rot_file,dda,flips,freqs,bonus_fids);

%% Design matching proton
nucleus = '1h';
kdt = 18;
fname = true;
sampfac = [2, 5000];
order = [3, 5]; %golden angle - tiny(5)
do_rot_file = true;

[k,dcf,t,grad,ind,fname] = design_radial(fov,mtx,kdt,fname,...
    gmax,smax,nucleus,dim,typ,sampfac,rewinder,order,do_rot_file,dda);


%DESIGN_INTERLEAVED_RADIAL  Design isotropic 2D or 3D radial sampling trajectory
% [k,dcf,t,grad,ind,fname] = design_interleaved_radial(fov,mtx,kdt,fname,gmax,smax,...
%            nucleus,dim,typ,sampfac,rewinder,order,do_rot_file,dda)
%                                                     [unit]    (default)
%        fov  Field-of-view                           [mm]
%        mtx  Matrix size                             [1]
%        kdt  Sampling dwell time (1/bw)              [us]      (16)
%      fname  Filename (string, or logical)                     (false)
%       gmax  Maximum gradient strength               [mT/m]    (33)
%       smax  Maximum slewrate                        [T/m/s]   (120)
%    nucleus  Nucleus                                           ('1H')
%        dim  Spatial dimensions (2=2D,3=3D,4=2D+z-PE)          (3)
%        typ  Radial variant:                                   (0)
%             0=centre-out constant; 1=centre-out density-adapted
%             2=full-spoke constant
%    sampfac  Sampling factor (>1 over-, <1 under-sampled)      (1.5 1)
%             1=frequency encoding direction (speed in k-space)
%             2=spokes spacing@kmax=sqrt(opnex/pi) in 3dradial
%               if sampfac(2)>10 -> #spokes=sampfac(2)
%   rewinder  Gradient ramp down method:                        (0)
%             0=fast, 1=return to centre of k-space, 2=crushing
%             2->[2 area-factor(default=1)]
%      order  Spokes ordering:                                  (1)
%             equidistant:1=subsequent,2=random; 3=golden angle
%             if 3->order(2)=type (input to golden_angle.m)     ([... 1])
%                   order(4)=0->tangential correct; 1->equidistant dcf
%             if dim==4->order(3)=ordering in z-PE
%do_rot_file  Write single static spoke into gradient waveform file and 
%             rotate via vap_phiXX.fdl and vap_thetaXX.fdl      (true)
%        dda  Dummy scans to discard                            ([0 1])
%     intlvs  Interleaves of same trajectory                    (1)
%
%          k  k-Space trajectory (normalised to +-0.5)
%             #acq-pts #interleaves #dim
%        dcf  Density compensation function (analytical)
%             assumes equal spokes distribution; some error for 3D mtx<32
%          t  Time                                    [s]
%       grad  Gradient trajectory                     [T/m/s]
%        ind  Index for acquired data
%
% Literature: Handbook of MRI Pulse Sequence; Bernstein, King, Zhou
%             mrm65_1091_2011_konstandin.pdf,mrm62_1565_2009_nagel.pdf
% 11/2022 Rolf Schulte