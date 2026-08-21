%% Base Specs/Params
gr_max = 33; %A6 Artist
s_max = 120; %A6 Artist

%% Design Low Res Test UTEs
fov=300;
mtx=64;
nuc='1H';
fname=true; % autoname

% Density-adapted Tiny-Golden-Angle Crusher-Spoiling Radial UTE
% design_radial(fov,mtx,kdt,fname,gmax,smax,nucleus,dim,typ,sampfac,rewinder,order,do_rot_file,dda)
kdt=10; % controls readout duration/bandwidth
dim=3; % kooshball
typ=1; % density adapted radial UTE
ro_samp=1.25; % read-out oversampling
proj_samp=1.0; % projection sampling
rewinder=2; % crusher
order=[3,34,0,0]; % golden angle, tiny angle 34, N/A, no equidistant correction
design_radial(fov,mtx,kdt,fname,gr_max,s_max,nuc,dim,typ,[ro_samp,proj_samp],rewinder,order);
%results in 12868 projs, 16382 is max per rep, 16 reps max = 1 rep of 12868 projs

% FLORET UTE
% linear 1 Hub
%[k,dcf,t,grad,ind,out]=design_floret(fov,mtx,tacq,hub,fname,gmax,smax,samp,nucleus,balanced,calc_dcf,plt,gdt2,order,do_rot_file)
tacq=2; % desired acquisition time in ms
hubs=1;
max_ang=90; % alpha in degrees
ro_samp=1.0;% read-out oversampling
proj_samp=1.0;% projection sampling
balanced=false;
order=0; % linear
design_floret(fov,mtx,tacq,[hubs,max_ang/90],fname,gr_max,s_max,[ro_samp, proj_samp],nuc,balanced,true,true,16,order,false);
%results in 604 projs, 16382 is max per rep, 16 reps max = 1 rep of 604 projs

% linear 3 Hub
%[k,dcf,t,grad,ind,out]=design_floret(fov,mtx,tacq,hub,fname,gmax,smax,samp,nucleus,balanced,calc_dcf,plt,gdt2,order,do_rot_file)
tacq=2; % desired acquisition time in ms
hubs=3;
max_ang=36; % alpha in degrees
ro_samp=1.0;% read-out oversampling
proj_samp=1.0;% projection sampling
balanced=false;
order=0; % linear
design_floret(fov,mtx,tacq,[hubs,max_ang/90],fname,gr_max,s_max,[ro_samp, proj_samp],nuc,balanced,true,true,16,order,false);
%results in 726 projs, 16382 is max per rep, 16 reps max = 1 rep of 726 projs

% linear 3 Hub - vap
%[k,dcf,t,grad,ind,out]=design_floret(fov,mtx,tacq,hub,fname,gmax,smax,samp,nucleus,balanced,calc_dcf,plt,gdt2,order,do_rot_file)
tacq=2; % desired acquisition time in ms
hubs=3;
max_ang=45; % alpha in degrees
ro_samp=1.0;% read-out oversampling
proj_samp=1.0;% projection sampling
balanced=false;
order=0; % linear
design_floret(fov,mtx,tacq,[hubs,max_ang/90],fname,gr_max,s_max,[ro_samp, proj_samp],nuc,balanced,true,true,16,order,true);
%results in 906 projs, 16382 is max per rep, 16 reps max = 1 rep of 906 projs

% interleaved 3 Hub - vap 
%[k,dcf,t,grad,ind,out]=design_floret(fov,mtx,tacq,hub,fname,gmax,smax,samp,nucleus,balanced,calc_dcf,plt,gdt2,order,do_rot_file)
tacq=2; % desired acquisition time in ms
hubs=3;
max_ang=45; % alpha in degrees
ro_samp=1.0;% read-out oversampling
proj_samp=1.033;% projection sampling results in 936 projs which is divisible by 3 (hubs), 13 (interleaves), 1 (reps)
balanced=false;
order=[2, 13]; % skip 13, force no rot file
design_floret(fov,mtx,tacq,[hubs,max_ang/90],fname,gr_max,s_max,[ro_samp, proj_samp],nuc,balanced,true,true,16,order,true);
%results in 936 projs, 16382 is max per rep, 16 reps max = 1 rep of 936 projs

%% Design High Res Test UTEs
fov=200;
mtx=256;
 
% FLORET
% linear 1 Hub
%[k,dcf,t,grad,ind,out]=design_floret(fov,mtx,tacq,hub,fname,gmax,smax,samp,nucleus,balanced,calc_dcf,plt,gdt2,order,do_rot_file)
tacq=2; % desired acquisition time in ms
hubs=1;
max_ang=90; % alpha in degrees
ro_samp=1.0;% read-out oversampling
proj_samp=1.00000636821;% projection sampling
balanced=false;
order=0; % linear
design_floret(fov,mtx,tacq,[hubs,max_ang/90],fname,gr_max,s_max,[ro_samp, proj_samp],nuc,balanced,true,true,16,order,false);
%results in 31408 projs, 16382 is max per rep, 16 reps max = 16 rep of 1963 projs

% linear 3 Hub
%[k,dcf,t,grad,ind,out]=design_floret(fov,mtx,tacq,hub,fname,gmax,smax,samp,nucleus,balanced,calc_dcf,plt,gdt2,order,do_rot_file)
tacq=2; % desired acquisition time in ms
hubs=3;
max_ang=36; % alpha in degrees
ro_samp=1.0;% read-out oversampling
proj_samp=1.001;% projection sampling
balanced=false;
order=0; % linear
design_floret(fov,mtx,tacq,[hubs,max_ang/90],fname,gr_max,s_max,[ro_samp, proj_samp],nuc,balanced,true,true,16,order,false);
%results in 37728 projs, 16382 is max per rep, 16 reps max = 16 rep of 2358 projs

% linear 3 Hub - vap
%[k,dcf,t,grad,ind,out]=design_floret(fov,mtx,tacq,hub,fname,gmax,smax,samp,nucleus,balanced,calc_dcf,plt,gdt2,order,do_rot_file)
tacq=2; % desired acquisition time in ms
hubs=3;
max_ang=45; % alpha in degrees
ro_samp=1.0;% read-out oversampling
proj_samp=1.00051;% projection sampling
balanced=false;
order=0; % linear
design_floret(fov,mtx,tacq,[hubs,max_ang/90],fname,gr_max,s_max,[ro_samp, proj_samp],nuc,balanced,true,true,16,order,true);
%results in 47136 projs, 16382 is max per rep, 16 reps max = 16 rep of 2946 projs

%% Design Neonatal Lung UTE
fov=200;
mtx=256;
nuc='1H';
fname=true; % autoname

% Density-adapted Tiny-Golden-Angle Crusher-Spoiling Radial UTE
% design_radial(fov,mtx,kdt,fname,gmax,smax,nucleus,dim,typ,sampfac,rewinder,order,do_rot_file,dda)
kdt=4; % controls readout duration/bandwidth
dim=3; % kooshball
typ=1; % density adapted radial UTE
ro_samp=1.25; % read-out oversampling
proj_samp=1.0; % projection sampling
rewinder=2; % crusher
order=[3,34,0,0]; % golden angle, tiny angle 34, N/A, no equidistant correction
[k1,dcf1,t1,grad1,ind1,fname1] = design_radial(fov,mtx,kdt,fname,gr_max,s_max,nuc,dim,typ,[ro_samp,proj_samp],rewinder,order);
%results in 205888 projs, 16382 is max per rep, 16 reps max = 16 reps of 12868 projs
%tr ~ 1.7 ms = ~6 minutes

% FLORET UTE
% %[k,dcf,t,grad,ind,out]=design_floret(fov,mtx,tacq,hub,fname,gmax,smax,samp,nucleus,balanced,calc_dcf,plt,gdt2,order,do_rot_file)
% tacq=2; % desired acquisition time in ms
% hubs=3;
% max_ang=90; % alpha in degrees
% ro_samp=1.25;% read-out oversampling
% proj_samp=1.36025472935;% projection sampling - results in 128160 projs which is divisible by 3 (hubs), 89 (interleaves), 16 (reps)
% balanced=false;
% order=[2, 89]; % skip 89, force no rot file
% [k2,dcf2,t2,grad2,ind2,fname2] = design_floret(fov,mtx,2,[hubs,max_ang/90],fname,gr_max,s_max,[ro_samp, proj_samp],nuc,balanced,true,true,4,order,false);
% %results in 128160 projs, 16382 is max per rep, 16 reps max = 16 reps of 8010 projs
% %tr ~ 2.3 ms = ~5 minutes