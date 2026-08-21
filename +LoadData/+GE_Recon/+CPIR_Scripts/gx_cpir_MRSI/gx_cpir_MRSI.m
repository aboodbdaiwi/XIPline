% Gas exchange, CPIR, MRSI
% Set constant params
tr_est = 10; %ms
voxel_size = 10;
fovs={[400, 400, 0.8*400],...
    [350, 350, 0.8*350],...
    [300, 300, 0.8*300],...
    [250, 250, 0.8*250]};
mtxs=cellfun(@(v) v./voxel_size,fovs,'UniformOutput',false);
scheme=4;
ktraj=3;
bw=20000;
npts=88;
gmax=70;% T1 limits from config
smax=150;% T1 limits from config
samps={[0.75,0.16],[0.75,0.18],[0.75,0.21],[0.75,0.24]}; % set so bh ~ 14,12,10,8
nucleus='129Xe';
npre=0;
verb=1;
crush=0;
ncok=[50 10]; %10 dummies at higher flip for sat
do_vap=false;
max_exc=0;

for f=1:length(fovs)
    fname=true;
    fov=fovs{f};
    mtx=round(mtxs{f});
    samp=samps{f};
    [ksca,dcf,t,grad,out,fname] = design_mrsi(fov,mtx,scheme,ktraj,bw,...
        npts,gmax,smax,fname,samp,nucleus,npre,verb,crush,ncok,do_vap);

    if (out.N.intl>max_exc)
        max_exc=out.N.intl;
    end

    fprintf('Estimated scan duration: %g seconds for %s\n',tr_est*out.N.intl*1d-3,fname);

    [snr,lw,sz,dur,mdcf,psf_vol] = analyse_trajectory([fname,'.mat']);
end

% Account for downstream magnetization
fa = acosd(1-(tr_est/249)); %fa tr equivalence to match 250ms TR90
sat_fas = linspace(3*fa,fa,ncok(2));
fas = [sat_fas,repmat(fa,max_exc,1)'];
filename = sprintf('vap_flip_gx_%gdeg_tr%g_%dsat',fa,tr_est,ncok(2));
filename = strrep(filename,'.','p');
filename = [filename, '.fdl'];
write_fdl(filename,fas,'flip'); %use desired (TR90~249ms) flip angle for acquisition 


%DESIGN_MRSI  Design trajectory for MRSI with arbitrary sampling patterns
% [k,dcf,t,grad,out,fname] = design_mrsi(fov,mtx,scheme,ktraj,bw,npts,...
%                 gmax,smax,fname,samp,nucleus,npre,verb,crush,ncok,do_vap)
%     fov  Field of view [x,y,z]                           [mm]
%          anisotropic FOV: reference: fov(1)=fov in freq encode dir
%     mtx  Matrix size: [x,y,z]->dim
%  scheme  Sampling scheme                                              (1)
%          1=uniform rectangular
%          2=uniform circular
%          3=accumulation-weighted Hanning
%          4=density-weighted Hanning
%   ktraj  Spatial k-space trajectory                                   (1)
%          1=consecutively
%         -1=reverse consecutively
%          2=random
%          3=centre to outside
%          4=outside to centre
%      bw  Spectral bandwidth (full)                       [Hz]      (5000)
%    npts  Spectral acquisition points                               (1024)
%    gmax  Maximum gradient amplitude                      [mT/m]      (33)
%    smax  Maximum slew rate                               [T/m/s]    (120)
%   fname  filename (if given->write_ak_wav) (true=generate)        (false)
%    samp  Sampling factor: 
%          scheme=3: #averages at k=0                                   (6)
%          scheme=4: >1=over, <1=under-sampling; 
%                    directions=[radial,circumferential]            ([1,1])
% nucleus  Nucleus                                                   ('1H')
%    npre  #acq points before phase encoding gradient                   (0)
%    verb  1 -> print output + plot gradients & k-space                 (1)
%          2 -> also run analyse_trajectory.m
%   crush  Crusher: 0=off                                               (0)
%                   1=same dir as phase encode;crush(2)=factor kmax
%                   2=refocusing (balanced)
%    ncok  Sample centre of k-space every ncok'th acq. (0=off)      ([0 0])
%          ncok(2)=#dummy scans acquired
%  do_vap  Save static trajectory and scale gradients via vap_agXXX (false)
%
%       k  k-space trajectory        [1,#exc,direction]
%     dcf  density compensation function
%       t  time
%    grad  gradient trajectory       [#pts,#exc,direction]
%          directions: x+y+z
%     out  output structure of write_ak_wav.m
%
%  Literature: mrm50_1266_2003_greiser.pdf
% 11/2024 Rolf Schulte


