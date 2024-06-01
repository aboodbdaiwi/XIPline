function traj = siemensspiralcoords(twix_obj,Delay)

if nargin <2
    Delay = [0 0 0];
end

%% Function to calculate trajectories from information passed in the BNI txt file
% Largely copied from the SpiralCoords node in GPI
% This only calculates trajectories for FLORET (i.e. not good for
% other spiral sequences, 3D-radial, stack of stars, etc... it may be 
% useful to add to this and have that functionality in the long run, but
% keep it simple for now.
%
% Notes: need to compile bnispiralgen as mex via e.g. "mex bnispiralgen.cpp"
%
%delay    = trajectory delay. size (1,3 array) -> [x y z] for non-isotopric
%           delays; if size == 1 assumes isotrpic delays
%           If not supplied, zero is assumed ([0 0 0])
%paramfile = name of BNI_RC_*.txt file that contains values of importance
%           If not specified, user is directed to select the file
%=========================================================================

%% PJN - Need to do this smarter, but here is one method to pass the correct parameters 
%Add FOV - in mm
%Add Imaging Dwell Time - Needs to be in us
%Add Resolution (FOV/Matrix) mm
%Add SType -0 for 2D spirals, 3 for FLORET
%Add NumArms
%Add Hubs - set to whatever if doing 2D
%Add Alpha
%Add [USType US0 US1 USR]
%Add BW Scale
%Order (0 for linear, 1 for halton randomized, 2 for golden means)
%Add Delay [x y z] (in us)

%% Pull all needed parameters from twix object
ImSize = twix_obj.hdr.Dicom.lBaseResolution;
Dwell = twix_obj.hdr.MeasYaps.sRXSPEC.alDwellTime{1,1}/1000/2;
FOV = twix_obj.hdr.Config.ReadFoV;
Dim = twix_obj.hdr.Dicom.tMRAcquisitionType;
if contains(Dim,'3D')
    SType = 3;
else
    SType = 0;
end
NSpi = twix_obj.hdr.MeasYaps.sWipMemBlock.adFree{1};
US0_Loc = 9;
US1_Loc = 10;
USR_Loc = 11;
UType_Loc = 8;
US0 = twix_obj.hdr.MeasYaps.sWipMemBlock.adFree{US0_Loc};
US1 = twix_obj.hdr.MeasYaps.sWipMemBlock.adFree{US1_Loc};
USR = twix_obj.hdr.MeasYaps.sWipMemBlock.adFree{USR_Loc};
UType = twix_obj.hdr.MeasYaps.sWipMemBlock.alFree{UType_Loc};
if isempty(UType)
    UType = 0;
end
USamp = [UType US0 US1 USR];
Ord_Loc = 5;
Seq = twix_obj.hdr.Config.SequenceFileName;
if contains(Seq,'DIFF')
    Ord_Loc = 4;
end
Ordering = twix_obj.hdr.MeasYaps.sWipMemBlock.alFree{Ord_Loc};
Alpha_Loc = 3;
Hubs_Loc = 2;
Alpha = twix_obj.hdr.MeasYaps.sWipMemBlock.adFree{Alpha_Loc};
Hubs = twix_obj.hdr.MeasYaps.sWipMemBlock.alFree{Hubs_Loc};
%edit for XA20A
try
    Nucleus = twix_obj.hdr.Spice.ResonantNucleus;
catch
    if ~strcmp(twix_obj.hdr.Dicom.TransmittingCoil,'129Xe_Vest')
        Nucleus = '1H';
    else
        Nucleus = '129Xe';
    end
end
Res = FOV/ImSize;


%% Extract values of interest from Params for spirals
%convert units to ms, kHz, m, mT
dwell = 0.001 * Dwell;
xdely = 0.001 * Delay(1);
ydely = 0.001 * Delay(2);
zdely = 0.001 * Delay(3);
mslew = 100; %PJN - Hardcode to what is on scanner
mgrad = 22;  %PJN - Hardcode to what is on scanner
if strcmp(Nucleus,('1H'))
    gamma = 42.575575; %PJN Hardcode to what is on scanner - This is proton..
else
    gamma = 11.777953; %-xenon
end
fovxy = 0.001 * FOV;
fovz = 0.001 * FOV;
resxy = 0.001 * Res;
resz = 0.001 * Res;
stype = SType;
narms = NSpi;
taper = 0; % PJN - I don't think that I use this in my implementation Param.spTAPER;
us_0 = USamp(2);%Param.spUS0;
us_1 = USamp(3);%Param.spUS1;
us_r = USamp(4);%Param.spUSR;
utype = USamp(1);%Param.spUSTYPE;
mgfrq = 0;%PJN - Not used in my implementationParam.spMGFRQ;
t2mch = 0; %PJN - Not used in my implementation Param.spT2MATCH;
slper = 0; %PJN - currently hardcoded to 0 on the scannerParam.spSLOP_PER;
gtype = 1;% %PJN - currently hardcoded to 1 on the scanner Param.spGTYPE;
spinout = 0;%PJN Not used in my implementationParam.spinout;
if spinout > 1
    numCalPnts = 0;
else
    numCalPnts = 0;
end
%numREADPTS = int32(Param.spREADPTS);


%% Extract values of interest from Params for FLORET
res_sphere = 1;%(pi / 6) ^ (1 / 3);
try
    hubs = Hubs;
catch
    hubs = 3;
end

try
    alpha0 = pi * Alpha / 180.0;
catch
    alpha0 = 0.25 * pi;
end

try
    % Hardcode rebin to 1... I changed the way I reordered trajectories, so
    % easier just to reorder in matlab than to go through mex
    rebin = 0;
    if rebin == 1
        nspiral = round(narms * alpha0 * (fovz / (resz * res_sphere)));
        binfact = round(nspiral / 34.0);
        nspiral = round(34 * binfact);
        narms = nspiral / (alpha0 * (fovz / (resz * res_sphere)));
    end
catch
    rebin = 0;
end
    

%% Compute
%call c/c++ functions
[coords,~] = Recon.bnispiralgen(dwell, xdely, ydely, zdely, mslew, mgrad,...
                               gamma, fovxy, fovz, resxy, resz, stype,...
                               narms, taper, hubs, alpha0, rebin, us_0,...
                               us_1, us_r, utype, mgfrq, t2mch, slper,...
                               gtype, spinout, numCalPnts);
                           
%% Debug Gradients
% NGrad = zeros(size(Grads,1),size(Grads,2)+50,size(Grads,3));
% NGrad(:,51:end,:) = Grads;

% G1 = reshape(NGrad(1,:,:),1,[])';
% G2 = reshape(NGrad(2,:,:),1,[])';
% G3 = reshape(NGrad(3,:,:),1,[])';
% 
% figure('Name','Gradients')
% subplot(3,1,1)
% plot(G1)
% 
% subplot(3,1,2)
% plot(G2)
% 
% subplot(3,1,3)
% plot(G3)

%% Remove extra hubs if single hub
if hubs == 1
    coords = coords(:,:,:,1);
    traj = coords;
else
    traj = coords(:,:,:,1);
    for i = 2:hubs
        traj = cat(3,traj,coords(:,:,:,i));
    end
end

if Ordering == 1
    vals = zeros(1,size(traj,3));
    for i = 1:size(traj,3)
        vals(i) = Halton_rand(i-1,2);
    end
    [~,sort_ind] = sort(vals);
    traj = traj(:,:,sort_ind);
elseif Ordering == 2 && SType == 3
    seed = 0.46557123;
    nums = 0:(size(traj,3)-1);
    vals = (nums*seed - floor(nums*seed));
    [~,sort_ind] = sort(vals);
    traj = traj(:,:,sort_ind);
elseif Ordering == 2 && SType == 0
    seed = pi*(3-sqrt(5));
    nums = 0:(size(traj,3)-1);
    vals = (nums*seed - floor(nums*seed));
    [~,sort_ind] = sort(vals);
    traj = traj(:,:,sort_ind);
end

traj = -traj;
%% Confirm match 
% nsamp = size(coords,2);
% %match data points from parameter file or real data
% if (numREADPTS > 0) && (nsamp > numREADPTS)
%     coords = coords(:,1:numREADPTS,:,:);
% end


%traj = coords;




  
