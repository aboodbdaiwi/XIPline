function traj = gx_traj(NPts,NPro,twix)

%I'll have to come up with something slightly different for proton, but for
%xenon, the gradient shape and ADC parameters are:

% Parameters for xenon
if nargin < 3
    RU = 100;
    FT =  3000;
    Dw = 10;
    Seq_Name = 'Unknown';
else
    Seq_Name = twix.hdr.Config.SequenceFileName;
    if str2double(Seq_Name((end-3):end)) >= 2102
        try
            RU = twix.hdr.MeasYaps.sWipMemBlock.alFree{8};
        catch
            RU = twix.hdr.MeasYaps.sWiPMemBlock.alFree{8};
        end
    else
        try
            RU = twix.hdr.MeasYaps.sWipMemBlock.alFree{9};
        catch
            RU = twix.hdr.MeasYaps.sWiPMemBlock.alFree{9};
        end
    end
    if isempty(RU)
        RU = 100;
    end
    Dw = twix.hdr.MeasYaps.sRXSPEC.alDwellTime{1,1}/1000; %Dwell is in nanoseconds in twix. Need us    
    FT = 10000;%Doesn't matter how long flat time is - it will only be sampled for the correct number of points.
end
% GMax = 10.613;
% Res = .0063;
% gamma = 11.777;

%% Create Gradient Shape
Ramp_Shape = linspace(0,1,RU);
Flat_Shape = linspace(1,1,FT);

Grad_Shape = [Ramp_Shape Flat_Shape];

%% Resample at the dwell time
Dw_Pts = 0:Dw:((NPts-1)*Dw);

re_Grad_Shape = interp1((1:length(Grad_Shape)),Grad_Shape,Dw_Pts);

re_Grad_Shape(isnan(re_Grad_Shape)) = 0;
%% Integrate to get k-space location for one radial projection
k_loc = cumtrapz(Dw_Pts,re_Grad_Shape);

%k_loc = k_loc*gamma/1000;

%% Get rotation angles for projections

Angs = twix.image.iceParam(1:2,:);
necho = 1;
while Angs(1,1) == Angs(1,necho)
    necho = necho + 1;
end
necho = necho - 1; %subtract off the last one

Angs = Angs(:,1:necho:end);
Angs(:,(Angs(1,:)==0 & Angs(2,:) ==0)) = [];
Angs = Angs'/10000;


traj = zeros(3,NPts,NPro);
for i = 1:NPro
    traj(1,:,i) = k_loc*sin(Angs(i,2))*cos(Angs(i,1));
    traj(2,:,i) = k_loc*sin(Angs(i,2))*sin(Angs(i,1));
    traj(3,:,i) = k_loc*cos(Angs(i,2));
end

rad = sqrt(traj(1,:,:).^2 + traj(2,:,:).^2 + traj(3,:,:).^2);
traj = traj/max(rad(:))/2;

% kFOV_desired = 1./(Res);
% kMax_desired = kFOV_desired/2;
% max_k = max(kMax_desired); %Here, we are in 1/m
% 
% traj = traj/max_k/2;

