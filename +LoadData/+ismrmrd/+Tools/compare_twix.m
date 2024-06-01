function is_correct = compare_twix(image_twix,traj_twix)

%Compare the parameters in an imaging file to the parameters in a
%trajectory file to see if the trajectory file has the right values:

%Need all the WIP parameters + matrix size + FOV + nucleus. If all match,
%then we're golden:

%Easy ones first:
%Compare Sequence Name:
Im_Seq = image_twix.hdr.Config.SequenceFileName;
traj_Seq = traj_twix.hdr.Config.SequenceFileName;

if ~strcmp(Im_Seq,traj_Seq)
    is_correct = false;
    return;
end

%Compare Image Size
Im_Size = image_twix.hdr.Dicom.lBaseResolution;
Traj_Size = traj_twix.hdr.Dicom.lBaseResolution;

if Im_Size ~= Traj_Size
    is_correct = false;
    return;
end

%Compare Field of View
Im_FOV = image_twix.hdr.Config.ReadFoV;
Traj_FOV = traj_twix.hdr.Config.ReadFoV;

if Im_FOV ~= Traj_FOV
    is_correct = false;
    return;
end

%Here's where I need to start being careful - Different Versions of the
%Spiral sequences have things in different places on the sequence special
%tab (stupid) (actually, not as bad as I thought I was)

%For the base sequence:
Sp_Loc = 1;
Hubs_Loc = 2;
Alpha_Loc = 3;
UType_Loc = 8;
US0_Loc = 9;
US1_Loc = 10;
USR_Loc = 11;
MeasTraj_Loc = 7;
No_OS_Loc = 6;

%Now, fix these values for the Diffusion Sequence:
if contains(Im_Seq,'DIFF')
    MeasTraj_Loc = 14;
    No_OS_Loc = nan;
end

%Compare Nucleus - edit for XA20A
try
    Im_Nuc = image_twix.hdr.Spice.ResonantNucleus;
catch
    if ~strcmp(image_twix.hdr.Dicom.TransmittingCoil,'129Xe_Vest')
        Im_Nuc = '1H';
    else
        Im_Nuc = '129Xe';
    end
end
Traj_Nuc = traj_twix.hdr.MeasYaps.sWipMemBlock.alFree{MeasTraj_Loc};
if Traj_Nuc == 2
    Traj_Nuc = '129Xe';
elseif Traj_Nuc == 1
    Traj_Nuc = '1H';
end
if ~strcmp(Im_Nuc,Traj_Nuc)
    is_correct = false;
    return;
end

%Compare number of spirals
Im_Sp = image_twix.hdr.MeasYaps.sWipMemBlock.adFree{Sp_Loc};
Traj_Sp = traj_twix.hdr.MeasYaps.sWipMemBlock.adFree{Sp_Loc};

if Im_Sp ~= Traj_Sp
    is_correct = false;
    return;
end

%Compare number of Hubs
Im_Hubs = image_twix.hdr.MeasYaps.sWipMemBlock.alFree{Hubs_Loc};
Traj_Hubs = traj_twix.hdr.MeasYaps.sWipMemBlock.alFree{Hubs_Loc};

if Im_Hubs ~= Traj_Hubs
    is_correct = false;
    return;
end

%Compare Alpha
Im_Alpha = image_twix.hdr.MeasYaps.sWipMemBlock.adFree{Alpha_Loc};
Traj_Alpha = traj_twix.hdr.MeasYaps.sWipMemBlock.adFree{Alpha_Loc};

if Im_Alpha ~= Traj_Alpha
    is_correct = false;
    return;
end

%Compare US_0
Im_US0 = image_twix.hdr.MeasYaps.sWipMemBlock.adFree{US0_Loc};
Traj_US0 = traj_twix.hdr.MeasYaps.sWipMemBlock.adFree{US0_Loc};

if Im_US0 ~= Traj_US0
    is_correct = false;
    return;
end

%Compare US_1
Im_US1 = image_twix.hdr.MeasYaps.sWipMemBlock.adFree{US1_Loc};
Traj_US1 = traj_twix.hdr.MeasYaps.sWipMemBlock.adFree{US1_Loc};

if Im_US1 ~= Traj_US1
    is_correct = false;
    return;
end

%Compare US_R
Im_USR = image_twix.hdr.MeasYaps.sWipMemBlock.adFree{USR_Loc};
Traj_USR = traj_twix.hdr.MeasYaps.sWipMemBlock.adFree{USR_Loc};

if Im_USR ~= Traj_USR
    is_correct = false;
    return;
end

%Compare UType
Im_UType = image_twix.hdr.MeasYaps.sWipMemBlock.alFree{UType_Loc};
Traj_UType = traj_twix.hdr.MeasYaps.sWipMemBlock.alFree{UType_Loc};

if ~isempty(Im_UType)
    if Im_UType ~= Traj_UType
        is_correct = false;
        return;
    end
end

if ~isnan(No_OS_Loc)
    Im_OS = image_twix.hdr.MeasYaps.sWipMemBlock.alFree{No_OS_Loc};
    Traj_OS = traj_twix.hdr.MeasYaps.sWipMemBlock.alFree{No_OS_Loc};
    if Im_OS ~= Traj_OS
        is_correct = false;
        return;
    end
end

%If I've made it to here, everything checks out, and we have the proper
%trajectory.
is_correct = true;
    
    
    
    