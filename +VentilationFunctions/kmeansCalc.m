
function    [Ventilation] = kmeansCalc(Ventilation,Proton,MainInput)

% You'll need a 2D/3D Xe image (Xenon) and a segmentation (ProtonMask) to run this
XeMasked = Ventilation.Image.*Ventilation.LungMask;
[~,XeMask,~] = VentilationFunctions.CallClustering(XeMasked);
kventilated = XeMask.*double(Ventilation.LungMask);
kdefect = double(Ventilation.LungMask) - kventilated;
kdefectVol = sum(kdefect(:));
kventilatedVol = sum(kventilated(:));
kVDP = kdefectVol / (kdefectVol + kventilatedVol) * 100;
disp(['KmeansVDP= ' num2str(kVDP), '%'])
Ventilation.KmeansVDP = KmeansVDP;
% Cluster Map (2=defect, 1=normal, 0=background)
XeClustersShow = double(Ventilation.LungMask) + kdefect;


% imslice(XeClustersShow)


end