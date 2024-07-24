function [Proton,GasExchange] = GasExchange_RegisterProton_to_Xenon(Proton,GasExchange,MainInput)
%   Inputs: 
%      ProtonImage - 3D array:
%      UncorrectedVentImage:
%      VentImage:
%      H_RecMatrix:
%      ProtonSignal:
%      LungSignal:
%      LungStd:
%      NoiseMask:
%      ProtonMax:
%          
%   Outputs:
%      LungMask
%
%   Example: 
%   LoadData_Proton_GasExchange_Philips_Sin('C:\Users\mwillmering\Documents\Subject1')
%
%   Package: https://github.com/cchmc-cpir/CCHMC-Gas-Exchange-Processing-Package
%
%   Author: Matthew Willmering
%   Work email: matthew.willmering@cchmc.org
%   Personal email: matt.willmering@gmail.com
%   Website: https://www.cincinnatichildrens.org/research/divisions/c/cpir

%% Register Proton to Xenon
ProtonImage = Proton.Image;
try
    ProtonImageHR = Proton.ProtonImageHR;
catch
    ProtonImageHR = Proton.Image;
end
UncorrectedVentImage = GasExchange.UncorrectedVentImage;
VentImage = GasExchange.VentImage;
H_RecMatrix = Proton.H_RecMatrix;
ProtonMax = max(abs(ProtonImage(:))); %Proton.ProtonMax;

cd(MainInput.XeDataLocation)
mkdir([MainInput.XeDataLocation '\Gas Exchange Analysis']);
outputpath = [MainInput.XeDataLocation '\Gas Exchange Analysis'];
cd(outputpath) 

NewMask=1;
[LungMask,ProtonSignal,LungSignal,LungStd,NoiseMask] = Registration.GasExchange_MakeProtonMask(ProtonImage,NewMask);
%View
% ProtonMaskMontage = figure('Name','Lung Mask');set(ProtonMaskMontage,'WindowState','minimized');
% montage(LungMask(H_RecMatrix:2*H_RecMatrix,H_RecMatrix:2*H_RecMatrix,H_RecMatrix:2*H_RecMatrix),'DisplayRange',[])%unregistered for these
% disp('Lung Mask Completed.')
% savefig('ProtonMaskMontage.fig')
% close(gcf)


NumPlotSlices = 7;%must be odd to work as expected
%move to imregdemons?
%Set up registration
disp('Registering Proton Image to Gas Image...')
[optimizer,metric] = imregconfig('multimodal');

%Modify images for better registration
VentReg = medfilt3(abs(UncorrectedVentImage), [3 3 3]);
ProtonReg = ProtonImage;
ProtonReg = medfilt3(abs(ProtonReg), [5 5 5]);
ProtonReg(ProtonImage>ProtonSignal) = ProtonMax;%make anything clearly proton the max
ProtonReg(ProtonImage<LungSignal) = LungSignal-(3*LungStd);%make anything clearly lung the min
ProtonReg(NoiseMask) = ProtonMax;%make anything clearly noise the max
ProtonReg = abs(ProtonReg - ProtonMax);%Invert to hopefully help

%create base tform since using 3x overgridded H image
tformBase = affine3d([1,0,0,0;0,1,0,0;0,0,1,0;-H_RecMatrix,-H_RecMatrix,-H_RecMatrix,1]);

%start with similarity transform
optimizer.GrowthFactor = 1.025;
optimizer.Epsilon = 1e-10;
optimizer.InitialRadius = 0.00005;
optimizer.MaximumIterations = 150;
tformSimilarity = imregtform(abs(ProtonReg), abs(VentReg), 'similarity', optimizer, metric,'InitialTransformation',tformBase); %initial similarity transform

%affine transform
optimizer.GrowthFactor = 1.01;
optimizer.Epsilon = 1e-12;
optimizer.InitialRadius = 0.00005;
optimizer.MaximumIterations = 500;
tform = imregtform(abs(ProtonReg), abs(VentReg), 'affine', optimizer, metric,'InitialTransformation',tformSimilarity,'PyramidLevels',4); %affine transform using inital similarity transform
ProtonRegistered = imwarp(abs(ProtonImage), tform,'OutputView',imref3d(size(abs(VentImage))));
ProtonHRRegistered = imwarp(abs(ProtonImageHR), tform,'OutputView',imref3d(size(abs(VentImage))));
ProtonMaskRegistred = logical(imwarp(abs(LungMask), tform,'OutputView',imref3d(size(abs(VentImage)))));

disp('Registering Proton Image to Gas Image Completed.')

%Determine Slices to Plot
%Coronal
Co_MaskVoxels = squeeze(sum(ProtonMaskRegistred,[1,2]));
Co_StartIndex = find(Co_MaskVoxels,1,'first');
Co_EndIndex = find(Co_MaskVoxels,1,'last');
Co_MiddleIndex = round((Co_StartIndex+Co_EndIndex)/2);
Co_Step = floor((Co_EndIndex-Co_StartIndex)/NumPlotSlices);
Slices_Co = (Co_MiddleIndex-Co_Step*(NumPlotSlices-1)/2):Co_Step:(Co_MiddleIndex+Co_Step*(NumPlotSlices-1)/2);
%Axial
Ax_MaskVoxels = squeeze(sum(ProtonMaskRegistred,[2,3]));
Ax_StartIndex = find(Ax_MaskVoxels,1,'first');
Ax_EndIndex = find(Ax_MaskVoxels,1,'last');
Ax_MiddleIndex = round((Ax_StartIndex+Ax_EndIndex)/2);
Ax_Step = floor((Ax_EndIndex-Ax_StartIndex)/NumPlotSlices);
Slices_Ax = (Ax_MiddleIndex-Ax_Step*(NumPlotSlices-1)/2):Ax_Step:(Ax_MiddleIndex+Ax_Step*(NumPlotSlices-1)/2);
%% 

%View
Registrationfig = figure('Name','Registration','units','normalized','outerposition',[0 0 1 4/NumPlotSlices]);set(Registrationfig,'WindowState','minimized');
set(Registrationfig,'color','white','Units','inches','Position',[0.25 0.25 2*NumPlotSlices 4*2+1])
tiledlayout(4,NumPlotSlices,'TileSpacing','none','Padding','compact');
for slice=1:NumPlotSlices
    nexttile
    imshowpair(abs(ProtonImage(H_RecMatrix:2*H_RecMatrix,H_RecMatrix:2*H_RecMatrix,H_RecMatrix+Slices_Co(slice)))/max(abs(ProtonImage(:))),abs(VentImage(:,:,Slices_Co(slice)))/max(abs(VentImage(ProtonMaskRegistred(:)))),'Scaling','none','ColorChannels','red-cyan')
    if slice == round(NumPlotSlices/2)
        title('Before Registration','FontSize',24)
    end
end
for slice=1:NumPlotSlices
    nexttile
    imshowpair(fliplr(rot90(squeeze(abs(ProtonImage(H_RecMatrix+Slices_Ax(slice),H_RecMatrix:2*H_RecMatrix,H_RecMatrix:2*H_RecMatrix))),-1))/max(abs(ProtonImage(:))),fliplr(rot90(squeeze(abs(VentImage(Slices_Ax(slice),:,:))),-1))/max(abs(VentImage(ProtonMaskRegistred(:)))),'Scaling','none','ColorChannels','red-cyan')
end
for slice=1:NumPlotSlices
    nexttile
    imshowpair(abs(ProtonRegistered(:,:,Slices_Co(slice)))/max(abs(ProtonRegistered(:))),abs(VentImage(:,:,Slices_Co(slice)))/max(abs(VentImage(ProtonMaskRegistred(:)))),'Scaling','none','ColorChannels','red-cyan')
    if slice == round(NumPlotSlices/2)
        title('After Registration','FontSize',24)
    end
end
for slice=1:NumPlotSlices
    nexttile
    imshowpair(fliplr(rot90(squeeze(abs(ProtonRegistered(Slices_Ax(slice),:,:))),-1))/max(abs(ProtonRegistered(:))),fliplr(rot90(squeeze(abs(VentImage(Slices_Ax(slice),:,:))),-1))/max(abs(VentImage(ProtonMaskRegistred(:)))),'Scaling','none','ColorChannels','red-cyan')
end
savefig('Registerationfig.fig')
close(gcf)
%% stor data
for slice =1:size(VentImage,3)
    A = ProtonRegistered(:,:,slice);
    B = VentImage(:,:,slice);
    ProtonRegisteredColored(:,:,:,slice) = imfuse(A,B,'falsecolor','ColorChannels','green-magenta'); %'red-cyan'|'green-magenta'
end

Proton.LungMask = LungMask;
Proton.tform = tform;
Proton.Slices_Co = Slices_Co; 
Proton.Slices_Ax = Slices_Ax; 
Proton.optimizer = optimizer;
Proton.ProtonRegistered = ProtonRegistered;
Proton.ProtonHRRegistered = ProtonHRRegistered;
Proton.ProtonMaskRegistred = ProtonMaskRegistred;
Proton.ProtonRegisteredColored = permute(ProtonRegisteredColored,[1 2 4 3]);
GasExchange.LungMask = LungMask;
GasExchange.ProtonRegistered = ProtonRegistered;
GasExchange.ProtonMaskRegistred = ProtonMaskRegistred;
GasExchange.ProtonRegisteredColored = permute(ProtonRegisteredColored,[1 2 4 3]);

%% wite tiff images

Global.write_imshowpair(ProtonRegistered,VentImage,MainInput.XeDataLocation)


end
