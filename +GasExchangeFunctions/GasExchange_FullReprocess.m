function thresholds = GasExchange_FullReprocess(DataFolderLocation,HealthyUpdate)

%Obtain location location to parent folder of all subjects
if exist('DataFolderLocation','var') == 0
    folder = uigetdir('C:\', 'Select Folder containing Gas Exchange Data Folders');
else
    folder = DataFolderLocation;
end
if exist('HealthyUpdate','var') == 0
  HealthyUpdate = 1; %1 process all for a new healthy cohort only
end

if HealthyUpdate == 1 %full reprocess for new healthy cohort
    NewMask = 0; %No need for new mask, keep modifications
    NewImages =  0; %No need for new images
    NewReg = 0; %No need for new registration
else %truly full reprocess
    NewMask = 1; %Note: any external modifications made will be lost!
    NewImages =  1;
    NewReg = 1;
end

%GasExchange_ProcessBatch(DataFolderLocation,NewMask,NewImages,NewReg,HealthyOnly)

%Run healthy and determine initial thresholds
GasExchange_ProcessBatch(folder,NewMask,NewImages,NewReg,1);
thresholds.firstPass = GasExchange_CreateHealthyCohort(folder);

%Check again to make sure thresholds didnt move
GasExchange_ProcessBatch(folder,NewMask,NewImages,NewReg,1);
thresholds.secondPass = GasExchange_CreateHealthyCohort(folder);

%Run all subjects with determined thresholds
GasExchange_ProcessBatch(folder,NewMask,NewImages,NewReg,0);

end