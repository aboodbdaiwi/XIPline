function ...
    [R_map,h_map,r_map,Lm_map,SVR_map,Na_map,DT_map,Dan_map,DL_map,Dmean_map,So_map,MSR_map,Residual_map,LungCMMorphometrySummary]...
    =CM_Lung_Morphometry(...
    Images,...
    final_mask,...
    noise_mask,...
    fitType,...
    bvalues,...
    Do,...
    delta,...
    Datapath,...
    WinBUGSPath)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%% 129Xe Cylinder Model (CM) Lung Morphometry  %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Requirements: WinBUGS program
%   (https://www.mrc-bsu.cam.ac.uk/software/bugs/the-bugs-project-winbugs/)
% 
%   Inputs:
%       Images = 4D images (x,y,slices,bvalue)
%       final_mask = lung mask, 3D (x,y,slices)
%       noise_mask = background mask without airways and artifacts, 3D (x,y,slices)
%       fitType = string: 'human' or 'animals'
%       bvalues = vector: has all b-values (s/cm2)
%       Do = Xenon free diffusion (cm2/s): 0.14
%       delta = The diffusion time (s)
%       Datapath = path to the location of the data
%       WinBUGSPath = location to the winbugs.exe file
%
%   Outputs:
%       R_map = external Acinar duct radius 
%       h_map = alveolar sleeve depth 
%       r_map = internal Acinar duct radius map
%       Lm_map = mean linear intercept map
%       SVR_map = surface to volume ratio map
%       Na_map = alveolar density map
%       DT_map = transverse diffusion coefficient map
%       Dan_map = anisotropic diffusion coefficient map
%       DL_map = longitudinal diffusion coefficient map
%       Dmean_map = mean diffusion coefficient map
%       So_map 
%   Package: 
%
%   Author: Abdullah S. Bdaiwi
%   Work email: abdullah.bdaiwi@cchmc.org
%   Personal email: abdaiwi89@gmail.com
%   Website: https://cpir.cchmc.org/
%
% This code is based on the cylinder model by A. L. Sukstanskii1 and D. A. Yablonskiy
% (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3184317/)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
%normalize images
Images = (Images./max(Images(:))).*100;
Images (Images < 0.0001) = 0.0001; % to avoid having zeros 
Datapath = char(Datapath);
% get rid of any slice that has no signal (un-segmented)
slices=size(Images,3);
slice_checker=zeros(1,slices);
for check_slice=1:slices
    sl=final_mask(:,:,check_slice);
    if sum(sl(:))>1
        slice_checker(check_slice)=1;
    elseif sum(sl(:))==0
        slice_checker(check_slice)=0;
    end
end
find_slice_checker=find(slice_checker);
Images1= zeros(size(Images,1),size(Images,2),length(find_slice_checker), size(Images,4));
Nfinal_mask= zeros(size(Images,1),size(Images,2),length(find_slice_checker));
Nnoise_mask= zeros(size(Images,1),size(Images,2),length(find_slice_checker));
for find_sl=1:length(find_slice_checker)
    choose_sl=find_slice_checker(find_sl);
    Images1(:,:,find_sl,:)=Images(:,:,choose_sl,:);
    Nfinal_mask(:,:,find_sl)=final_mask(:,:,choose_sl);
    Nnoise_mask(:,:,find_sl)=noise_mask(:,:,choose_sl);
end
% figure; imslice(Ndiffimg)
% figure; imslice(Nfinal_mask)
Images=Images1;
final_mask=Nfinal_mask;
noise_mask=Nnoise_mask;

cd(Datapath);
nbvalues=size(Images,4); % number of b values
Runit=1e-4; % R will be calculated in centimeter 

delta = delta/1000; % in us
% calculate L1 and L2 which are the characteristic diffusion lengths for one- and two-dimensional diffusion
L1=sqrt(2*delta*Do); 
L2=sqrt(4*delta*Do);

%error function constant terms
% a1 = 0.0705230784;  a2 = 0.0422820123; a3 = 0.0092705272; a4 = 0.0001520143; a5 = 0.0002765672; a6 = 0.0000430638;
a1 = 0.278393; a2 = 0.230389; a3 = 0.000972; a4 = 0.078108;

Bayes_h=zeros(size(Images,3),size(Images,1)*size(Images,2));
Bayes_So_from_R=zeros(size(Images,3),size(Images,1)*size(Images,2));
Bayes_R=zeros(size(Images,3),size(Images,1)*size(Images,2));
Bayes_DT=zeros(size(Images,3),size(Images,1)*size(Images,2));
Bayes_Dan=zeros(size(Images,3),size(Images,1)*size(Images,2));
Bayes_So_from_D=zeros(size(Images,3),size(Images,1)*size(Images,2));
Bayes_Residual=zeros(size(Images,3),size(Images,1)*size(Images,2),size(Images,4));
Bayes_MSR=zeros(size(Images,3),size(Images,1)*size(Images,2));

for i=1:size(Images,3)
    Image_toFit= squeeze(Images(:,:,i,:)); 
%     Image_toFit=reshape(Image_toFit,[size(Image_toFit,1),size(Image_toFit,2),nbvalues]);
    Image_toFit=(Image_toFit./max(Image_toFit(:))).*100;    

    final_mask_slice=final_mask(:,:,i); 
    noise_mask_slice=noise_mask(:,:,i);
    bvalues(1)=1e-10; % set first b-value close to zero to avoild Inf values when fitting for R,h,DL,DT

    % calculate SNR and sigma
    SNR_vec = zeros(1, nbvalues);
    signal_n_avg = zeros(1, nbvalues);
    std_noise_n = zeros(1, nbvalues); 

    for n = 1:nbvalues            
        signal_vec = Image_toFit(:,:,n).*final_mask_slice; %Calculating the siganl vector
        signal_vec(final_mask_slice==0)=[];
        signal_n_avg(n) = mean(signal_vec);

        noise_vec=Image_toFit(:,:,n).*noise_mask_slice; %Calculating the noise vector
        noise_vec(noise_mask_slice==0)=[];
        std_noise_n (n)= std(noise_vec);            %the standard deviation of the noise
        SNR_vec(n) = signal_n_avg(n) / std_noise_n(n); %signal to noise ratio
        sigma(n) = std_noise_n(1);
    end

    %%%%%%%%%%%%%%%%%%%%%   R and h fit  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    So=Image_toFit(:,:,1).*final_mask_slice; 
    So=reshape(So,[size(Image_toFit,1)*size(Image_toFit,2),1]);
    So_prior = So(any(So,2),:);

    M1 = zeros(size(Image_toFit,1),size(Image_toFit,2),nbvalues);
    for n=1:nbvalues
        M1(:,:,n)=((Image_toFit(:,:,n).*final_mask_slice).^2)/(std_noise_n(1)^2);
    end

    M_Data_vec=reshape(M1,[size(Image_toFit,1)*size(Image_toFit,2),nbvalues]);
    M = M_Data_vec(any(M_Data_vec,2),:);
    num_ones=final_mask_slice;
    num_ones=sum(num_ones(:)); % for loopping    
    
    if strcmp(fitType,'human') == 1
        Fit_Model_R='LM_Model_R_So_Human_1D_h.txt';
        DiffusionFunctions.LM_Model_R_So_Human_1D_h(num_ones,nbvalues); 
        %initializing
        lower_limit=120e-4;
        upper_limit=160e-4;
        h_prior = lower_limit + (upper_limit-lower_limit).*rand(num_ones,1);

        lower_limit=280e-4;
        upper_limit=300e-4;
        R_prior = lower_limit + (upper_limit-lower_limit).*rand(num_ones,1);

    elseif strcmp(fitType,'animals') == 1
        Fit_Model_R='R_So_Model_Animal_1D.txt';
        DiffusionFunctions.LM_Model_R_So_Animal_1D(num_ones,nbvalues);
        %initializing
        lower_limit=30e-4;
        upper_limit=60e-4;
        h_prior = lower_limit + (upper_limit-lower_limit).*rand(num_ones,1);

        lower_limit=80e-4;
        upper_limit=120e-4;
        R_prior = lower_limit + (upper_limit-lower_limit).*rand(num_ones,1);
    end
        
    dataStruct= struct('M', M, 'b', bvalues, 'sigma', sigma, 'Do', Do,...
        'L1', L1, 'L2', L2,'a1',a1,'a2',a2,'a3',a3,'a4',a4);

    initStructs =struct('So',So_prior,'h', h_prior, 'R',R_prior);

    % Call WindBuges
    nchains  = 1; % How Many Chains?
    nburnin  = 1000; % How Many Burn-in Samples? 1000
    nsamples = 500;  % How Many Recorded Samples? 500
    fprintf( 'Running WinBUGS...\n' );

    [samples, stats] = DiffusionFunctions.matbugs(dataStruct, ...
        fullfile(Datapath , Fit_Model_R), ...
        'init', initStructs, ...
        'nChains', nchains, ...
        'view', 0, 'nburnin', nburnin, 'nsamples', nsamples, ...
        'thin', 1, 'DICstatus', 0, 'refreshrate',100,  ...
        'monitorParams', {'So','h','R','MSR','residual'}, ...
        'Bugdir', WinBUGSPath);

    Bayes_So_from_R(i,1:num_ones) =stats.mean.So;
    Bayes_h(i,1:num_ones) =stats.mean.h./Runit;
    Bayes_R(i,1:num_ones) =stats.mean.R./Runit;
    Bayes_Residual(i,1:num_ones,:)=stats.mean.residual;
    Bayes_MSR(i,1:num_ones)=stats.mean.MSR;

    %%%%%%%%%%%%%%%%%%%%% DT and DL fitting  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmp(fitType,'human') == 1
        Fit_Model_D='LM_Model_D_Human_1D.txt';
        DiffusionFunctions.LM_Model_D_Human_1D(num_ones,nbvalues);

        %initializing
        lower_limit=0.001;
        upper_limit=0.1;
        DT_prior = lower_limit + (upper_limit-lower_limit).*rand(num_ones,1);

        lower_limit=0.001;
        upper_limit=0.1;
        Dan_prior = lower_limit + (upper_limit-lower_limit).*rand(num_ones,1);

    elseif strcmp(fitType,'animals') == 1
        Fit_Model_D='LM_Model_D_Animal_1D.txt';
        DiffusionFunctions.LM_Model_D_Animal_1D(num_ones,length(bvalues));
        
        %initializing
        lower_limit=0.01;
        upper_limit=0.04;
        DT_prior = lower_limit + (upper_limit-lower_limit).*rand(num_ones,1);

        lower_limit=0.01;
        upper_limit=0.04;
        Dan_prior = lower_limit + (upper_limit-lower_limit).*rand(num_ones,1);
    end
        bvalues(1)=0;
    dataStruct= struct('M', M, 'b', bvalues, 'sigma', sigma,...
                'a1',a1,'a2',a2,'a3',a3,'a4',a4);

    initStructs =struct('So',So_prior,'DT', DT_prior, 'Dan',Dan_prior);

    % Call WindBuges
    nchains  = 1; % How Many Chains?
    nburnin  = 2000; % How Many Burn-in Samples?
    nsamples = 2000;  % How Many Recorded Samples?
    fprintf( 'Running WinBUGS...\n' );

    [samples, stats] = DiffusionFunctions.matbugs(dataStruct, ...
        fullfile(Datapath , Fit_Model_D), ...
        'init', initStructs, ...
        'nChains', nchains, ...
        'view', 0, 'nburnin', nburnin, 'nsamples', nsamples, ...
        'thin', 1, 'DICstatus', 0, 'refreshrate',100,  ...
        'monitorParams', {'So','DT','Dan'}, ...
        'Bugdir', WinBUGSPath);

    Bayes_DT(i,1:num_ones) =stats.mean.DT;
    Bayes_Dan(i,1:num_ones) =stats.mean.Dan;
    Bayes_So_from_D(i,1:num_ones)=stats.mean.So;
    
end
 
%%%%%%%%%%%%%%%%%%%%convert vector to 2D image  %%%%%%%%%%%%%%%%%%%%%%%%%%% 
pixel_location=find(final_mask); %to find pixel location 
num_ones2=sum(final_mask(:)); 

Bayes_h_2=Bayes_h';
Bayes_h_2=reshape(Bayes_h_2,[size(final_mask,1)*size(final_mask,2)*size(final_mask,3),1]);
Bayes_h_2 = Bayes_h_2(any(Bayes_h_2,2),:);

Bayes_So_from_R2=Bayes_So_from_R';
Bayes_So_from_R2=reshape(Bayes_So_from_R2,[size(final_mask,1)*size(final_mask,2)*size(final_mask,3),1]);
Bayes_So_from_R2 = Bayes_So_from_R2(any(Bayes_So_from_R2,2),:);

Bayes_R_2=Bayes_R';
Bayes_R_2=reshape(Bayes_R_2,[size(final_mask,1)*size(final_mask,2)*size(final_mask,3),1]);
Bayes_R_2 = Bayes_R_2(any(Bayes_R_2,2),:);

for ii=1:size(Images,4)
    b_Residual2 =Bayes_Residual(:,:,ii);
    b_Residual2=b_Residual2';
    Bayes_Residual2=reshape(b_Residual2,[size(final_mask,1)*size(final_mask,2)*size(Images,3),1]);
    Bayes_Residual2 = Bayes_Residual2(any(Bayes_Residual2,2),:);
    Bayes_Residual3(:,ii)=Bayes_Residual2;
end

Bayes_MSR=Bayes_MSR';
Bayes_MSR=reshape(Bayes_MSR,[size(final_mask,1)*size(final_mask,2)*size(final_mask,3),1]);
Bayes_MSR = Bayes_MSR(any(Bayes_MSR,2),:);

Bayes_DT_2=Bayes_DT';
Bayes_DT_2=reshape(Bayes_DT_2,[size(final_mask,1)*size(final_mask,2)*size(final_mask,3),1]);
Bayes_DT_2 = Bayes_DT_2(any(Bayes_DT_2,2),:);

Bayes_Dan_2=Bayes_Dan';
Bayes_Dan_2=reshape(Bayes_Dan_2,[size(final_mask,1)*size(final_mask,2)*size(final_mask,3),1]);
Bayes_Dan_2 = Bayes_Dan_2(any(Bayes_Dan_2,2),:);

Bayes_So_from_D2=Bayes_So_from_D';
Bayes_So_from_D2=reshape(Bayes_So_from_D2,[size(final_mask,1)*size(final_mask,2)*size(final_mask,3),1]);
Bayes_So_from_D2 = Bayes_So_from_D2(any(Bayes_So_from_D2,2),:);

h_map=zeros(size(final_mask));
Bayes_So_from_R_map=zeros(size(final_mask));
R_map=zeros(size(final_mask));
DT_map=zeros(size(final_mask));
Dan_map=zeros(size(final_mask));
Bayes_So_from_D_map=zeros(size(final_mask)); 
Bayes_Residual_map=zeros(size(Images)); 
Bayes_MSR_map=zeros(size(final_mask)); 
b_Residual=zeros(size(final_mask)); 

for bayes_loop=1:num_ones2
    choosen_pixel=pixel_location(bayes_loop);
    h_map(choosen_pixel) = Bayes_h_2(bayes_loop);
    Bayes_So_from_R_map(choosen_pixel) = Bayes_So_from_R2(bayes_loop);
    R_map(choosen_pixel) = Bayes_R_2(bayes_loop);
    DT_map(choosen_pixel) = Bayes_DT_2(bayes_loop);
    Dan_map(choosen_pixel) = Bayes_Dan_2(bayes_loop);
    Bayes_So_from_D_map(choosen_pixel) = Bayes_So_from_D2(bayes_loop);    
    Bayes_MSR_map(choosen_pixel) = Bayes_MSR(bayes_loop);

end 

for jj=1:size(Images,4)
    for bayes_loop=1:num_ones2
        choosen_pixel=pixel_location(bayes_loop);
        chosen_b=Bayes_Residual3(:,jj);
        b_Residual(choosen_pixel) = chosen_b(bayes_loop); 
        Bayes_Residual_map(:,:,:,jj) = b_Residual; 
    end
end    


So_map=Bayes_So_from_R_map;
% So_map=Bayes_So_from_D_map;

%%%%%%%%%%%%%%%%%% DL, r, L, SVR, Lm and Na  %%%%%%
DL_map=DT_map+Dan_map;
Dmean_map = DL_map-((2/3).*Dan_map);
r_map=R_map-h_map;
L_map = 0.765.*R_map; 

SVR_map=((2*pi.*R_map.*L_map)+2*pi.*((R_map.^2)-(r_map.^2))+16.*h_map.*L_map)./(pi.*(R_map.^2).*L_map);
Lm_map =4./SVR_map;
Na_map = 1./(pi.*(R_map.^2).*L_map);

SVR_map=SVR_map./Runit;% (in cm-1) 
Na_map=Na_map.*10.0e+9; % in mm^-3

SVR_map(isinf(SVR_map)|isnan(SVR_map))=0;
Na_map(isinf(Na_map)|isnan(Na_map))=0;

MSR_map=Bayes_MSR_map;
Residual_map = Bayes_Residual_map;    

%% Find mean and std of all maps
DT_mean = round(mean(DT_map(DT_map>0)),3); 
DL_mean = round(mean(DL_map(DL_map>0)),3); 
Dan_mean = round(mean(Dan_map(Dan_map>0)),3); 
Dmean_mean = round(mean(Dmean_map(Dmean_map>0)),3); 
R_mean = round(mean(R_map(R_map>0)));
h_mean = round(mean(h_map(h_map>0))); 
r_mean = round(mean(r_map(r_map>0))); 
Lm_mean = round(mean(Lm_map(Lm_map>0))); 
SVR_mean = round(mean(SVR_map(SVR_map>0)));
Na_mean = round(mean(Na_map(Na_map>0))); 

DT_std = round(std(DT_map(DT_map>0)),3); 
DL_std = round(std(DL_map(DL_map>0)),3);
Dan_std = round(std(Dan_map(Dan_map>0)),3); 
Dmean_std = round(std(Dmean_map(Dmean_map>0)),3); 
R_std = round(std(R_map(R_map>0))); 
h_std = round(std(h_map(h_map>0)));
r_std = round(std(r_map(r_map>0))); 
Lm_std = round(std(Lm_map(Lm_map>0))); 
SVR_std = round(std(SVR_map(SVR_map>0))); 
Na_std = round(std(Na_map(Na_map>0))); 

LungCMMorphometrySummary = table({'Mean';'SD'},[DT_mean;DT_std],[DL_mean;DL_std],[Dan_mean;Dan_std],[Dan_mean;Dmean_std],...
    [R_mean;R_std],[h_mean;h_std],[r_mean;r_std],[Lm_mean;Lm_std],...
    [SVR_mean;SVR_std],[Na_mean;Na_std]);
headers=({'summary','DT(cm^2/s)', 'DL(cm^2/s)', 'Dan(cm^2/s)', 'Dmean(cm^2/s)',...
    'R(um)', 'h(um)', 'r(um)', 'Lm(um)',...
    'SVR(cm^-1)', 'Na(mm^-3)'});
LungCMMorphometrySummary.Properties.VariableNames = headers;
writetable(LungCMMorphometrySummary,[Datapath,'CM_MorphometrySummary.xlsx'],'Sheet',1)
%% write tiff maps
cmap = jet;
cmap(1,:) = 0;
close all;

tiff = figure('MenuBar','none','ToolBar','none','DockControls','off','Resize','off','WindowState','minimized');%figure for tiffs
ax1 = axes('Parent',tiff);ax2 = axes('Parent',tiff);%make axis for both images
set(ax1,'Visible','off');set(ax2,'Visible','off');%turn off axis 
set(ax1,'units','inches');set(ax2,'units','inches');%make axis units inches
set(ax1,'position',[0 0 2 2]);set(ax2,'position',[0 0 2 2]);%make axis same as image
set(gcf,'units','inches'); % set the figure units to pixels
set(gcf,'position',[1 1 2 2])% set the position of the figure to axes
disp('Saving Morphometrymays into Tiff...')
Proton_Image = zeros(size(DT_map));

%DT_map 
for slice=1:size(DT_map,3) %repeat for rest of slices
    [~,~] = Global.imoverlay(squeeze(abs(Proton_Image(:,:,slice))),squeeze(DT_map(:,:,slice)),[],[0,0.99*max(Proton_Image(:))],cmap,1,gca);
    colormap(gca,cmap); caxis([0 0.14]);   
%     X = print('-RGBImage',['-r',num2str(size(DT_map,2)/2)]);%2 inches
     Xdata = getframe(gcf);
     X = Xdata.cdata;     
    if (slice == 1)
        imwrite(X,[Datapath,'\DT_map.tif'],'Description',strcat('Package Version: ', '1','; Cohort: ', 'test'));%write new/ overwrite tiff
    else
        imwrite(X,[Datapath,'\DT_map.tif'],'WriteMode','append','Description',strcat('Package Version: ', '1','; Cohort: ', 'test'));%append tiff
    end
end

%DL_map 
for slice=1:size(DL_map,3) %repeat for rest of slices
    [~,~] = Global.imoverlay(squeeze(abs(Proton_Image(:,:,slice))),squeeze(DL_map(:,:,slice)),[],[0,0.99*max(Proton_Image(:))],cmap,1,gca);
    colormap(gca,cmap); caxis([0 0.14]);   
%     X = print('-RGBImage',['-r',num2str(size(DL_map,2)/2)]);%2 inches
     Xdata = getframe(gcf);
     X = Xdata.cdata;    
    if (slice == 1)
        imwrite(X,[Datapath,'\DL_map.tif'],'Description',strcat('Package Version: ', '1','; Cohort: ', 'test'));%write new/ overwrite tiff
    else
        imwrite(X,[Datapath,'\DL_map.tif'],'WriteMode','append','Description',strcat('Package Version: ', '1','; Cohort: ', 'test'));%append tiff
    end
end

%Dan_map 
for slice=1:size(Dan_map,3) %repeat for rest of slices
    [~,~] = Global.imoverlay(squeeze(abs(Proton_Image(:,:,slice))),squeeze(Dan_map(:,:,slice)),[],[0,0.99*max(Proton_Image(:))],cmap,1,gca);
    colormap(gca,cmap); caxis([0 0.14]);   
%     X = print('-RGBImage',['-r',num2str(size(Dan_map,2)/2)]);%2 inches
     Xdata = getframe(gcf);
     X = Xdata.cdata;     
    if (slice == 1)
        imwrite(X,[Datapath,'\Dan_map.tif'],'Description',strcat('Package Version: ', '1','; Cohort: ', 'test'));%write new/ overwrite tiff
    else
        imwrite(X,[Datapath,'\Dan_map.tif'],'WriteMode','append','Description',strcat('Package Version: ', '1','; Cohort: ', 'test'));%append tiff
    end
end

%Dmean_map 
for slice=1:size(Dmean_map,3) %repeat for rest of slices
    [~,~] = Global.imoverlay(squeeze(abs(Proton_Image(:,:,slice))),squeeze(Dmean_map(:,:,slice)),[],[0,0.99*max(Proton_Image(:))],cmap,1,gca);
    colormap(gca,cmap); caxis([0 0.14]);   
%     X = print('-RGBImage',['-r',num2str(size(Dmean_map,2)/2)]);%2 inches
     Xdata = getframe(gcf);
     X = Xdata.cdata;    
    if (slice == 1)
        imwrite(X,[Datapath,'\Dmean_map.tif'],'Description',strcat('Package Version: ', '1','; Cohort: ', 'test'));%write new/ overwrite tiff
    else
        imwrite(X,[Datapath,'\Dmean_map.tif'],'WriteMode','append','Description',strcat('Package Version: ', '1','; Cohort: ', 'test'));%append tiff
    end
end

%R_map 
for slice=1:size(R_map,3) %repeat for rest of slices
    [~,~] = Global.imoverlay(squeeze(abs(Proton_Image(:,:,slice))),squeeze(R_map(:,:,slice)),[],[0,0.99*max(Proton_Image(:))],cmap,1,gca);
    colormap(gca,cmap); caxis([0 400]);   
%     X = print('-RGBImage',['-r',num2str(size(R_map,2)/2)]);%2 inches
     Xdata = getframe(gcf);
     X = Xdata.cdata;     
    if (slice == 1)
        imwrite(X,[Datapath,'\R_map.tif'],'Description',strcat('Package Version: ', '1','; Cohort: ', 'test'));%write new/ overwrite tiff
    else
        imwrite(X,[Datapath,'\R_map.tif'],'WriteMode','append','Description',strcat('Package Version: ', '1','; Cohort: ', 'test'));%append tiff
    end
end

%h_map 
for slice=1:size(h_map,3) %repeat for rest of slices
    [~,~] = Global.imoverlay(squeeze(abs(Proton_Image(:,:,slice))),squeeze(h_map(:,:,slice)),[],[0,0.99*max(Proton_Image(:))],cmap,1,gca);
    colormap(gca,cmap); caxis([0 400]);   
%     X = print('-RGBImage',['-r',num2str(size(h_map,2)/2)]);%2 inches
     Xdata = getframe(gcf);
     X = Xdata.cdata;    
    if (slice == 1)
        imwrite(X,[Datapath,'\h_map.tif'],'Description',strcat('Package Version: ', '1','; Cohort: ', 'test'));%write new/ overwrite tiff
    else
        imwrite(X,[Datapath,'\h_map.tif'],'WriteMode','append','Description',strcat('Package Version: ', '1','; Cohort: ', 'test'));%append tiff
    end
end

%r_map 
for slice=1:size(r_map,3) %repeat for rest of slices
    [~,~] = Global.imoverlay(squeeze(abs(Proton_Image(:,:,slice))),squeeze(r_map(:,:,slice)),[],[0,0.99*max(Proton_Image(:))],cmap,1,gca);
    colormap(gca,cmap); caxis([0 400]);   
%     X = print('-RGBImage',['-r',num2str(size(r_map,2)/2)]);%2 inches
     Xdata = getframe(gcf);
     X = Xdata.cdata;     
    if (slice == 1)
        imwrite(X,[Datapath,'\r_mapp.tif'],'Description',strcat('Package Version: ', '1','; Cohort: ', 'test'));%write new/ overwrite tiff
    else
        imwrite(X,[Datapath,'\r_mapp.tif'],'WriteMode','append','Description',strcat('Package Version: ', '1','; Cohort: ', 'test'));%append tiff
    end
end

%Lm_map 
for slice=1:size(Lm_map,3) %repeat for rest of slices
    [~,~] = Global.imoverlay(squeeze(abs(Proton_Image(:,:,slice))),squeeze(Lm_map(:,:,slice)),[],[0,0.99*max(Proton_Image(:))],cmap,1,gca);
    colormap(gca,cmap); caxis([0 400]);   
%     X = print('-RGBImage',['-r',num2str(size(Lm_map,2)/2)]);%2 inches
     Xdata = getframe(gcf);
     X = Xdata.cdata;     
    if (slice == 1)
        imwrite(X,[Datapath,'\Lm_map.tif'],'Description',strcat('Package Version: ', '1','; Cohort: ', 'test'));%write new/ overwrite tiff
    else
        imwrite(X,[Datapath,'\Lm_map.tif'],'WriteMode','append','Description',strcat('Package Version: ', '1','; Cohort: ', 'test'));%append tiff
    end
end

%SVR_map 
for slice=1:size(SVR_map,3) %repeat for rest of slices
    [~,~] = Global.imoverlay(squeeze(abs(Proton_Image(:,:,slice))),squeeze(SVR_map(:,:,slice)),[],[0,0.99*max(Proton_Image(:))],cmap,1,gca);
    colormap(gca,cmap); caxis([0 400]);   
%     X = print('-RGBImage',['-r',num2str(size(SVR_map,2)/2)]);%2 inches
     Xdata = getframe(gcf);
     X = Xdata.cdata;   
    if (slice == 1)
        imwrite(X,[Datapath,'\SVR_map.tif'],'Description',strcat('Package Version: ', '1','; Cohort: ', 'test'));%write new/ overwrite tiff
    else
        imwrite(X,[Datapath,'\SVR_map.tif'],'WriteMode','append','Description',strcat('Package Version: ', '1','; Cohort: ', 'test'));%append tiff
    end
end

%Na_map 
for slice=1:size(Na_map,3) %repeat for rest of slices
    [~,~] = Global.imoverlay(squeeze(abs(Proton_Image(:,:,slice))),squeeze(Na_map(:,:,slice)),[],[0,0.99*max(Proton_Image(:))],cmap,1,gca);
    colormap(gca,cmap); caxis([0 400]);   
%     X = print('-RGBImage',['-r',num2str(size(Na_map,2)/2)]);%2 inches
     Xdata = getframe(gcf);
     X = Xdata.cdata;     
    if (slice == 1)
        imwrite(X,[Datapath,'\Na_map.tif'],'Description',strcat('Package Version: ', '1','; Cohort: ', 'test'));%write new/ overwrite tiff
    else
        imwrite(X,[Datapath,'\Na_map.tif'],'WriteMode','append','Description',strcat('Package Version: ', '1','; Cohort: ', 'test'));%append tiff
    end
end

disp('Saving DT_map Tiff Completed.')
close all;

%% plot and save results
cmap = colormap (jet);
cmap(1,:)=0;
close all

%DT maps
DT_Montage = figure('Name','DT Images');set(DT_Montage,'WindowState','minimized');
montage(reshape(DT_map,[size(DT_map,1), size(DT_map,2), 1, size(DT_map,3)]),...
    'Size',[1 size(DT_map,3)],'DisplayRange',[0 1]);
set(gca,'units','pixels'); % set the axes units to pixels
x = get(gca,'position'); % get the position of the axes
set(gcf,'units','pixels'); % set the figure units to pixels
y = get(gcf,'position'); % get the figure position
set(gcf,'position',[y(1) y(2) x(3) x(4)])% set the position of the figure to the length and width of the axes
set(gca,'units','normalized','position',[0 0 1 1]) % set the axes units to pixels
if strcmp(fitType,'human') == 1
    colormap(cmap); caxis([0 0.14]); 
elseif strcmp(fitType,'animals') == 1
    colormap(cmap); caxis([0 0.06]); 
end
xx = colorbar;
xx.Label.String = ' (cm^2.S^-^1)';
DT_MontagePosition=get(gcf,'position');

%DL maps
DL_Montage = figure('Name','DL Images');set(DL_Montage,'WindowState','minimized');
montage(reshape(DL_map,[size(DL_map,1), size(DL_map,2), 1, size(DL_map,3)]),...
    'Size',[1 size(DL_map,3)],'DisplayRange',[0 1]);
set(gca,'units','pixels'); % set the axes units to pixels
x = get(gca,'position'); % get the position of the axes
set(gcf,'units','pixels'); % set the figure units to pixels
y = get(gcf,'position'); % get the figure position
set(gcf,'position',[y(1) y(2) x(3) x(4)])% set the position of the figure to the length and width of the axes
set(gca,'units','normalized','position',[0 0 1 1]) % set the axes units to pixels
if strcmp(fitType,'human') == 1
    colormap(cmap); caxis([0 0.14]); 
elseif strcmp(fitType,'animals') == 1
    colormap(cmap); caxis([0 0.06]); 
end
xx = colorbar;
xx.Label.String = ' (cm^2.S^-^1)';
DL_MontagePosition=get(gcf,'position');

%Dan maps
Dan_Montage = figure('Name','Dan Images');set(Dan_Montage,'WindowState','minimized');
montage(reshape(Dan_map,[size(Dan_map,1), size(Dan_map,2), 1, size(Dan_map,3)]),...
    'Size',[1 size(Dan_map,3)],'DisplayRange',[0 1]);
set(gca,'units','pixels'); % set the axes units to pixels
x = get(gca,'position'); % get the position of the axes
set(gcf,'units','pixels'); % set the figure units to pixels
y = get(gcf,'position'); % get the figure position
set(gcf,'position',[y(1) y(2) x(3) x(4)])% set the position of the figure to the length and width of the axes
set(gca,'units','normalized','position',[0 0 1 1]) % set the axes units to pixels
if strcmp(fitType,'human') == 1
    colormap(cmap); caxis([0 0.14]); 
elseif strcmp(fitType,'animals') == 1
    colormap(cmap); caxis([0 0.06]); 
end
xx = colorbar;
xx.Label.String = ' (cm^2.S^-^1)';
Dan_MontagePosition=get(gcf,'position');

%Dmean maps
Dmean_Montage = figure('Name','Dmean Images');set(Dmean_Montage,'WindowState','minimized');
montage(reshape(Dmean_map,[size(Dmean_map,1), size(Dmean_map,2), 1, size(Dmean_map,3)]),...
    'Size',[1 size(Dmean_map,3)],'DisplayRange',[0 1]);
set(gca,'units','pixels'); % set the axes units to pixels
x = get(gca,'position'); % get the position of the axes
set(gcf,'units','pixels'); % set the figure units to pixels
y = get(gcf,'position'); % get the figure position
set(gcf,'position',[y(1) y(2) x(3) x(4)])% set the position of the figure to the length and width of the axes
set(gca,'units','normalized','position',[0 0 1 1]) % set the axes units to pixels
if strcmp(fitType,'human') == 1
    colormap(cmap); caxis([0 0.14]); 
elseif strcmp(fitType,'animals') == 1
    colormap(cmap); caxis([0 0.06]); 
end
xx = colorbar;
xx.Label.String = ' (cm^2.S^-^1)';
Dmean_MontagePosition=get(gcf,'position');

%R maps
R_Montage = figure('Name','R Images');set(R_Montage,'WindowState','minimized');
montage(reshape(R_map,[size(R_map,1), size(R_map,2), 1, size(R_map,3)]),...
    'Size',[1 size(R_map,3)],'DisplayRange',[0 1]);
set(gca,'units','pixels'); % set the axes units to pixels
x = get(gca,'position'); % get the position of the axes
set(gcf,'units','pixels'); % set the figure units to pixels
y = get(gcf,'position'); % get the figure position
set(gcf,'position',[y(1) y(2) x(3) x(4)])% set the position of the figure to the length and width of the axes
set(gca,'units','normalized','position',[0 0 1 1]) % set the axes units to pixels
if strcmp(fitType,'human') == 1
    colormap(cmap); caxis([0 400]); colorbar;
elseif strcmp(fitType,'animals') == 1
    colormap(cmap); caxis([0 150]); colorbar;
end
xx = colorbar;
xx.Label.String = ' (\mum)';
R_MontagePosition=get(gcf,'position');

%h maps
h_Montage = figure('Name','h Images');set(h_Montage,'WindowState','minimized');
montage(reshape(h_map,[size(h_map,1), size(h_map,2), 1, size(h_map,3)]),...
    'Size',[1 size(h_map,3)],'DisplayRange',[0 1]);
set(gca,'units','pixels'); % set the axes units to pixels
x = get(gca,'position'); % get the position of the axes
set(gcf,'units','pixels'); % set the figure units to pixels
y = get(gcf,'position'); % get the figure position
set(gcf,'position',[y(1) y(2) x(3) x(4)])% set the position of the figure to the length and width of the axes
set(gca,'units','normalized','position',[0 0 1 1]) % set the axes units to pixels
if strcmp(fitType,'human') == 1
    colormap(cmap); caxis([0 400]); colorbar;
elseif strcmp(fitType,'animals') == 1
    colormap(cmap); caxis([0 150]); colorbar;
end
xx = colorbar;
xx.Label.String = ' (\mum)';
h_MontagePosition=get(gcf,'position');

%r maps
r_Montage = figure('Name','r Images');set(r_Montage,'WindowState','minimized');
montage(reshape(r_map,[size(r_map,1), size(r_map,2), 1, size(r_map,3)]),...
    'Size',[1 size(r_map,3)],'DisplayRange',[0 1]);
set(gca,'units','pixels'); % set the axes units to pixels
x = get(gca,'position'); % get the position of the axes
set(gcf,'units','pixels'); % set the figure units to pixels
y = get(gcf,'position'); % get the figure position
set(gcf,'position',[y(1) y(2) x(3) x(4)])% set the position of the figure to the length and width of the axes
set(gca,'units','normalized','position',[0 0 1 1]) % set the axes units to pixels
if strcmp(fitType,'human') == 1
    colormap(cmap); caxis([0 400]); colorbar;
elseif strcmp(fitType,'animals') == 1
    colormap(cmap); caxis([0 150]); colorbar;
end
xx = colorbar;
xx.Label.String = ' (\mum)';
r_MontagePosition=get(gcf,'position');

%Lm maps
Lm_map(isnan(Lm_map))=0;
Lm_Montage = figure('Name','Lm Images');set(Lm_Montage,'WindowState','minimized');
montage(reshape(Lm_map,[size(Lm_map,1), size(Lm_map,2), 1, size(Lm_map,3)]),...
    'Size',[1 size(Lm_map,3)],'DisplayRange',[0 1]);
set(gca,'units','pixels'); % set the axes units to pixels
x = get(gca,'position'); % get the position of the axes
set(gcf,'units','pixels'); % set the figure units to pixels
y = get(gcf,'position'); % get the figure position
set(gcf,'position',[y(1) y(2) x(3) x(4)])% set the position of the figure to the length and width of the axes
set(gca,'units','normalized','position',[0 0 1 1]) % set the axes units to pixels
if strcmp(fitType,'human') == 1
    colormap(cmap); caxis([0 400]); colorbar;
elseif strcmp(fitType,'animals') == 1
    colormap(cmap); caxis([0 150]); colorbar;
end
xx = colorbar;
xx.Label.String = ' (\mum)';
Lm_MontagePosition=get(gcf,'position');

%SVR maps
SVR_map(isinf(SVR_map)|isnan(SVR_map))=0;
SVR_Montage = figure('Name','SVR Images');set(SVR_Montage,'WindowState','minimized');
montage(reshape(SVR_map,[size(SVR_map,1), size(SVR_map,2), 1, size(SVR_map,3)]),...
    'Size',[1 size(SVR_map,3)],'DisplayRange',[]);
set(gca,'units','pixels'); % set the axes units to pixels
x = get(gca,'position'); % get the position of the axes
set(gcf,'units','pixels'); % set the figure units to pixels
y = get(gcf,'position'); % get the figure position
set(gcf,'position',[y(1) y(2) x(3) x(4)])% set the position of the figure to the length and width of the axes
set(gca,'units','normalized','position',[0 0 1 1]) % set the axes units to pixels
if strcmp(fitType,'human') == 1
    colormap(cmap); caxis([0 400]); colorbar;
elseif strcmp(fitType,'animals') == 1
    colormap(cmap); caxis([0 2000]); colorbar;
end
xx = colorbar;
xx.Label.String = ' (cm^-^1)';
SVR_MontagePosition=get(gcf,'position');

%Na maps
Na_map(isinf(Na_map)|isnan(Na_map))=0;
Na_Montage = figure('Name','Na Images');set(Na_Montage,'WindowState','minimized');
montage(reshape(Na_map,[size(Na_map,1), size(Na_map,2), 1, size(Na_map,3)]),...
    'Size',[1 size(Na_map,3)],'DisplayRange',[]);
set(gca,'units','pixels'); % set the axes units to pixels
x = get(gca,'position'); % get the position of the axes
set(gcf,'units','pixels'); % set the figure units to pixels
y = get(gcf,'position'); % get the figure position
set(gcf,'position',[y(1) y(2) x(3) x(4)])% set the position of the figure to the length and width of the axes
set(gca,'units','normalized','position',[0 0 1 1]) % set the axes units to pixels
if strcmp(fitType,'human') == 1
    colormap(cmap); caxis([0 400]); colorbar;
elseif strcmp(fitType,'animals') == 1
    colormap(cmap); caxis([0 20000]); colorbar;
end
xx = colorbar;
xx.Label.String = ' (mm^-^3)';
Na_MontagePosition=get(gcf,'position');

% save maps in ppt file
nslices = size(Na_map,3);
if nslices > 8
    nslices = 8;
end
scale_fac = 2*nslices;
% nslices=size(Na_map,3);

ReportTitle = 'Diffusion_Analysis';
cd(Datapath);
%Start new presentation
isOpen  = Global.exportToPPTX();
if ~isempty(isOpen) %If PowerPoint already started, then close first and then open a new one
    Global.exportToPPTX('close');
end
ppt_file_name = 'Diffusion_Analysis.pptx';
if exist(ppt_file_name, 'file') == 2            
    Global.exportToPPTX('open',ppt_file_name);
    Global.exportToPPTX('switchslide',1);
else            
    Global.exportToPPTX('new','Dimensions',[16 9], ...
        'Title',ReportTitle, ...
        'Author','CPIR @ CCHMC');
end 
%Add slides
%%%%%%%%%%%%%%%%%%%%%%%%%% DT/DL/Dan/Dmean %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Global.exportToPPTX('addslide'); % DT/DL/Dan/Dmean
Global.exportToPPTX('addtext',sprintf('CM-Morphometry Analysis'),'Position',[6 0 5 1],'Color','b','FontSize',25);
Global.exportToPPTX('addpicture',DT_Montage,'Position',[0 .9 scale_fac scale_fac*(DT_MontagePosition(4)/DT_MontagePosition(3))]);
Global.exportToPPTX('addpicture',DL_Montage,'Position',[0 2.9 scale_fac scale_fac*(DL_MontagePosition(4)/DL_MontagePosition(3))]);
Global.exportToPPTX('addpicture',Dan_Montage,'Position',[0 4.9 scale_fac scale_fac*(Dan_MontagePosition(4)/Dan_MontagePosition(3))]);
Global.exportToPPTX('addpicture',Dmean_Montage,'Position',[0 6.9 scale_fac scale_fac*(Dmean_MontagePosition(4)/Dmean_MontagePosition(3))]);

DT_title = sprintf('DT-Mean= %1.3f ± %1.3f (cm^2/s)',DT_mean,DT_std);
Global.exportToPPTX('addtext',DT_title,'Position',[0 0.5 6 1],'Color','r','FontSize',20);
DL_title = sprintf('DL-Mean= %1.3f ± %1.3f (cm^2/s)',DL_mean,DL_std);    
Global.exportToPPTX('addtext',DL_title,'Position',[0 2.5 6 1],'Color','r','FontSize',20);
Dan_title = sprintf('Dan-Mean= %1.3f ± %1.3f (cm^2/s)',Dan_mean,Dan_std);    
Global.exportToPPTX('addtext',Dan_title,'Position',[0 4.5 6 1],'Color','r','FontSize',20);
Dmean_title = sprintf('Dmean-Mean= %1.3f ± %1.3f (cm^2/s)',Dmean_mean,Dmean_std);    
Global.exportToPPTX('addtext',Dmean_title,'Position',[0 6.5 6 1],'Color','r','FontSize',20);

Global.exportToPPTX('addslide'); %DT/DL/Dan/Dmean Histogram
DT_histFig = figure; histogram(DT_map(DT_map>0),100); 
set(DT_histFig,'Name','DT_histFig'); title(DT_histFig.CurrentAxes,DT_title);ylabel(DT_histFig.CurrentAxes,'Pixel Cont.');xlabel(DT_histFig.CurrentAxes,'DT (cm^2/s)');
set(gca,'FontSize',13,'FontWeight','bold'); xlim([0 0.14]);
% print(DT_histFig,[Datapath , '\','DT_histFig'],'-dpng','-r300');
saveas(gca,'DT_histFig.png');
Global.exportToPPTX('addpicture',DT_histFig,'Position',[0 0 7 4.5]);

DL_histFig = figure; histogram(DL_map(DL_map>0),100); 
set(DL_histFig,'Name','DL_histFig'); title(DL_histFig.CurrentAxes,DL_title);ylabel(DL_histFig.CurrentAxes,'Pixel Cont.');xlabel(DL_histFig.CurrentAxes,'DL (cm^2/s)');
set(gca,'FontSize',13,'FontWeight','bold'); xlim([0 0.14]);
% print(DL_histFig,[Datapath , '\','DL_histFig'],'-dpng','-r300');
saveas(gca,'DL_histFig.png');
Global.exportToPPTX('addpicture',DL_histFig,'Position',[7.5 0 7 4.5]);

Dan_histFig = figure; histogram(Dan_map(Dan_map>0),100); 
set(Dan_histFig,'Name','Dan_histFig'); title(Dan_histFig.CurrentAxes,Dan_title);ylabel(Dan_histFig.CurrentAxes,'Pixel Cont.');xlabel(Dan_histFig.CurrentAxes,'Dan (cm^2/s)');
set(gca,'FontSize',13,'FontWeight','bold'); xlim([0 0.14]);
% print(Dan_histFig,[Datapath , '\','Dan_histFig'],'-dpng','-r300');
saveas(gca,'Dan_histFig.png');
Global.exportToPPTX('addpicture',Dan_histFig,'Position',[0 4.5 7 4.5]);

Dmean_histFig = figure; histogram(Dmean_map(Dmean_map>0),100); 
set(Dmean_histFig,'Name','Dmean_histFig'); title(Dmean_histFig.CurrentAxes,Dmean_title);ylabel(Dmean_histFig.CurrentAxes,'Pixel Cont.');xlabel(Dmean_histFig.CurrentAxes,'Dmean (cm^2/s)');
set(gca,'FontSize',13,'FontWeight','bold'); xlim([0 0.14]);
% print(Dmean_histFig,[Datapath , '\','Dmean_histFig'],'-dpng','-r300');
saveas(gca,'Dmean_histFig.png');
Global.exportToPPTX('addpicture',Dmean_histFig,'Position',[7.5 4.5 7 4.5]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% R/h/r/Lm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Global.exportToPPTX('addslide'); % R/h/r/Lm
Global.exportToPPTX('addpicture',R_Montage,'Position',[0 0.4 scale_fac scale_fac*(R_MontagePosition(4)/R_MontagePosition(3))]);
Global.exportToPPTX('addpicture',h_Montage,'Position',[0 2.4 scale_fac scale_fac*(h_MontagePosition(4)/h_MontagePosition(3))]);
Global.exportToPPTX('addpicture',r_Montage,'Position',[0 4.4 scale_fac scale_fac*(r_MontagePosition(4)/r_MontagePosition(3))]);
Global.exportToPPTX('addpicture',Lm_Montage,'Position',[0 6.4 scale_fac scale_fac*(Lm_MontagePosition(4)/Lm_MontagePosition(3))]);

R_title = sprintf('R-Mean= %3.2f ± %2.2f (um)',R_mean,R_std);
Global.exportToPPTX('addtext',R_title,'Position',[0 0 6 1],'Color','r','FontSize',20);
h_title = sprintf('h-Mean= %3.2f ± %2.2f (um)',h_mean,h_std);    
Global.exportToPPTX('addtext',h_title,'Position',[0 2 6 1],'Color','r','FontSize',20);
r_title = sprintf('r-Mean= %3.2f ± %2.2f (um)',r_mean,r_std);    
Global.exportToPPTX('addtext',r_title,'Position',[0 4 6 1],'Color','r','FontSize',20);
Lm_title = sprintf('Lm-Mean= %3.2f ± %2.2f (um)',Lm_mean,Lm_std);    
Global.exportToPPTX('addtext',Lm_title,'Position',[0 6 6 1],'Color','r','FontSize',20);

Global.exportToPPTX('addslide'); %R/h/r/Lm Histogram
R_histFig = figure; histogram(R_map(R_map>0),100); 
set(R_histFig,'Name','R_histFig'); title(R_histFig.CurrentAxes,R_title);ylabel(R_histFig.CurrentAxes,'Pixel Cont.');xlabel(R_histFig.CurrentAxes,'R (um)');
set(gca,'FontSize',13,'FontWeight','bold'); xlim([0 400]);
% print(R_histFig,[Datapath , '\','R_histFig'],'-dpng','-r300');
saveas(gca,'R_histFig.png');
Global.exportToPPTX('addpicture',R_histFig,'Position',[0 0 7 4.5]);

h_histFig = figure; histogram(h_map(h_map>0),100); 
set(h_histFig,'Name','h_histFig'); title(h_histFig.CurrentAxes,h_title);ylabel(h_histFig.CurrentAxes,'Pixel Cont.');xlabel(h_histFig.CurrentAxes,'h (um)');
set(gca,'FontSize',13,'FontWeight','bold'); xlim([0 400]);
% print(h_histFig,[Datapath , '\','h_histFig'],'-dpng','-r300');
saveas(gca,'h_histFig.png');
Global.exportToPPTX('addpicture',h_histFig,'Position',[7.5 0 7 4.5]);

r_histFig = figure; histogram(r_map(r_map>0),100); 
set(r_histFig,'Name','r_histFig'); title(r_histFig.CurrentAxes,r_title);ylabel(r_histFig.CurrentAxes,'Pixel Cont.');xlabel(r_histFig.CurrentAxes,'r (um)');
set(gca,'FontSize',13,'FontWeight','bold'); xlim([0 400]);
% print(r_histFig,[Datapath , '\','r_histFig'],'-dpng','-r300');
saveas(gca,'r_histFig.png');
Global.exportToPPTX('addpicture',r_histFig,'Position',[0 4.5 7 4.5]);

Lm_histFig = figure; histogram(Lm_map(Lm_map>0),100); 
set(Lm_histFig,'Name','Lm_histFig'); title(Lm_histFig.CurrentAxes,Lm_title);ylabel(Lm_histFig.CurrentAxes,'Pixel Cont.');xlabel(Lm_histFig.CurrentAxes,'Lm (um)');
set(gca,'FontSize',13,'FontWeight','bold'); xlim([0 400]);
% print(Lm_histFig,[Datapath , '\','Lm_histFig'],'-dpng','-r300');
saveas(gca,'Lm_histFig.png');
Global.exportToPPTX('addpicture',Lm_histFig,'Position',[7.5 4.5 7 4.5]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SVR/Na %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Global.exportToPPTX('addslide'); % SVR/Na
Global.exportToPPTX('addpicture',SVR_Montage,'Position',[0 0.4 scale_fac scale_fac*(SVR_MontagePosition(4)/SVR_MontagePosition(3))]);
Global.exportToPPTX('addpicture',Na_Montage,'Position',[0 2.4 scale_fac scale_fac*(Na_MontagePosition(4)/Na_MontagePosition(3))]);

SVR_title = sprintf('SVR-Mean= %3.2f ± %2.2f (cm^-1)',SVR_mean,SVR_std);
Global.exportToPPTX('addtext',SVR_title,'Position',[0 0 6 1],'Color','r','FontSize',20);
Na_title = sprintf('Na-Mean= %3.2f ± %2.2f (mm^-3)',Na_mean,Na_std);    
Global.exportToPPTX('addtext',Na_title,'Position',[0 2 6 1],'Color','r','FontSize',20);

SVR_histFig = figure; histogram(SVR_map(SVR_map>0),100); 
set(SVR_histFig,'Name','SVR_histFig'); title(SVR_histFig.CurrentAxes,SVR_title);ylabel(SVR_histFig.CurrentAxes,'Pixel Cont.');xlabel(SVR_histFig.CurrentAxes,'SVR (cm^-1)');
set(gca,'FontSize',13,'FontWeight','bold');
if strcmp(fitType,'human') == 1
    xlim([0 400]);
elseif strcmp(fitType,'animals') == 1
    xlim([0 2000]);
end
% print(SVR_histFig,[Datapath , '\','SVR_histFig'],'-dpng','-r300');
saveas(gca,'SVR_histFig.png');
Global.exportToPPTX('addpicture',SVR_histFig,'Position',[0 4.5 7 4.5]);

Na_histFig = figure; histogram(Na_map(Na_map>0),100); 
set(Na_histFig,'Name','Na_histFig'); title(Na_histFig.CurrentAxes,Na_title);ylabel(Na_histFig.CurrentAxes,'Pixel Cont.');xlabel(Na_histFig.CurrentAxes,'Na (mm^-3)');
set(gca,'FontSize',13,'FontWeight','bold'); 
if strcmp(fitType,'human') == 1
    xlim([0 400]);
elseif strcmp(fitType,'animals') == 1
    xlim([0 20000]);
end
% print(Na_histFig,[Datapath , '\','Na_histFig'],'-dpng','-r300');
saveas(gca,'Na_histFig.png');
Global.exportToPPTX('addpicture',Na_histFig,'Position',[7.5 4.5 7 4.5]);
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd(Datapath);
Global.exportToPPTX('save',fullfile(Datapath, ReportTitle));
Global.exportToPPTX('close');
fprintf('PowerPoint file has been saved\n');    
% test = [Datapath, ReportTitle];
%% plot residual
residual=Bayes_Residual3;
figure(1);
subplot(1,2,1);
chosen_res=sum(Bayes_Residual3,2);
chosen_res=chosen_res./max(chosen_res);
% chosen_res(chosen_res==0)=NaN;
histogram(chosen_res,100, 'normalization', 'pdf');
hold on

mu = mean(chosen_res);
sigma = std(chosen_res);
y = -max(chosen_res):0.01:max(chosen_res);
f = exp(-(y-mu).^2./(2*sigma^2))./(sigma*sqrt(2*pi));
plot(y,f,'LineWidth',1.5)
xlabel('Residual','FontSize',12,'FontWeight','bold','Color','k')
ylabel('Frequency','FontSize',12,'FontWeight','bold','Color','k')
set(gca,'FontSize',10)
xlim([-0.5 0.5])

%fit residual with data 
figure(1);
subplot(1,2,2);

fitted_data= Images.*final_mask;
fitted_data=reshape(fitted_data,[],1);
fitted_data = fitted_data(any(fitted_data,2),:);
fitted_residual=residual;
fitted_residual=reshape(fitted_residual,[],1);

xline=1:1:max(fitted_data);
yline=zeros(1,length(xline));
plot(fitted_data,fitted_residual,'b.','MarkerSize',10);
hold on 
plot(xline,yline,'r','MarkerSize',5);
xlim([0 max(fitted_data)])
ylim([-100 100])
xlabel('Fitted Data','FontSize',12,'FontWeight','bold','Color','k')
ylabel('Residual','FontSize',12,'FontWeight','bold','Color','k')
set(gca,'FontSize',10)
% print('CM_Residual','-r150','-dpng');
saveas(gca,'CM_Residual.png');

close all;
%% save maps in mat file
save_data=[Datapath , 'CM_Morphometry.mat'];
save(save_data,'DT_map','DL_map','Dan_map','Dmean_map','R_map','h_map','r_map','Lm_map','SVR_map','Na_map','So_map','MSR_map','Residual_map');   

%% 
fprintf('Lung Morphometry fit compeleted\n'); 
Run_Time=toc/60
end


% Please add updates here
% 
% 
% 
%
