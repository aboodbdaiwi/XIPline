
function [Diffusion] = CM_MorphometryFit(Diffusion,MainInput)

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

%   Package: https://github.com/aboodbdaiwi/HP129Xe_Analysis_App
%
%   Author: Abdullah S. Bdaiwi
%   Work email: abdullah.bdaiwi@cchmc.org
%   Personal email: abdaiwi89@gmail.com
%   Website: https://www.cincinnatichildrens.org/research/divisions/c/cpir
%
%   Please add updates at the end. Ex: 3/10/24 - ASB: updated flip angle caluclation for spiral diffusion

% This code is based on the cylinder model by A. L. Sukstanskii1 and D. A. Yablonskiy
% (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3184317/)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Images = Diffusion.Image;
final_mask =  Diffusion.lung_mask;
noise_mask = Diffusion.noise_mask;
fitType = Diffusion.MorphometryAnalysisType;
bvalues = Diffusion.b_values;
Do = Diffusion.Do;
delta = Diffusion.Delta;
Datapath = Diffusion.outputpath;
WinBUGSPath = Diffusion.WinBUGSPath;

tic
% denoise images using bm3d
% Images = (Images - min(Images(:)))/(max(Images(:)) - min(Images(:)));
% for i = 1:size(Images,3)
%     for j = 1:size(Images,4)
%         Images(:,:,i,j) = Global.bm3d.BM3D(squeeze(Images(:,:,i,j)), 0.02);
%     end
% end
%  figure; imslice(Images)
%normalize images
Images = (Images - min(Images(:)))/(max(Images(:)) - min(Images(:))).*100;
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

bvalues(1) = 1e-5;

So=zeros(size(Images,3),size(Images,1)*size(Images,2));
R=zeros(size(Images,3),size(Images,1)*size(Images,2));
r=zeros(size(Images,3),size(Images,1)*size(Images,2));
%% 

disp('start morphometery fitting, please wait......')
for i=1:size(Images,3)
    disp(['slice = ',num2str(i)]);
    Image_toFit= squeeze(Images(:,:,i,:)); 
%     Image_toFit=reshape(Image_toFit,[size(Image_toFit,1),size(Image_toFit,2),nbvalues]);
    Image_toFit=(Image_toFit./max(Image_toFit(:)));    

    final_mask_slice=final_mask(:,:,i); 
    noise_mask_slice=noise_mask(:,:,i);

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

    end
    SNR_vec = signal_n_avg ./ mean(std_noise_n); %signal to noise ratio
    sigma = mean(std_noise_n);
    %%%%%%%%%%%%%%%%%%%%%   R and h fit  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    M1 = zeros(size(Image_toFit,1),size(Image_toFit,2),nbvalues);
    for n=1:nbvalues
        M1(:,:,n)=(Image_toFit(:,:,n).*final_mask_slice);
    end

    M_Data_vec=reshape(M1,[size(Image_toFit,1)*size(Image_toFit,2),nbvalues]);
    M = M_Data_vec(any(M_Data_vec,2),:);

    num_ones=final_mask_slice;
    num_ones=sum(num_ones(:)); % for loopping    
    S0_est = zeros(num_ones,1);
    R_est = size(S0_est);
    r_est = size(S0_est);

    weights=1; 
    parfor j= 1:num_ones 
        if strcmp(fitType,'human') == 1
            lb=[0.0200 0.0090 0.5*M(j,1)]; ub=[0.0650 0.0500 1.5*M(j,1)];
            initialvalues = [0.0300 0.0140 M(j,1)];
            [eestm,fval(j)] = DiffusionFunctions.MLfitRr2(M(j,:),bvalues,sigma^2,initialvalues,"cylRandr",weights,lb,ub,Do,delta);       % This is the only code line needed + function MLfitconloc1 
        elseif strcmp(fitType,'animals') == 1
            % lb=[0.00001 0.005 0.5]; ub=[M(j,1)+3*sqrt(s2) 0.12 1.1];
            lb=[0.060 0.0010 0.5*M(j,1)]; ub=[0.0140 0.0090 1.5*M(j,1)];
            initialvalues=[0.0100 0.0050 M(j,1)];
            [eestm,fval(j)] = DiffusionFunctions.MLfitRr2(M(j,:),bvalues,sigma^2,initialvalues,"cylRandrAnimal",weights,lb,ub,Do,delta);       % This is the only code line needed + function MLfitconloc1 
        end
        R_est(j)=eestm(1);
        r_est(j)=eestm(2);
        S0_est(j)=eestm(3);
    end

    So(i,1:num_ones) = S0_est;
    R(i,1:num_ones) = R_est*10000;
    r(i,1:num_ones) = r_est*10000;
end

%%%%%%%%%%%%%%%%%%%%convert vector to 2D image  %%%%%%%%%%%%%%%%%%%%%%%%%%% 
pixel_location=find(final_mask); %to find pixel location 
num_ones2=sum(final_mask(:)); 

Bayes_r_2=r';
Bayes_r_2=reshape(Bayes_r_2,[size(final_mask,1)*size(final_mask,2)*size(final_mask,3),1]);
Bayes_r_2 = Bayes_r_2(any(Bayes_r_2,2),:);

Bayes_So_from_R2=So';
Bayes_So_from_R2=reshape(Bayes_So_from_R2,[size(final_mask,1)*size(final_mask,2)*size(final_mask,3),1]);
Bayes_So_from_R2 = Bayes_So_from_R2(any(Bayes_So_from_R2,2),:);

Bayes_R_2=R';
Bayes_R_2=reshape(Bayes_R_2,[size(final_mask,1)*size(final_mask,2)*size(final_mask,3),1]);
Bayes_R_2 = Bayes_R_2(any(Bayes_R_2,2),:);

r_map=zeros(size(final_mask));
So_map=zeros(size(final_mask));
R_map=zeros(size(final_mask));

for loc = 1:numel(pixel_location)
    if loc <= num_ones2
        choosen_pixel = pixel_location(loc);
        r_map(choosen_pixel) = Bayes_r_2(loc);
        So_map(choosen_pixel) = Bayes_So_from_R2(loc);
        R_map(choosen_pixel) = Bayes_R_2(loc);
    end
end 
%h_map = h_map.*10;
% figure; imslice(R_map); colormap(jet)
% figure; imslice(r_map); colormap(jet)

% Check if a parallel pool is open
if ~isempty(gcp('nocreate'))
    % If a parallel pool is open, delete it
    delete(gcp('nocreate'));
end
% figure; imslice(R_map); colormap(jet); colorbar; clim([0 600]);
%%%%%%%%%%%%%%%%%% DL, r, L, SVR, Lm and Na  %%%%%%

h_map=R_map-r_map;
L_map = 0.765.*R_map; 

SVR_map=((2*pi.*R_map.*L_map)+2*pi.*((R_map.^2)-(r_map.^2))+16.*h_map.*L_map)./(pi.*(R_map.^2).*L_map);
Lm_map=pi.*L_map.*R_map.^2./(0.5.*pi.*R_map.*L_map+0.5.*pi.*h_map.*(2.*R_map-h_map)+4.*h_map.*L_map);
%Lm_map =4./SVR_map;
Na_map = 1./(pi.*(R_map.^2).*L_map);

SVR_map=SVR_map./Runit;% (in cm-1) 
Na_map=Na_map.*10.0e+9; % in mm^-3
SVR_map(SVR_map>1000) = 0;
Na_map(Na_map>1000) = 0;

SVR_map(isinf(SVR_map)|isnan(SVR_map))=0;
Na_map(isinf(Na_map)|isnan(Na_map))=0;


%% Find mean and std of all maps
R_mean = round(mean(R_map(R_map>0)));
h_mean = round(mean(h_map(h_map>0))); 
r_mean = round(mean(r_map(r_map>0))); 
Lm_mean = round(mean(Lm_map(Lm_map>0))); 
SVR_mean = round(mean(SVR_map(SVR_map>0)));
Na_mean = round(mean(Na_map(Na_map>0))); 
 
R_std = round(std(R_map(R_map>0))); 
h_std = round(std(h_map(h_map>0)));
r_std = round(std(r_map(r_map>0))); 
Lm_std = round(std(Lm_map(Lm_map>0))); 
SVR_std = round(std(SVR_map(SVR_map>0))); 
Na_std = round(std(Na_map(Na_map>0))); 

LungCMMorphometrySummary = table({'Mean';'SD'},[R_mean;R_std],[h_mean;h_std],[r_mean;r_std],[Lm_mean;Lm_std],...
    [SVR_mean;SVR_std],[Na_mean;Na_std]);
headers=({'summary','R(um)', 'h(um)', 'r(um)', 'Lm(um)',...
    'SVR(cm^-1)', 'Na(mm^-3)'});
LungCMMorphometrySummary.Properties.VariableNames = headers;
writetable(LungCMMorphometrySummary,[Datapath,'\CM_MorphometrySummary.xlsx'],'Sheet',1)
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
Proton_Image = zeros(size(R_map));

%R_map 
for slice=1:size(R_map,3) %repeat for rest of slices
    [~,~] = Global.imoverlay(squeeze(abs(Proton_Image(:,:,slice))),squeeze(R_map(:,:,slice)),[],[0,0.99*max(Proton_Image(:))],cmap,1,gca);
    colormap(gca,cmap); clim([0 400]);   
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
    colormap(gca,cmap); clim([0 400]);   
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
    colormap(gca,cmap); clim([0 400]);   
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
    colormap(gca,cmap); clim([0 400]);   
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
    colormap(gca,cmap); clim([0 400]);   
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
    colormap(gca,cmap); clim([0 400]);   
%     X = print('-RGBImage',['-r',num2str(size(Na_map,2)/2)]);%2 inches
     Xdata = getframe(gcf);
     X = Xdata.cdata;     
    if (slice == 1)
        imwrite(X,[Datapath,'\Na_map.tif'],'Description',strcat('Package Version: ', '1','; Cohort: ', 'test'));%write new/ overwrite tiff
    else
        imwrite(X,[Datapath,'\Na_map.tif'],'WriteMode','append','Description',strcat('Package Version: ', '1','; Cohort: ', 'test'));%append tiff
    end
end

close all;

%% plot and save results
cmap = colormap (jet);
cmap(1,:)=0;
close all
colorbarlimit = 400;
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
    colormap(cmap); clim([0 colorbarlimit]); colorbar;
elseif strcmp(fitType,'animals') == 1
    colormap(cmap); clim([0 400]); colorbar;
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
    colormap(cmap); clim([0 colorbarlimit]); colorbar;
elseif strcmp(fitType,'animals') == 1
    colormap(cmap); clim([0 150]); colorbar;
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
    colormap(cmap); clim([0 colorbarlimit]); colorbar;
elseif strcmp(fitType,'animals') == 1
    colormap(cmap); clim([0 150]); colorbar;
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
    colormap(cmap); clim([0 colorbarlimit]); colorbar;
elseif strcmp(fitType,'animals') == 1
    colormap(cmap); clim([0 150]); colorbar;
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
    colormap(cmap); clim([0 colorbarlimit]); colorbar;
elseif strcmp(fitType,'animals') == 1
    colormap(cmap); clim([0 2000]); colorbar;
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
    colormap(cmap); clim([0 colorbarlimit]); colorbar;
elseif strcmp(fitType,'animals') == 1
    colormap(cmap); clim([0 20000]); colorbar;
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% R/h/r/Lm %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Global.exportToPPTX('addslide'); % R/h/r/Lm
Global.exportToPPTX('addtext',sprintf('CM-Morphometry Analysis'),'Position',[6 0 5 1],'Color','b','FontSize',25);

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
set(gca,'FontSize',13,'FontWeight','bold'); xlim([0 colorbarlimit]);
% print(R_histFig,[Datapath , '\','R_histFig'],'-dpng','-r300');
saveas(gca,'R_histFig.png');
Global.exportToPPTX('addpicture',R_histFig,'Position',[0 0 7 4.5]);

h_histFig = figure; histogram(h_map(h_map>0),100); 
set(h_histFig,'Name','h_histFig'); title(h_histFig.CurrentAxes,h_title);ylabel(h_histFig.CurrentAxes,'Pixel Cont.');xlabel(h_histFig.CurrentAxes,'h (um)');
set(gca,'FontSize',13,'FontWeight','bold'); xlim([0 colorbarlimit]);
% print(h_histFig,[Datapath , '\','h_histFig'],'-dpng','-r300');
saveas(gca,'h_histFig.png');
Global.exportToPPTX('addpicture',h_histFig,'Position',[7.5 0 7 4.5]);

r_histFig = figure; histogram(r_map(r_map>0),100); 
set(r_histFig,'Name','r_histFig'); title(r_histFig.CurrentAxes,r_title);ylabel(r_histFig.CurrentAxes,'Pixel Cont.');xlabel(r_histFig.CurrentAxes,'r (um)');
set(gca,'FontSize',13,'FontWeight','bold'); xlim([0 colorbarlimit]);
% print(r_histFig,[Datapath , '\','r_histFig'],'-dpng','-r300');
saveas(gca,'r_histFig2.png');
Global.exportToPPTX('addpicture',r_histFig,'Position',[0 4.5 7 4.5]);

Lm_histFig = figure; histogram(Lm_map(Lm_map>0),100); 
set(Lm_histFig,'Name','Lm_histFig'); title(Lm_histFig.CurrentAxes,Lm_title);ylabel(Lm_histFig.CurrentAxes,'Pixel Cont.');xlabel(Lm_histFig.CurrentAxes,'Lm (um)');
set(gca,'FontSize',13,'FontWeight','bold'); xlim([0 colorbarlimit]);
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


SVR_histFig = figure; histogram(SVR_map(SVR_map>0)); 
set(SVR_histFig,'Name','SVR_histFig'); title(SVR_histFig.CurrentAxes,SVR_title);ylabel(SVR_histFig.CurrentAxes,'Pixel Cont.');xlabel(SVR_histFig.CurrentAxes,'SVR (cm^-1)');
set(gca,'FontSize',13,'FontWeight','bold');
if strcmp(fitType,'human') == 1
    xlim([0 1000]);
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
    xlim([0 1000]);
elseif strcmp(fitType,'animals') == 1
    xlim([0 2000]);
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

close all;

% store result
Diffusion.R_map = R_map;
Diffusion.h_map = h_map;
Diffusion.r_map = r_map;
Diffusion.Lm_map = Lm_map;
Diffusion.SVR_map = SVR_map;
Diffusion.Na_map = Na_map;
Diffusion.So_map = So_map;

Diffusion.R_mean = LungCMMorphometrySummary{1,2};
Diffusion.h_mean = LungCMMorphometrySummary{1,3};
Diffusion.r_mean = LungCMMorphometrySummary{1,4};
Diffusion.Lm_mean = LungCMMorphometrySummary{1,5};
Diffusion.SVR_mean = LungCMMorphometrySummary{1,6};
Diffusion.Na_mean = LungCMMorphometrySummary{1,7};

Diffusion.R_std = LungCMMorphometrySummary{2,2};
Diffusion.h_std = LungCMMorphometrySummary{2,3};
Diffusion.r_std = LungCMMorphometrySummary{2,4};
Diffusion.Lm_std = LungCMMorphometrySummary{2,5};
Diffusion.SVR_std = LungCMMorphometrySummary{2,6};
Diffusion.Na_std = LungCMMorphometrySummary{2,7};
%% save maps in mat file
save_data=[Datapath , '\CM_Morphometry.mat'];
save(save_data,'R_map','h_map','r_map','Lm_map','SVR_map','Na_map','So_map');   
%% 
fprintf('Lung Morphometry fit compeleted\n'); 
Run_Time=toc/60
end


% Please add updates here
% 
% 
% 
%
