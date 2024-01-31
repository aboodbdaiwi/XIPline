function [DDC_map,alpha_map,So_map,LmD_map,LungSEMMorphometrySummary]=SEM_Morphometry(Images,BNmask,noise_mask,bvalues,Do,Delta,Datapath,WinBUGSPath)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%% 129Xe stretched exponential model (SEM) Lung Morphometry  %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Requirements: WinBUGS program
%   (https://www.mrc-bsu.cam.ac.uk/software/bugs/the-bugs-project-winbugs/)
% 
%   Inputs:
%       Images = 4D images (x,y,slices,bvalue)
%       final_mask = lung mask, 3D (x,y,slices)
%       noise_mask = background mask without airways and artifacts, 3D (x,y,slices)
%       bvalues = vector: has all b-values (s/cm2)
%       Do = Xenon free diffusion (cm2/s): 0.14
%       delta = The diffusion time (s)
%       Datapath = path to the location of the data
%       WinBUGSPath = location to the winbugs.exe file
%
%   Outputs:
%       DDC_map = external Acinar duct radius 
%       alpha_map = alveolar sleeve depth 
%       So_map = internal Acinar duct radius map
%       LmD_map = mean linear intercept map

%   Package: https://github.com/aboodbdaiwi/HP129Xe_Analysis_App
%
%   Author: Abdullah S. Bdaiwi
%   Work email: abdullah.bdaiwi@cchmc.org
%   Personal email: abdaiwi89@gmail.com
%   Website: https://www.cincinnatichildrens.org/research/divisions/c/cpir
%
%   Please add updates at the end. Ex: 3/10/24 - ASB: update .... 
%
% This code is based on the stretched exponential model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
Datapath = char(Datapath);
Images=(Images./max(Images(:))).*100;
Images (Images < 0.0001) = 0.0001; % to avoid having zeros 

% get rid of any slice that has no signal (un-segmented)
slices=size(Images,3);
slice_checker=zeros(1,slices);
for check_slice=1:slices
    sl=BNmask(:,:,check_slice);
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
    Nfinal_mask(:,:,find_sl)=BNmask(:,:,choose_sl);
    Nnoise_mask(:,:,find_sl)=noise_mask(:,:,choose_sl);
end
% figure; imslice(Ndiffimg)
% figure; imslice(Nfinal_mask)
Ndiffimg = Images1;
BNmask = Nfinal_mask;
noise_mask = Nnoise_mask;

%normalize images
Ndiffimg=(Ndiffimg./max(Ndiffimg(:))).*100;
Ndiffimg (Ndiffimg < 0.0001) = 0.0001; % to avoid having zeros 

%% Bayesain DDC and alpha fit
% if strcmp(Bayesain_fit,'yes') == 1    
    fprintf( 'performing Bayesain Fitting for DDC and alpha...\n' );
    cd(Datapath);
    DDC_Model='DDC_alpha_Bayesian_Model_1D.txt';
    So_map=zeros(size(BNmask));
    DDC_map=zeros(size(BNmask));
    alpha_map=zeros(size(BNmask));
    sigma_n=zeros(size(BNmask,3),length(bvalues));
    Bayes_So = zeros(size(BNmask,3),size(BNmask,1)*size(BNmask,2));
    Bayes_DDC = zeros(size(BNmask,3),size(BNmask,1)*size(BNmask,2));
    Bayes_alpha = zeros(size(BNmask,3),size(BNmask,1)*size(BNmask,2));
    for slice_n=1:size(BNmask,3)
        M = zeros(size(BNmask,1),size(BNmask,2),length(bvalues));
        for n = 1:length(bvalues)
            noise_vec=Ndiffimg(:,:,slice_n,n);
            noise_vec(noise_mask(:,:,slice_n)==0)=[];
            std_slice_n(n) = std(noise_vec(:));
            sigma(n) = std_slice_n(1);
%             if std_slice_n(n) < 0.5
%                 std_slice_n(n) = 0.5;
%             end
            M(:,:,n)=((Ndiffimg(:,:,slice_n,n).*BNmask(:,:,slice_n)).^2)/(std_slice_n(1)^2);
        end
        M_Data_vec=reshape(M,[size(BNmask,1)*size(BNmask,2),length(bvalues)]);
        M = M_Data_vec(any(M_Data_vec,2),:);
        dataStruct= struct('M', M, 'b', bvalues', 'sigma', sigma);
        num_ones=BNmask(:,:,slice_n);
        num_ones=sum(num_ones(:)); % for loopping 
        DDC_Model_1D = DiffusionFunctions.DDC_alpha_Bayesian_Model_1D(num_ones,length(bvalues));

        %initializing
        a=0; b_int=Do;
        DDC_prior = a + (b_int-a).*rand(num_ones,1); %DDC prior  

        a=0.000001; b_int=1;
        alpha_prior = a + (b_int-a).*rand(num_ones,1); %alpha prior  

        So_prior=reshape(Ndiffimg(:,:,slice_n,1).*BNmask(:,:,slice_n),[size(BNmask,1)*size(BNmask,2),1]); % So %adc prior
        So_prior = So_prior(any(So_prior,2),:);  %So prior 
        initStructs =struct('So', So_prior','DDC', DDC_prior','alpha',alpha_prior);

        % Call WindBuges
        nchains  = 1; % How Many Chains?
        nburnin  = 1500; % How Many Burn-in Samples? 1000
        nsamples = 1000;  % How Many Recorded Samples? 500
        fprintf( 'Running WinBUGS for Bayesian Fit...\n' );
        [samples, stats] = DiffusionFunctions.matbugs(dataStruct, ...
            fullfile(Datapath , DDC_Model), ...
            'init', initStructs, ...
            'nChains', nchains, ...
            'view', 0, 'nburnin', nburnin, 'nsamples', nsamples, ...
            'thin', 1, 'DICstatus', 0, 'refreshrate',100,  ...
            'monitorParams', {'So','DDC','alpha'}, ...
            'Bugdir', WinBUGSPath);

        Bayes_So(slice_n,1:num_ones) =stats.mean.So;
        Bayes_DDC(slice_n,1:num_ones) =stats.mean.DDC;
        Bayes_alpha(slice_n,1:num_ones) =stats.mean.alpha;
        
    end

    pixel_location=find(BNmask); %to find pixel location to convert the array back to 3D image
    num_ones=sum(BNmask(:)); % for loopping 

    Bayes_So_2=Bayes_So';
    Bayes_So_2=reshape(Bayes_So_2,[size(BNmask,1)*size(BNmask,2)*size(BNmask,3),1]);
    Bayes_So_2 = Bayes_So_2(any(Bayes_So_2,2),:);

    Bayes_DDC_2=Bayes_DDC';
    Bayes_DDC_2=reshape(Bayes_DDC_2,[size(BNmask,1)*size(BNmask,2)*size(BNmask,3),1]);
    Bayes_DDC_2 = Bayes_DDC_2(any(Bayes_DDC_2,2),:);

    Bayes_alpha_2=Bayes_alpha';
    Bayes_alpha_2=reshape(Bayes_alpha_2,[size(BNmask,1)*size(BNmask,2)*size(BNmask,3),1]);
    Bayes_alpha_2 = Bayes_alpha_2(any(Bayes_alpha_2,2),:);

    %convert vector to 2D image 
    for bayes_loop=1:num_ones
        choosen_pixel=pixel_location(bayes_loop);
        So_map(choosen_pixel) = Bayes_So_2(bayes_loop);
        DDC_map(choosen_pixel) = Bayes_DDC_2(bayes_loop);
        alpha_map(choosen_pixel) = Bayes_alpha_2(bayes_loop);
    end 
% end
fprintf('Completed\n');  


%%  % DDC and alpha fit (Least square fit)
% DDC_map=zeros(size(BNmask));
% alpha_map=zeros(size(BNmask));
% So_map=zeros(size(BNmask));
% fun=@(x,xdata)abs((x(1)).*exp((-xdata.*x(2)).^x(3)));
% x0 = [100,0.035,0.5];
% lb = [0.01,0.01,0];
% ub = [120,0.14,1];
% for Kzs=size(BNmask,3)
%     for kys=1:size(BNmask,1)
%         for kxs=1:size(BNmask,2)      
%             Signal_inPixel = squeeze(Images(kys,kxs,Kzs,:));  
%             Signal_inPixel=Signal_inPixel';
%             H_D = lsqcurvefit(fun,x0,bvalues,Signal_inPixel,lb,ub); 
%             % H_D = lsqcurvefit(fun,x0,bvalues,S); 
%             So_map(kys,kxs,Kzs)= H_D(1);
%             DDC_map(kys,kxs,Kzs) = H_D(2);
%             alpha_map(kys,kxs,Kzs) = H_D(3);
%         end
%     end
% end
%% set up vectors distributions for LmD fitting 
fprintf('Calculating P(D) and LmD maps...\n');  
Do_vec = 0.001:0.0001:Do;
Delta = Delta/1000; 
LmD_vec=(sqrt(2.*Delta.*Do_vec)).*1e4;
alpha_vec=0.1:0.1:0.9;
C_vec=[0.89, 0.50, 0.35, 0.25, 0, 0.13, 0.22, 0.4015, 0.33];
B_vec=[0.145, 0.197, 0.243, 0.285, 0.382, 0.306, 0.360, 0.435, 0.700];

alpha_map_thurshed=alpha_map;
alpha_map_thurshed(alpha_map_thurshed<=0.1)=0.1;
alpha_map_thurshed(alpha_map_thurshed>=0.9)=0.9;
alpha_map_thurshed=alpha_map_thurshed.*BNmask;
alpha_map_thurshed=alpha_map_thurshed.*BNmask;

fD_delta1=(alpha_map_thurshed.*(0.5-alpha_map_thurshed))./(1-alpha_map_thurshed);
fD_delta1=fD_delta1.*BNmask;
fD_delta2=(alpha_map_thurshed.*(alpha_map_thurshed-0.5))./(1-alpha_map_thurshed);
fD_delta2=fD_delta2.*BNmask;
PD_power1=(1-alpha_map_thurshed./2)./(1-alpha_map_thurshed);
PD_power1=PD_power1.*BNmask;
PD_power2=alpha_map_thurshed./(1-alpha_map_thurshed);
PD_power2=PD_power2.*BNmask;

tau=1./DDC_map;
tau=tau.*BNmask;
tau(isinf(tau)|isnan(tau)) = 0;

% Preallocate memory  
B=zeros(size(BNmask));
C=zeros(size(BNmask));
fD=zeros(size(BNmask,1),size(BNmask,2),size(BNmask,3),length(Do_vec));
PD=zeros(size(fD));
% maxPD_loc=zeros(size(BNmask));
LmD_map=zeros(size(BNmask));
for k = 1:size(BNmask,3)
    for i = 1:size(BNmask,1)
        for j = 1:size(BNmask,2)
            if BNmask(i,j,k) == 1
                B(i,j,k) = interp1(alpha_vec,B_vec,alpha_map_thurshed(i,j,k));
                C(i,j,k) = interp1(alpha_vec,C_vec,alpha_map_thurshed(i,j,k));
                if alpha_map_thurshed(i,j,k)<=0.5
                    fD(i,j,k,:)=1./(  1+C(i,j,k).*((Do_vec.*tau(i,j,k)).^(fD_delta1(i,j,k)))  );
                elseif alpha_map_thurshed(i,j,k)>0.5
                    fD(i,j,k,:) = (1+C(i,j,k)*((Do_vec.*tau(i,j,k)).^(fD_delta2(i,j,k))));                
                end
                PD(i,j,k,:) = tau(i,j,k).*B(i,j,k).*(1./((Do_vec.* tau(i,j,k)).^PD_power1(i,j,k))).*...
                           (exp((-(1-alpha_map_thurshed(i,j,k)).*alpha_map_thurshed(i,j,k).^PD_power2(i,j,k))./((Do_vec.* tau(i,j,k)).^PD_power2(i,j,k)))).*...
                            squeeze(fD(i,j,k,:))';
    %             maxPD_loc(i,j)=find(PD(i,j,:)==max(PD(i,j,:)));       
    %             LmD_map(i,j)=LmD_vec(maxPD_loc(i,j));
                PD(isinf(PD)|isnan(PD)) = 0;
                PD_pixel = squeeze(PD(i,j,k,:));
                counts = PD_pixel(1:end-1);
                histfig=figure;
                H_LmD=histogram('BinEdges',LmD_vec,'BinCounts',counts); 
                counts = H_LmD.Values; close(histfig);
                LmD_map(i,j,k) = sum(LmD_vec(1:end-1) .* counts) / sum(counts);
            end      
        end
    end
end 
fprintf('Completed\n');  

%% Find mean and std of all maps
DDC_mean = round(mean(DDC_map(DDC_map>0)),3); 
alpha_mean = round(mean(alpha_map(alpha_map>0)),2); 
LmD_mean = round(mean(LmD_map(LmD_map>0))); 

DDC_std = round(std(DDC_map(DDC_map>0)),3); 
alpha_std = round(std(alpha_map(alpha_map>0)),2);
LmD_std = round(std(LmD_map(LmD_map>0))); 

LungSEMMorphometrySummary = table({'Mean';'SD'},[DDC_mean;DDC_std],[alpha_mean;alpha_std],[LmD_mean;LmD_std]);

%% 
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
%DDC_map 
Proton_Image = zeros(size(DDC_map));
for slice=1:size(DDC_map,3) %repeat for rest of slices
    [~,~] = Global.imoverlay(squeeze(abs(Proton_Image(:,:,slice))),squeeze(DDC_map(:,:,slice)),[],[0,0.99*max(Proton_Image(:))],cmap,1,gca);
    colormap(gca,cmap); caxis([0 0.14]);   
%     X = print('-RGBImage',['-r',num2str(size(DDC_map,2)/2)]);%2 inches
     Xdata = getframe(gcf);
     X = Xdata.cdata;     
    if (slice == 1)
        imwrite(X,[Datapath,'\DDC_map.tif'],'Description',strcat('Package Version: ', '1','; Cohort: ', 'test'));%write new/ overwrite tiff
    else
        imwrite(X,[Datapath,'\DDC_map.tif'],'WriteMode','append','Description',strcat('Package Version: ', '1','; Cohort: ', 'test'));%append tiff
    end
end

%alpha_map 
for slice=1:size(alpha_map,3) %repeat for rest of slices
    [~,~] = Global.imoverlay(squeeze(abs(Proton_Image(:,:,slice))),squeeze(alpha_map(:,:,slice)),[],[0,0.99*max(Proton_Image(:))],cmap,1,gca);
    colormap(gca,cmap); caxis([0 1]);   
%     X = print('-RGBImage',['-r',num2str(size(alpha_map,2)/2)]);%2 inches
     Xdata = getframe(gcf);
     X = Xdata.cdata;   
    if (slice == 1)
        imwrite(X,[Datapath,'\alpha_map.tif'],'Description',strcat('Package Version: ', '1','; Cohort: ', 'test'));%write new/ overwrite tiff
    else
        imwrite(X,[Datapath,'\alpha_map.tif'],'WriteMode','append','Description',strcat('Package Version: ', '1','; Cohort: ', 'test'));%append tiff
    end
end

%LmD_map 
for slice=1:size(LmD_map,3) %repeat for rest of slices
    [~,~] = Global.imoverlay(squeeze(abs(Proton_Image(:,:,slice))),squeeze(LmD_map(:,:,slice)),[],[0,0.99*max(Proton_Image(:))],cmap,1,gca);
    colormap(gca,cmap); caxis([0 400]);   
%     X = print('-RGBImage',['-r',num2str(size(LmD_map,2)/2)]);%2 inches
     Xdata = getframe(gcf);
     X = Xdata.cdata;     
    if (slice == 1)
        imwrite(X,[Datapath,'\LmD_map.tif'],'Description',strcat('Package Version: ', '1','; Cohort: ', 'test'));%write new/ overwrite tiff
    else
        imwrite(X,[Datapath,'\LmD_map.tif'],'WriteMode','append','Description',strcat('Package Version: ', '1','; Cohort: ', 'test'));%append tiff
    end
end

disp('Saving Tiff Completed.')
close all;

%% plot and save results
cmap = colormap (jet);
cmap(1,:)=0;
close all

fprintf('Save results in PowerPoint\n');   
%DDC maps
DDC_Montage = figure('Name','DDC Images');set(DDC_Montage,'WindowState','minimized');
montage(reshape(DDC_map,[size(DDC_map,1), size(DDC_map,2), 1, size(DDC_map,3)]),...
    'Size',[1 size(DDC_map,3)],'DisplayRange',[0 1]);
set(gca,'units','pixels'); % set the axes units to pixels
x = get(gca,'position'); % get the position of the axes
set(gcf,'units','pixels'); % set the figure units to pixels
y = get(gcf,'position'); % get the figure position
set(gcf,'position',[y(1) y(2) x(3) x(4)])% set the position of the figure to the length and width of the axes
set(gca,'units','normalized','position',[0 0 1 1]) % set the axes units to pixels
colormap(cmap); caxis([0 Do]); 
xx = colorbar;
xx.Label.String = ' (cm^2/S)';
DDC_MontagePosition=get(gcf,'position');

%Alpha maps
Alpha_Montage = figure('Name','Alpha Images');set(Alpha_Montage,'WindowState','minimized');
montage(reshape(alpha_map,[size(alpha_map,1), size(alpha_map,2), 1, size(alpha_map,3)]),...
    'Size',[1 size(alpha_map,3)],'DisplayRange',[0 1]);
set(gca,'units','pixels'); % set the axes units to pixels
x = get(gca,'position'); % get the position of the axes
set(gcf,'units','pixels'); % set the figure units to pixels
y = get(gcf,'position'); % get the figure position
set(gcf,'position',[y(1) y(2) x(3) x(4)])% set the position of the figure to the length and width of the axes
set(gca,'units','normalized','position',[0 0 1 1]) % set the axes units to pixels
colormap(cmap); caxis([0 1]);  colorbar;
Alpha_MontagePosition=get(gcf,'position');
 
%LmD maps
LmD_Montage = figure('Name','LmD Images');set(LmD_Montage,'WindowState','minimized');
montage(reshape(LmD_map,[size(LmD_map,1), size(LmD_map,2), 1, size(LmD_map,3)]),...
    'Size',[1 size(LmD_map,3)],'DisplayRange',[0 1]);
set(gca,'units','pixels'); % set the axes units to pixels
x = get(gca,'position'); % get the position of the axes
set(gcf,'units','pixels'); % set the figure units to pixels
y = get(gcf,'position'); % get the figure position
set(gcf,'position',[y(1) y(2) x(3) x(4)])% set the position of the figure to the length and width of the axes
set(gca,'units','normalized','position',[0 0 1 1]) % set the axes units to pixels
colormap(cmap); caxis([0 LmD_vec(end)]); colorbar;
xx = colorbar;
xx.Label.String = ' (\mum)';
LmD_MontagePosition=get(gcf,'position');

nslices = size(LmD_map,3);
if nslices > 8
    nslices = 8;
end
scale_fac = 2*nslices;
% save maps in ppt file
ReportTitle='Diffusion_Analysis';
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
%%%%%%%%%%%%%%%%%%%%%%%%%% DDC/alpha/LmD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Global.exportToPPTX('addslide'); % DDC/alpha/LmD
Global.exportToPPTX('addpicture',DDC_Montage,'Position',[0.5 0.8 scale_fac scale_fac*(DDC_MontagePosition(4)/DDC_MontagePosition(3))]);
Global.exportToPPTX('addpicture',Alpha_Montage,'Position',[0.5 3.3 scale_fac scale_fac*(Alpha_MontagePosition(4)/Alpha_MontagePosition(3))]);
Global.exportToPPTX('addpicture',LmD_Montage,'Position',[0.5 5.8 scale_fac scale_fac*(LmD_MontagePosition(4)/LmD_MontagePosition(3))]);

Global.exportToPPTX('addtext',sprintf('SEM-Morphometry Analysis'),'Position',[6 0 5 1],'Color','b','FontSize',25);
DDC_title = sprintf('DDC-Mean= %1.3f ± %1.3f (cm^2/s)',DDC_mean,DDC_std);
Global.exportToPPTX('addtext',DDC_title,'Position',[0.5 0.5 6 1],'Color','r','FontSize',20);
alpha_title = sprintf('Alpha-Mean= %1.3f ± %1.3f ',alpha_mean,alpha_std);    
Global.exportToPPTX('addtext',alpha_title,'Position',[0.5 3 6 1],'Color','r','FontSize',20);
LmD_title = sprintf('LmD-Mean= %1.1f ± %1.1f (um)',LmD_mean,LmD_std);    
Global.exportToPPTX('addtext',LmD_title,'Position',[0.5 5.5 6 1],'Color','r','FontSize',20);

Global.exportToPPTX('addslide'); % DDC/alpha/LmD histograms
DDC_histFig = figure; histogram(DDC_map(DDC_map>0),100); 
set(DDC_histFig,'Name','DDC_histFig'); title(DDC_histFig.CurrentAxes,DDC_title);ylabel(DDC_histFig.CurrentAxes,'Pixel Cont.');xlabel(DDC_histFig.CurrentAxes,'DDC (cm^2/s)');
set(gca,'FontSize',13,'FontWeight','bold'); xlim([0 0.14]);
% print(DDC_histFig,[Datapath , '\','DCC_histFig'],'-dpng','-r300');
saveas(gca,'DDC_histFig.png');
Global.exportToPPTX('addpicture',DDC_histFig,'Position',[0 0 7 4.5]);

alpha_histFig = figure; histogram(alpha_map(alpha_map>0),100); 
set(alpha_histFig,'Name','alpha_histFig'); title(alpha_histFig.CurrentAxes,alpha_title);ylabel(alpha_histFig.CurrentAxes,'Pixel Cont.');xlabel(alpha_histFig.CurrentAxes,'alpha');
set(gca,'FontSize',13,'FontWeight','bold'); xlim([0 1]);
% print(alpha_histFig,[Datapath , '\','alpha_histFig'],'-dpng','-r300');
saveas(gca,'alpha_histFig.png');
Global.exportToPPTX('addpicture',alpha_histFig,'Position',[7.5 0 7 4.5]);

LmD_histFig = figure; histogram(LmD_map(LmD_map>0),100); 
set(LmD_histFig,'Name','LmD_histFig'); title(LmD_histFig.CurrentAxes,LmD_title);ylabel(LmD_histFig.CurrentAxes,'Pixel Cont.');xlabel(LmD_histFig.CurrentAxes,'LmD (um)');
set(gca,'FontSize',13,'FontWeight','bold'); xlim([0 400]);
% print(LmD_histFig,[Datapath , '\','LmD_histFig'],'-dpng','-r300');
saveas(gca,'LmD_histFig.png');
Global.exportToPPTX('addpicture',LmD_histFig,'Position',[0 4.5 7 4.5]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd(Datapath);
Global.exportToPPTX('save',fullfile(Datapath , ReportTitle));
Global.exportToPPTX('close');
fprintf('PowerPoint file has been saved\n');  

close all;
%% save maps in mat file
save_data=[Datapath , 'SEM_Morphometry.mat'];
save(save_data,'DDC_map','alpha_map','So_map','LmD_map');   

%% 
fprintf('SEM Lung Morphometry fit compeleted\n'); 
Run_Time=toc/60
end

% Please add updates here
% 
% 
% 
%