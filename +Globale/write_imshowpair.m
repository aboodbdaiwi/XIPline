
function write_imshowpair (Image1,Image2,DataLocation)
%Tiffs 
% mkdir([DataLocation '\Gax Exchange Analysis']);
% outputpath = [DataLocation '\Gax Exchange Analysis'];
outputpath = DataLocation;
cd(outputpath) 

tiff = figure('MenuBar','none','ToolBar','none','DockControls','off','Resize','off','WindowState','minimized');%figure for tiffs
ax1 = axes('Parent',tiff);ax2 = axes('Parent',tiff);%make axis for both images
set(ax1,'Visible','off');set(ax2,'Visible','off');%turn off axis
set(ax1,'units','inches');set(ax2,'units','inches');%make axis units inches
set(ax1,'position',[0 0 2 2]);set(ax2,'position',[0 0 2 2]);%make axis same as image
set(gcf,'units','inches'); % set the figure units to pixels
set(gcf,'position',[1 1 2 2])% set the position of the figure to axes
disp('Saving Images Tiff...')

for slice=1:size(Image2,3) %repeat for rest of slices
%     Xe = imshow(fixedVolume(:,:,slice),[1,6],'Parent',ax2);%plot Xe image
%     set(Xe, 'AlphaData', ProtonMaskRegistered(:,:,slice));%mask Xe image
%     colormap(ax2,SixBinMap);%change colors
%     imshowpair(movingRegisteredVolume(:,:,10), fixedVolume(:,:,10));
    imshowpair(Image1(:,:,slice), Image2(:,:,slice),'ColorChannels','green-magenta'); %'red-cyan' or 'green-magenta'
    set(gca,'units','pixels'); % set the axes units to pixels
    x = get(gca,'position'); % get the position of the axes
    set(gcf,'units','pixels'); % set the figure units to pixels
    y = get(gcf,'position'); % get the figure position
    set(gcf,'position',[y(1) y(2) x(3) x(4)])% set the position of the figure to the length and width of the axes
    set(gca,'units','normalized','position',[0 0 1 1]) % set the axes units to pixels
    X = print('-RGBImage',['-r',num2str(0)]);%2 inches
    X = imresize(X,[size(Image1,1) size(Image1,2)]);    
    if (slice == 1)
        imwrite(X,[outputpath,'\ProtonRegistered.tif']);%write new/ overwrite tiff
    else
        imwrite(X,[outputpath,'\ProtonRegistered.tif'],'WriteMode','append');%append tiff
    end

end
close all;
disp('Completed...')