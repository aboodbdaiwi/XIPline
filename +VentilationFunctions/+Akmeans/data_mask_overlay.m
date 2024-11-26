function data_mask_overlay(background_data,mask_data,title_string)

% Wei Zha @ Feb, 2016
if nargin<3
    title_string='';
end
[nRows,nCols,nSlices]=size(background_data);
plotting_range=get_plotting_range(mask_data);
color_list={'red','yellow','blue','pink','orange','purple'};
coloreddata=zeros(nRows,nCols,3,nSlices);
for nc=1:numel(color_list)
    color.(color_list{nc})=zeros(nRows,nCols);
end
nSlices = min(nSlices,size(mask_data,3));
fig_col=5;
fig_row=ceil(nSlices/fig_col);
for sl=1:nSlices
    
    greydata =background_data(:,:,sl);
    greydata(isnan(greydata))=0;
    
    slicemask =logical(mask_data(:,:,sl));
    minIntensity = min(greydata(:));
    maxIntensity = max(greydata(:));
    if maxIntensity==0
        maxIntensity=1;
    end
    
    color.green=edge(slicemask);


    coloreddata(:,:,:,sl)=plentyofcolors(((greydata-minIntensity)./(maxIntensity-minIntensity)),...
        color.red,color.blue,color.green,color.yellow,color.pink,color.purple,color.orange);

    subaxis(fig_row,fig_col,sl,'SpacingVert',0,'MR',0,'SpacingHoriz',0,'Spacing',0); 

    imshow(coloreddata(plotting_range.y,plotting_range.x,:,sl),[]);brighten(.3);%colormap('jet');%axis off;

    if sl==3
        title(title_string,'FontSize',16,'FontWeight','bold')
    end
end

