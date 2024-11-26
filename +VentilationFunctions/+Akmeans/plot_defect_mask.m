function plot_defect_mask(background_data,mask_data,temp_handles,plotting_range,title_string)
%
% Create image overlay montage
% inputs:
% background_data  - bottom layer typically in grayscale
% mask_data        - upper layer typically in colored contours, i.e.,
%                    defect by lobes
% temp_handles     .lungmask  -binary mask required;
%                  .lobemask  - optional
%                  .lobemaskRUL - optional                
%                  .lobemaskRML - optional                
%                  .lobemaskRLL - optional                
%                  .lobemaskLUL - optional                
%                  .lobemaskLLL  - optional               
% plotting_range (optional) .x -numerical vector
%                           .y -numerical vector
%                           
% title_string (optional)
%
% Copyright W. Zha@2014

[temp_handles,~,lobe_label_list,valid_lobe_labels]=identify_lobes(temp_handles);
[nRows,nCols,nSlices] = size(background_data);
if nargin<4
    plotting_range.x =1:nRows;
    plotting_range.y =1:nCols;
    title_string=[];
    
elseif nargin<5
    title_string=[];
end
fig_col=5;
fig_row=ceil(nSlices/fig_col);

lobe_color_list = {'green','yellow','pink','purple','orange'};
coloreddata = zeros(nRows,nCols,3,nSlices);


for sl=1:nSlices
    
    greydata =background_data(:,:,sl);
    greydata(isnan(greydata))=0;
    
    this_defect_mask =logical(mask_data(:,:,sl));
    minIntensity = min(greydata(:));
    maxIntensity = max(greydata(:));
    if maxIntensity==0
        maxIntensity=1;
    end
    
    color.red = zeros(nRows,nCols);
    for n_lobe = 1:numel(lobe_color_list)
        this_color = lobe_color_list{n_lobe};
        color.(this_color)=zeros(nRows,nCols);
    end
    if isempty(valid_lobe_labels)
        color.blue = edge(uint8(this_defect_mask));

    else
        color.blue = zeros(nRows,nCols);
        
        
        for n_lobe = 1:numel(lobe_label_list)
            this_color = lobe_color_list{n_lobe};
            this_lobe = lobe_label_list{n_lobe};
            lobemask_name = sprintf('lobemask%s',upper(this_lobe));
            this_lobe_defect_mask = zeros(nRows,nCols);
            if isfield(temp_handles,'slice_range')
                lobe.(lobemask_name)=temp_handles.(lobemask_name)(:,:,temp_handles.slice_range);
            else
                lobe.(lobemask_name)=temp_handles.(lobemask_name);
            end
            sl_no=size(lobe.(lobemask_name),3);
            if sl_no<nSlices
                lobemask = zeros(nRows,nCols,nSlices);
                lobemask(:,:,1:sl_no)=lobe.(lobemask_name);
            else
                lobemask = lobe.(lobemask_name);
            end
            this_lobe_defect_mask(this_defect_mask & lobemask(:,:,sl))=1;
            color.(this_color) = edge(uint8(this_lobe_defect_mask));
        end
        
    end
    coloreddata(:,:,:,sl)=plentyofcolors(((greydata-minIntensity)./(maxIntensity-minIntensity)),...
        color.red,color.blue,color.green,color.yellow,color.pink,color.purple,color.orange);
    subaxis(fig_row,fig_col,sl,'SpacingVert',0,'MR',0,'SpacingHoriz',0,'Spacing',0); 
    imshow(coloreddata(plotting_range.y,plotting_range.x,:,sl),[]);

    if sl==3
        title(title_string,'FontSize',16,'FontWeight','bold')
    end
end

