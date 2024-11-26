function [h_figure,overlay_data] = proton_gas_overlay(proton,gas,title_string,slice_range,plotting_range,type)
% Wei Zha @ 3/7/2015
[nRows,nCols,nSlices]=size(proton);
% Check if the size of proton and gas are different
if ~isequal(size(proton), size(gas))
    % Resize proton to match the size of gas in all three dimensions
    proton = imresize3(proton, size(gas));
else
    % If sizes are already the same, no need to resize
    proton = proton;
end

if nargin<3
    slice_range = 1:nSlices;
    title_string =' ';
    plotting_range=get_plotting_range(gas);
    type = 'overlay';
    
elseif nargin <4
    plotting_range=get_plotting_range(gas);
    type = 'overlay';
    
    if isfield(plotting_range,'z')
        slice_range = plotting_range.z;
    else
        slice_range = 1:nSlices;
    end
elseif nargin == 4
    plotting_range=get_plotting_range(gas);
    type = 'overlay';
    %     if isfield(plotting_range,'z')
    %         slice_range = plotting_range.z;
    %     else
    %         slice_range = 1:nSlices;
    %     end
    
end
if isempty(slice_range)
    slice_range = 1:nSlices;
end
if isempty(plotting_range)
    plotting_range=get_plotting_range(gas)
    
end
if size(slice_range,1)~=1
    slice_range = slice_range.';
end
nPlots=length(slice_range);
nCols_subplot = min(7,nPlots);
fig_col=min(nCols_subplot,nPlots);
fig_row=ceil(nPlots/fig_col);
if fig_row*fig_col == nPlots
    fig_col_real=fig_col+1; 
else
    fig_col_real=fig_col;
end
% rgb=zeros(nRows,nCols,3);
rgb = zeros(length(plotting_range.y),length(plotting_range.x),3);
overlay_data = zeros(nRows, nCols, 3, nSlices);
plot_count=0;

gas(isnan(gas))=0;
is_ventilation_mask = length(unique(gas))>2 & length(unique(gas))<=7;
is_defect_mask = length(unique(gas))<=2;

h_figure = figure;
h_figure.Color = 'black';
% h_figure.Position = [51,506,1700,400];
h_figure.Position = [51,506,800,400];
% h_figure.Position = [20,100,600,400];

VentValues = unique(gas);
VentValues(VentValues == 0) = [];
VentLevels = max(length(VentValues),max(VentValues));
% if length(VentLevels)<max(VentLevels)
%     VentLevels = 1:max(VentLevels);
% end

if VentLevels >4 && VentLevels<7
    vent_levels = {'Red','DarkOrange','YellowGreen','GreenYellow','SteelBlue','MediumBlue'};
    if is_defect_mask
        vent_levels = {'Red'};
    end
elseif VentLevels ==4
    vent_levels = {'Crimson','Orange','LawnGreen','Blue'};
elseif VentLevels ==3
    vent_levels = {'Crimson','LawnGreen','Blue'};
    
end

for sl=slice_range
    plot_count = plot_count + 1;
    greydata = proton(plotting_range.y,plotting_range.x,sl);%proton(:,:,sl);
    greydata(isnan(greydata))=0;
    slice_gas=gas(plotting_range.y,plotting_range.x,sl);%gas(:,:,sl);
    %     if ~strcmp(type,'CoV')
    minIntensity = min(greydata(:));
    maxIntensity = max(greydata(:));
    if maxIntensity==0
        maxIntensity=1;
    end
    if ~is_ventilation_mask && ~is_defect_mask
        BlackAndWhite = (greydata-minIntensity)./(maxIntensity-minIntensity);
        rgb(:,:,1) = BlackAndWhite;
        rgb(:,:,2) = BlackAndWhite;
        rgb(:,:,3) = BlackAndWhite + slice_gas./max(slice_gas(:));
    else
        BlackAndWhite = (greydata-minIntensity)./(maxIntensity-minIntensity);
        rgb = repmat(BlackAndWhite,[1,1,3]);
        colorVec = zeros(VentLevels,3);
        
        %         if VentLevels ==6
        
        for n_level = 1:VentLevels
            vent_color = vent_levels{n_level};
            vent.(vent_color) = zeros(length(plotting_range.y),length(plotting_range.x));%zeros(nRows,nCols);
            vent.(vent_color)(slice_gas==n_level)=1;
            colorVec(n_level,:) = rgb_color_to_code(vent_color);
            rgb(:,:,1) = rgb(:,:,1) + colorVec(n_level,1).*vent.(vent_color);
            rgb(:,:,2) = rgb(:,:,2) + colorVec(n_level,2).*vent.(vent_color);
            rgb(:,:,3) = rgb(:,:,3) + colorVec(n_level,3).*vent.(vent_color);
        end
        %         elseif VentLevels==4
        %             % defect-red; most ventilated -blue;
        %             for n_level = 1:numel(vent_levels)
        %                 vent_color = vent_levels{n_level};
        %                 colorVec(n_level,:) = rgb_color_to_code(vent_color);
        %                 vent.(vent_color) = zeros(nRows,nCols);
        %                 vent.(vent_color)(slice_gas==n_level)=1;
        %             end
        %             BlackAndWhite = (greydata-minIntensity)./(maxIntensity-minIntensity);
        %             rgb(:,:,1) = BlackAndWhite + vent.red + .7*vent.orange;
        %             rgb(:,:,2) = BlackAndWhite + vent.green;
        %             rgb(:,:,3) = BlackAndWhite + vent.blue + .3*vent.orange;
        %
        %         end
    end
    
    overlay_data(plotting_range.y,plotting_range.x, :, sl) = rgb;
    if fig_col_real>fig_col && plot_count == fig_col_real
        continue;
    end
    if plot_count ==nPlots-1
        sub1 = subaxis(fig_row,fig_col_real,plot_count,'SpacingVert',0,'MR',0,...
            'SpacingHoriz',0,'Spacing',0);
    elseif plot_count == nPlots
        sub2 = subaxis(fig_row,fig_col_real,plot_count,'SpacingVert',0,'MR',0, ...
            'SpacingHoriz',0,'Spacing',0);
        
    else
        subaxis(fig_row,fig_col_real,plot_count,'SpacingVert',0,'MR',0, ...
            'SpacingHoriz',0,'Spacing',0);
%             'SpacingHoriz',0,'Spacing',0,'MarginLeft',0.05,'MarginRight',0.1);
    end
    bkg = repmat(zeros(length(plotting_range.y),length(plotting_range.x)),[1,1,3]);
 imshow(bkg,[]);axis image off;
 switch type
     case 'overlay'
         %     if ~strcmp(type,'CoV')
         
         hold on;
         %         imagesc(rgb(plotting_range.y,plotting_range.x,:));brighten(.3);
         imshow(rgb,[]);%brighten(.3);
         % else
     case 'CoV'
         if sum(slice_gas(:)>0)
             %             hold on;imagesc(slice_gas,[0,0.6]); hold off;
             hold on;
%              imshow(slice_gas,[0,0.60]);
             imshow(slice_gas,[0,0.5]);
         end
         %             imshow(rgb(plotting_range.y,plotting_range.x,:),[0,2]);%axis image off;
         %             h_im = imagesc(slice_gas(plotting_range.y,plotting_range.x,:),[0,2]);
         %             alphadata = ones(size(slice_gas));
         %             alphadata(slice_gas<0)=0;
         %             set(h_im,'AlphaData',alphadata);
         %             cMapName = 'jet2';
         %             cmapSize = 120;
         colors = {'Black','Indigo','MediumBlue','DodgerBlue','Cyan','Lime','GreenYellow','Yellow','Goldenrod','DarkOrange','Red','DarkRed'};
         %         colorVec = jet2(12);
         colorVec = zeros(numel(colors),3);
         for n = 1:numel(colors)
             colorVec(n,:) = rgb_color_to_code(colors{n});
         end
     case 'ADC'
         if sum(slice_gas(:)>0)
             %             hold on;imagesc(slice_gas,[0,0.6]); hold off;
             hold on;
%              imshow(slice_gas,[0,0.60]);
             imshow(slice_gas,[0,0.05]);
         end
         %             imshow(rgb(plotting_range.y,plotting_range.x,:),[0,2]);%axis image off;
         %             h_im = imagesc(slice_gas(plotting_range.y,plotting_range.x,:),[0,2]);
         %             alphadata = ones(size(slice_gas));
         %             alphadata(slice_gas<0)=0;
         %             set(h_im,'AlphaData',alphadata);
         %             cMapName = 'jet2';
         %             cmapSize = 120;
         colors = {'Black','Indigo','MediumBlue','DodgerBlue','Cyan','Lime','GreenYellow','Yellow','Goldenrod','DarkOrange','Red','DarkRed'};
         %         colorVec = jet2(12);
         colorVec = zeros(numel(colors),3);
         for n = 1:numel(colors)
             colorVec(n,:) = rgb_color_to_code(colors{n});
         end
     case 'raw'
         hold on;
         imshow(rgb,[]);colormap('gray');
 end
    %%
    
    if sl == slice_range(round(nCols_subplot/2))
        title(title_string,'FontSize',16,'FontWeight','bold','Color','w')
    elseif plot_count == nPlots && is_ventilation_mask && ~is_defect_mask
        %             if ~strcmp(type,'CoV')
        hc = colorbar('location','eastoutside');
        colormap(colorVec);
        %             hc.Limits = [1,VentLevels+1]; % commented out in 2016a
        hc.Ticks = linspace(hc.Limits(1),hc.Limits(2),VentLevels);
        
        if VentLevels == 4
            %             ventilation_colormap = [1,0,0;1,.6445,0;0,1,0;0,0,1];
            %             hc.TickLabels = {'1 (defect)','2','3','4'};
%             hc.TickLabels = {'VDR','LVR','MVR','HVR'};
             hc.TickLabels = {'VDP','LVP','MVP','HVP'};
       elseif VentLevels == 5
            hc.TickLabels = {'1 (defect)','2','3','4','5'};
            
        elseif VentLevels == 6
            hc.TickLabels = {'1 (defect)','2','3','4','5','6'};
        elseif VentLevels == 3
            hc.TickLabels = {'1 (defect)','2','3'};
        end
        
        hc.Color = 'w';
        hc.FontSize = 20;
        hc.FontWeight = 'bold';
        s1Pos = get(sub1,'position');
        s2Pos = get(sub2,'position');
        s2Pos(3:4) = s1Pos(3:4);
        set(sub2,'position',s2Pos);
        %             end
    elseif plot_count == nPlots && (strcmp(type,'CoV') || strcmp(type,'ADC'))
        hc = colorbar('location','eastoutside');
        colormap(colorVec);
        hc.Ticks = linspace(hc.Limits(1),hc.Limits(2),3);
        hc.Color = 'w';
        hc.FontSize = 30;
        hc.FontWeight = 'bold';
        s1Pos = get(sub1,'position');
        s2Pos = get(sub2,'position');
        s2Pos(3:4) = s1Pos(3:4);
        set(sub2,'position',s2Pos);
    end
    
    %     else
    %         if plot_count == nPlots
    %             sub2 = subaxis(fig_row,fig_col,plot_count,'SpacingVert',0,'MR',0, ...
    %                 'SpacingHoriz',0,'Spacing',0,'MarginLeft',0.05,'MarginRight',0.1);
    %         else
    %                     subaxis(fig_row,fig_col,plot_count,'SpacingVert',0,'MR',0, ...
    %                 'SpacingHoriz',0,'Spacing',0,'MarginLeft',0.05,'MarginRight',0.1);
    %         end
    %                 imshow(slice_gas(plotting_range.y,plotting_range.x,:),[0,2]);colormap('jet')
    %     if plot_count == nPlots
    %                     hc = colorbar('location','eastoutside');
    %                         hc.Color = 'w';
    %             hc.FontSize = 20;
    %             hc.FontWeight = 'bold';
    %
    %             set(sub2,'position',s2Pos);
    %
    %     end
end


