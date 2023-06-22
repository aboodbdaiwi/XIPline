function my_pixinfo_tool(im)
% Create figure, setting up properties
fig = figure('Toolbar','none', ...
              'Menubar','none', ...
              'Name','My Pixel Info Tool', ...
              'NumberTitle','off', ...
              'IntegerHandle','off','position',[500 200 600 800]);
% im = im./max(im(:));
% Create axes and reposition the axes
% to accommodate the Pixel Region tool panel
ax = axes('Units','normalized', ...
           'Position',[0 .5 1 .5]);

% Display image in the axes
img = imshow(im,[min(im(:)) max(im(:))]); 

if length(size(im))== 2
    colorbar; caxis([0 max(im(:))])
end
% Add Distance tool, specifying axes as parent
distool = imdistline(ax);

% Add Pixel Information tool, specifying image as parent
pixinfo = impixelinfo(img);

% Add Display Range tool, specifying image as parent
drange = imdisplayrange(img);

% Add Pixel Region tool panel, specifying figure as parent
% and image as target
pixreg = impixelregionpanel(fig,img);

% Reposition the Pixel Region tool to fit in the figure
% window, leaving room for the Pixel Information and
% Display Range tools
set(pixreg, 'units','normalized','position',[0 .08 1 .4])

end