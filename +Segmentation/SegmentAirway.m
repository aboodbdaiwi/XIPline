function airway_mask = SegmentAirway(Image)
%   Inputs:
%      
%   Outputs:
%                   
%   Package: https://github.com/aboodbdaiwi/HP129Xe_Analysis_App
%
%   Author: Abdullah S. Bdaiwi
%   Work email: abdullah.bdaiwi@cchmc.org
%   Personal email: abdaiwi89@gmail.com
%   Website: https://www.cincinnatichildrens.org/research/divisions/c/cpir
%
%   Please add updates at the end. Ex: 3/10/24 - ASB: update .... 

b1img = Image(:,:,:,1);
szimg = size(b1img);
airway_mask = zeros(szimg);
if length(szimg) < 3
    z_size = 1;
else
    z_size = szimg(3);
end
% Airway segmentation
for slice_num = 1:z_size

    airway = 1;
    airway_slice_mask = zeros(szimg(1),szimg(2));

   while airway == 1

        imshow(b1img(:,:,slice_num), 'InitialMagnification', 1000, 'DisplayRange', []);
        title(gca,['Image Slice ', int2str(slice_num)])  
        %Construct a quest with three options
        choice = questdlg('Is there an AIRWAY?', ...
        'ROI Choice','Yes','No','No');

        % Handle response
        switch choice
            case 'Yes'
            disp([choice '. Segmenting region containing airway.'])
            BW = roipoly;
            roi_mask = ones(szimg(1),szimg(2))-BW;
            airway_slice_mask = airway_slice_mask + BW;
            b1img(:,:,slice_num)= roi_mask.*double(b1img(:,:,slice_num));  

            case 'No'    
            airway = 0;       
            break
        end 
   end

    close

    % generate airway mask for current slice
    airway_slice_mask(airway_slice_mask>0)=1; 
    airway_mask(:,:,slice_num) = airway_slice_mask;

    if slice_num == z_size
        disp('Segmenting airways complete.'); 
    end
end

end

