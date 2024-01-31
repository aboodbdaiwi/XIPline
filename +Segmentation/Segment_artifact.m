function artifact_mask = Segment_artifact(airway_segment,awcof,diffimg);
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

b1img = diffimg(:,:,:,1);
szimg = size(b1img);
artifact_mask = zeros(szimg);

% Quick and dirty trick for removing airways
if strcmp(airway_segment,'quicknEZ') == 1;
    b5img = diffimg(:,:,:,5);
    bratio = b1img./b5img;
    artifact_mask(bratio>awcof)=1;
    
else % manual airway segmentation

    for slice_num = 1:szimg(3)
   
        artifact = 1;
        airway_slice_mask = zeros(szimg(1),szimg(2));

       while artifact == 1

        imshow(b1img(:,:,slice_num), 'InitialMagnification', 1000, 'DisplayRange', []);
        title(gca,['Image Slice ', int2str(slice_num)])  
        %Construct a quest with three options
        choice = questdlg('Is there an artifact?', ...
        'ROI Choice','Yes','No','No');

        % Handle response
        switch choice
            case 'Yes'
            disp([choice '. Segmenting region containing artifacts.'])
            BW = roipoly;
            roi_mask = ones(szimg(1),szimg(2))-BW;
            airway_slice_mask = airway_slice_mask + BW;
            b1img(:,:,slice_num)= roi_mask.*double(b1img(:,:,slice_num));  

            case 'No'    
            artifact = 0;       
            break
        end 
       end

    close
    
    % generate airway mask for current slice
    airway_slice_mask(airway_slice_mask>0)=1; 
    artifact_mask(:,:,slice_num) = airway_slice_mask;
    
    if slice_num == szimg(3)
        disp('Segmenting artifacts complete.'); 
    end
    
    end
end
end

