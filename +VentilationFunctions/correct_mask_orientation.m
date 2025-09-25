function [CorrectedMask, Ventilation] = correct_mask_orientation(Ventilation)
    
    testmask = Segmentation.SegmentLungthresh(Ventilation.Image,3,1); 
    structuringElementSize = 2;
    se = strel('disk', structuringElementSize); % Create a disk-shaped structuring element
    testmask = double(imclose(testmask, se)); % Perform morphological closing
    structuringElementSize = 2; % Size of the structuring element
    % Create a structuring element for erosion
    se = strel('cube', structuringElementSize); % Create a cube-shaped structuring element
    % Perform erosion to shrink the binary mask
    testmask = imerode(testmask, se); % Perform erosion
    % figure; Global.imslice(testmask,'testmask')
    
    A = Ventilation.LungMask;
    sizeArray1 = size(Ventilation.Image);
    sizeArray2 = size(A);
    B = A;
    % Check if slices in the 3d dimension 
    if ~isequal(sizeArray1, sizeArray2)
        minIndex = find(sizeArray2 == min(sizeArray2));
        if minIndex == 1
            B = permute(B, [2, 3, minIndex]);
        elseif minIndex == 2
            B = permute(B, [1, 3, minIndex]);
        elseif minIndex == 3
            B = permute(B, [1, 2, minIndex]);
        end
    else
        % Display a message indicating that the sizes match
        disp('Array sizes are the same.');
    end
    
    % check if number of mask slices don't match number of vent image slices
    if size(B,3) ~= size(Ventilation.Image,3)
        tempmask = zeros(size(Ventilation.Image));
        sliceIndex = zeros(1,size(testmask, 3));
        for slice = 1:size(testmask, 3)
            % Check if any pixel has a value of 1 in the current slice
            if any(testmask(:, :, slice) == 1, 'all')
                sliceIndex(slice) = 1;
            end
        end
        startslice = find(sliceIndex);
        if length(startslice) == size(B,3)
            tempmask(:,:,startslice(1):startslice(end)) = B;
        else
            try
                tempmask(:,:,startslice(1):startslice(1)+size(B,3)-1) = B;
            catch
                disp('mask has to be fixed manually')
            end
        end
        B = tempmask;
    end
    %figure; Global.imslice(B,'B')
    % check if orientation between mask and vent image is matching
    [DiceScore] = Global.DiceCoeff(testmask, B);
    if DiceScore.DSC < 0.7 % play with this thershold if it's still not working
        DSC = zeros(3,4,3,2);
        for op = 1:3
            if op == 1
                for i = 1:4
                    C = imrotate(B,90*i);
                    [DiceScore] = Global.DiceCoeff(testmask, C);
                    DSC(op,i,1,1) = DiceScore.DSC;
                end
            elseif op == 2
                for i = 1:4
                    C = imrotate(B,90*i);
                    for j = 1:3
                        C = flip(C,j);
                        [DiceScore] = Global.DiceCoeff(testmask, C);
                        DSC(op,i,j,1) = DiceScore.DSC;
                    end
                end
            elseif op == 3
                for i = 1:4
                    C = imrotate(B,90*i);
                    for j = 1:3
                        C = flip(C,j);
                        for k = 1:2
                            C = flip(C,k+1);
                            [DiceScore] = Global.DiceCoeff(testmask, C);
                            DSC(op,i,j,k) = DiceScore.DSC;
                        end
                    end
                end
            end
    
        end
        [~, linearIndex] = max(DSC(:));
        [indx1, indx2, indx3, indx4] = ind2sub(size(DSC), linearIndex);
        if indx1 == 1
            D = imrotate(B,90*indx2); 
        elseif indx1 == 2
            D = imrotate(B,90*indx2); 
            D = flip(D,indx3);
        elseif indx1 == 3
            D = imrotate(B,90*indx2); 
            D = flip(D,indx3);
            D = flip(D,indx4);
        end
    %     figure; Global.imslice(D,'mask')
        B = D;
    end
    CorrectedMask = B;
    Ventilation.LungMask = B;
    Ventilation.AirwayMask = zeros(size(Ventilation.LungMask));

end



