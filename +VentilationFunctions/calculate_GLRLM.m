
function [GLRLM,GLRLM_average]  = calculate_GLRLM (Image,quantize,lungmask)
% [GLRLM,GLRLM_average]  = calculate_GLRLM (Image,quantize,lungmask)
% gray level run length matrix computation
%
% input:
% - img, an input grayscale image (RGB images are converted to grayscale)
% - quantize, quantization levels. Normally set to 16. Should be larger than 1.
% - mask, a binary mask to use with values of 1 at the ROI's. 
%
% output: texture features
% -GLRLM: includes (for each slice)
%    1. SHORT RUN EMPHASIS (SRE) 
%    2. LONG RUN EMPHASIS(LRE)
%    3. GRAY LEVEL NON-UNIFORMITY (GLN)
%    4. RUN PERCENTAGE (RP)
%    5. RUN LENGTH NON-UNIFORMITY (RLN)
%    6. LOW GRAY LEVEL RUN EMPHASIS (LGRE)
%    7. HIGH GRAY LEVEL RUN EMPHASIS (HGRE)
% - GLRLM_average
%       avreage of all 7 parameters mentioned above across slices
%
%
% (c) Wout Oude Elferink, 13-5-2015
% University Of Twente, The Netherlands
% modified: Abdullah S. Bdaiwi, 6/9/2021
% university of Cincinnati

%% 
GLRLM = zeros(7,size(Image,3));
for img_id = 1:size(Image,3)
    img = Image(:,:,img_id);
    mask = double(lungmask(:,:,img_id));
    if sum(mask(:)) > 5
        % if color => make gray scale
        if size(img,3)>1
           img = im2gray(img); 
        end
        img = im2double(img); % to double
        % crop the image to the mask bounds for faster processing
        stats = regionprops(mask,'BoundingBox');
        bx = int16(floor(stats.BoundingBox)) + int16(floor(stats.BoundingBox)==0);
        img = img(bx(2):bx(2)+bx(4)-1,bx(1):bx(1)+bx(3)-1);
        mask = mask(bx(2):bx(2)+bx(4)-1,bx(1):bx(1)+bx(3)-1);
        % adjust range
        mini = min(img(:));   % find minimum
        img = img-mini;       % let the range start at 0
        maxi = max(img(:));   % find maximum
        % quantize the image to discrete integer values in the range 1:quantize
        levels = maxi/quantize:maxi/quantize:maxi-maxi/quantize;
        img = imquantize(img,levels);
        % apply the mask
        img(~mask) = 0;
        % initialize glrlm: p(i,j)
        % -  with i the amount of bin values (quantization levels)
        % -  with j the maximum run length (because yet unknown, assume maximum length
        %    of image)
        % -  four different orientations are used (0, 45, 90 and 135 degrees)
        p0 = zeros(quantize,max(size(img)));
        p45 = zeros(quantize,max(size(img)));
        p90 = zeros(quantize,max(size(img)));
        p135 = zeros(quantize,max(size(img)));
        % initialize maximum value for j
        maximgS = max(size(img));
        % add zeros to the borders
        img = padarray(img,[1 1]);
        % initialize rotation
        img45 = imrotate(img,45);
        % find the run length for each quantization level
        for i = 1:quantize
            % find the pixels corresponding to the quantization level
            BW = int8(img == i);
            BWr = int8(img45 == i);    

            % find the start and end points of the run length
            G0e = (BW(2:end-1,2:end-1) - BW(2:end-1,3:end)) == 1;
            G0s = (BW(2:end-1,2:end-1) - BW(2:end-1,1:end-2)) == 1;
            G45e = (BWr(2:end-1,2:end-1) - BWr(2:end-1,3:end)) == 1;
            G45s = (BWr(2:end-1,2:end-1) - BWr(2:end-1,1:end-2)) == 1;
            G90e = (BW(2:end-1,2:end-1) - BW(3:end,2:end-1)) == 1;
            G90s = (BW(2:end-1,2:end-1) - BW(1:end-2,2:end-1)) == 1;
            G135e = (BWr(2:end-1,2:end-1) - BWr(3:end,2:end-1)) == 1;
            G135s = (BWr(2:end-1,2:end-1) - BWr(1:end-2,2:end-1)) == 1;

            % find the indexes
            G0s = G0s'; G0s = find(G0s(:));
            G0e = G0e'; G0e = find(G0e(:));
            G45s = G45s'; G45s = find(G45s(:));
            G45e = G45e'; G45e = find(G45e(:));
            G90s = find(G90s(:));
            G90e = find(G90e(:));
            G135s = find(G135s(:));
            G135e = find(G135e(:));

            % find the lengths
            lengths0 = G0e - G0s + 1;
            lengths45 = G45e - G45s + 1;
            lengths90 = G90e - G90s + 1;
            lengths135 = G135e - G135s + 1;

            % fill the matrix
            p0(i,:) = hist(lengths0,1:maximgS);
            p45(i,:) = hist(lengths45,1:maximgS);
            p90(i,:) = hist(lengths90,1:maximgS);
            p135(i,:) = hist(lengths135,1:maximgS);    
        end
        % add all orientations
        p = p0 + p45 + p90 + p135;
        % calculate the features
        totSum = sum(p(:));
        SRE = sum(sum(p,1) ./ ((1:maximgS).^2)) / totSum;
        LRE = sum(sum(p,1) .* ((1:maximgS).^2)) / totSum;
        RLN = sum(sum(p,1) .^2) / totSum;
        RP = totSum / sum(mask(:));
        GLN = sum(sum(p,2) .^2) / totSum;
        LGRE = sum(sum(p,2) .* ((1:quantize)'.^2)) / totSum;
        HGRE = sum(sum(p,2) .^2) / totSum;

        GLRLM(1,img_id) = SRE;    %SHORT RUN EMPHASIS (SRE) 
        GLRLM(2,img_id) = LRE;    %LONG RUN EMPHASIS(LRE)
        GLRLM(3,img_id) = RLN;    %GRAY LEVEL NON-UNIFORMITY (GLN) 
        GLRLM(4,img_id) = RP;     %RUN PERCENTAGE (RP) 
        GLRLM(5,img_id) = GLN;    %RUN LENGTH NON-UNIFORMITY (RLN)
        GLRLM(6,img_id) = LGRE;   %LOW GRAY LEVEL RUN EMPHASIS (LGRE)
        GLRLM(7,img_id) = HGRE;   %HIGH GRAY LEVEL RUN EMPHASIS (HGRE)
    end 
end

% find the mean across slices
[ii,~,v] = find(GLRLM);
GLRLM_average = accumarray(ii,v,[],@mean);

disp('GLRLM analysis completed')

end
    