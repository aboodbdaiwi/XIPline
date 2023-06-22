function lung_mask = SegmentLungOstuthresh(Image,nzcof)

[counts,~] = imhist(Image,16);
T = otsuthresh(counts);
lung_mask = imbinarize(Image,T.*nzcof);


% lung_mask = imbinarize(NVentImage,'adaptive','ForegroundPolarity','bright','Sensitivity',0.4);
%figure; imslice(lungmask)
end