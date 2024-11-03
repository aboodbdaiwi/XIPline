function DeformLabels3D(fileNameOriginalLabel,fileNameDefUx,fileNameDefUy,fileNameDefUz,fileNameDeformedLabel)

    OriginalLabel = metaImageRead(fileNameOriginalLabel);
    hdr_OriginalLabel = metaImageInfo(fileNameOriginalLabel);
    DefUx = metaImageRead(fileNameDefUx);
    DefUy = metaImageRead(fileNameDefUy);
    DefUz = metaImageRead(fileNameDefUz);
    
    DeformedLabel = volWarp(OriginalLabel,DefUx,DefUy,DefUz,'nearest');
    
    metaImageWrite(DeformedLabel, fileNameDeformedLabel, hdr_OriginalLabel);
end