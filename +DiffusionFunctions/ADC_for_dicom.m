

data_type = 'multi_dcm';
if strcmp(data_type, 'single_dcm')==1;
    [parentFile,parentPath] = uigetfile('*.dcm', 'Select Ventilation Dicom Image file');
    cd(parentPath)
    info=dicominfo([parentPath parentFile]);
    A = dicomread(info);
    A = double(squeeze(A));
   diffimg=A;
    
else
    diffimg = DICOM_Load; % now allow for import of extensionless file names and list of .dcm
end
size_diffimg=size(diffimg);
diffimg= reshape(diffimg,[size_diffimg(1), size_diffimg(2),5,4]);
diffimg=permute(diffimg, [1 2 4 3]);















