function pfile2mat(fname)
%DICOM2MAT Read in dicom file and save as matlab .mat

[d,h] = read_p(fname);
save([fname '.mat'],'-mat','-v7.3','d','h');
