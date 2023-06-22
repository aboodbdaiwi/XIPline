function    VDPHelp(font_size)
            c1='Preforming VDP Xenon Analysis!\n';
            c2='To Preforming Original Analysis:\n';
            c3='-Choose Original from analysis type, then thresholds and general set up.\n';
            c4='-Choose file type (single or multiple).\n';
            c5='-Load data and writeout images.\n';
            c6='-Segment Lung Parenchyma (Use proton images to help you segemnt the lungs accurately) then writeout mask.\n';
            c7='-Start VDP Analysis.\n\n';
            
            cc1='To Preforming the N4 Bias Analysis:\n';
            cc2='-Save the original images in ITK Snap as a compressed NIFTI file using conventional file naming pattern.\n';
            cc3='-Load the maskarray.dcm into ITK Snap, binarise it (turn all the values into 1 or 0), make any edits, and save as compressed NIFTI.\n'; 
            cc4='-Save the two compressed NIFTI files into the cpir-z4 workstation using WinSCP (see login details). Choose the correct study folder .\n';
            cc5='-login: IP Address:  10.8.40.49, User id: cpir-z4, Password: xenon129 \n';
            cc6='-Open PUTTY and login, Type the following commands in order to produce an N4 Bias Field Corrected Image Set:\n';
            cc7='$ ls  -this will open up a list of files/commands that you can copy\n'; 
            cc8='$ cd biasFieldcorrection  -this will access the correct folder\n';
            cc9='$ workon virt3-py3   -activates virtual environment for running deep learning software\n';
            cc10='$ python3 cinci\\_bias\\_correction.py   -to begin running the iterations\n';
            cc11='You should now see a compressed NIFTI file saved in the same folder as the other files.\n'; 
            cc12='Copy the three NIFTI files into your own workstation and save them in a suitable place.\n';            
            cc13='-Back to HPXenon App and choose N4 Bias from analysis type, then thresholds and general set up.\n.\n';
            cc14='Load data (Ventilation Nifti file) and writeout images.\n'; 
            cc15='Load mask (Ventilation Mask Nifti file) and writeout mask.\n'; 
            cc16='-Start VDP Analysis.\n\n';
            
            cc17='All analysis result will be in the data loaction. \n\n';
            cc18='If the program stopped somewhere, use debugging tools and start from beginning.\n\n';
            cc19='To analyze new subject, click "Restart".\n';

            message1=[c1,c2,c3,c4,c5,c6,c7,cc1, cc2, cc3, cc4,cc5,cc6,cc7,cc8,cc9,cc10,cc11,cc12,cc13,cc14,cc15,cc16,cc17,cc18,cc19 ];
            message1 = sprintf(message1);  
%             strrep(sprintf(message1), '_', '\_');
            msgboxw(message1,font_size); 
end