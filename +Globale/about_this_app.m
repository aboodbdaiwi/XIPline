function    about_this_app(font_size)
            c1='Preforming Diffusion Xenon Analysis!\n\n';
            c2='To begin:\n';
            c3='-Choose the type of filter then load .data file.\n';
            c4='-Write out Images.\n';
            c5='If you alrady have the masks, then click on "Load mask". Four masks have to be included: lung\\_mask,airway\\_mask, final\\_mask,noise\\_mask.\n';
            c6='-If you do not have masks, then segment Lung Parenchyma and airways.\n';
            c7='-Write out final mask.\n';   
            c8='-Calcualte SNR and weighting vector\n';
            c9='-Choose type of fitting then start analysis.\n\n';
            c10='All analysis result will be in the data loaction. \n\n';
            c11='If the program stopped somewhere, use debugging tools and start from beginning.\n\n';
            c12='To analyze new subject, click "Restart".\n';

            message1=[c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12];
            message1 = sprintf(message1);
            Globale.msgboxw(message1,font_size); 
end