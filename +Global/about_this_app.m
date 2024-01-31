function    about_this_app(font_size)
    c1 = 'Hyperpolarized 129Xe MRI Analysis App\n';
    c2 = 'A MATLAB-based user interface application is developed\n';
    c3 = 'to analyze hyperpolarized 129Xe MRI data. \n';
    c4 = 'This open source, research tool provides basic\n';
    c5 = 'analysis of hyperpolarized 129Xe calibration, \n';
    c6 = 'ventilation, diffusion, and gas exchange \n';
    c7 = 'images in a single, customizable platform.\n';
    c8 = '\n';   
    c9 = 'For any questions, please email Abdullah Bdaiwi\n';
    c10 = 'abdullah.bdaiwi@cchmc.org\n';

    %message1=[c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12];
    message1=[c1,c2,c3,c4,c5,c6,c7,c8,c9, c10];
    message1 = sprintf(message1);
    Global.msgboxw(message1,font_size); 
end