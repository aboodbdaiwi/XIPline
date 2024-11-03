function CO_segmentation(filename_H_3,filename_H_seeds,filename_He_3,filename_He_seeds,dir,flag)

h = fspecial('gaussian',5,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H_seeds = metaImageRead(filename_H_seeds);
H_img3D = double(metaImageRead(filename_H_3));
[rows,cols,heights] = size(H_img3D);
H_hdr = metaImageInfo(filename_H_3);


He_seeds = metaImageRead(filename_He_seeds);
He_img3D = double(metaImageRead(filename_He_3));
He_img3D(isnan(He_img3D)) = 0;
He_hdr = metaImageInfo(filename_He_3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

min_gray = min(H_img3D(:));
max_gray = max(H_img3D(:));
H_img3D_norm = (H_img3D-min_gray)/(max_gray-min_gray);

min_gray = min(He_img3D(:));
max_gray = max(He_img3D(:));
He_img3D_norm = (He_img3D-min_gray)/(max_gray-min_gray);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ur1_smooth = imfilter(H_img3D_norm,h,'same');
[dx,dy,dz] = gradient(255*ur1_smooth);
grad = sqrt(dx.^2+dy.^2+dz.^2);
max_grad = max(grad(:));
min_grad = min(grad(:));
grad = (grad-min_grad)./(max_grad-min_grad);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% comment the below two lines to switch between co- (ture) and single
% modality (false) segmentation. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% flag = true;       %default is co-segmentation
% flag = false;    %default is co-segmentation
flag
if (flag)   % for co-segmentation
    alpha1 = 0.5 + 3.0*exp(-80*grad);
    alpha1 = alpha1;
    beta = 0.03
else        % for H-segmentation
    alpha1 = 0.06 + 2*exp(-40*grad);
    alpha1 = alpha1/0.4;
    beta = 0.0
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ur2_smooth = imfilter(He_img3D_norm,h,'same');
[dx,dy,dz] = gradient(255*ur2_smooth);
grad = sqrt(dx.^2+dy.^2+dz.^2);
max_grad = max(grad(:));
min_grad = min(grad(:));
grad = (grad-min_grad)./(max_grad-min_grad);

alpha2 = 0.2 + 2.0*exp(-50*grad);
alpha2 = alpha2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
big_num = 1e5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

samples = 0:1:255;
Width = 3;
constant1 = 0.02;
constant2 = 0.01;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IND = round(255*H_img3D_norm)+1;
Ivals = H_img3D_norm(find(H_seeds == 1));
Ivals = Ivals*255;
pdf_lung = ksdensity(Ivals,samples,'Width',Width);

look_lung = -constant1*log(pdf_lung./length(Ivals)+eps);

Ct1 = look_lung(IND);
Ct1(H_seeds == 2)= big_num;
Ct1(:,:,[1,end]) = big_num;
Ct1(:,[1,end],:) = big_num;
Ct1([1,end],:,:) = big_num;

Ivals = H_img3D_norm(find(H_seeds == 2));
Ivals = Ivals*255;
pdf_lung = ksdensity(Ivals,samples,'Width',Width);

look_lung = -constant1*log(pdf_lung./length(Ivals)+eps);
Cs1 = look_lung(IND);
Cs1(H_seeds == 1) = big_num;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IND = round(255*He_img3D_norm)+1;
Ivals = He_img3D_norm(find(He_seeds == 1));
Ivals = Ivals*255;
pdf_lung = ksdensity(Ivals,samples,'Width',Width);

look_lung = -constant2*log(pdf_lung./length(Ivals)+eps);

Ct2 = look_lung(IND);
Ct2(He_seeds == 2) = big_num;
Ct2(:,:,[1,end]) = big_num;
Ct2(:,[1,end],:) = big_num;
Ct2([1,end],:,:) = big_num;

Ivals = He_img3D_norm(find(He_seeds == 2));
Ivals = Ivals*255;
pdf_lung = ksdensity(Ivals,samples,'Width',Width);

look_lung = -constant2*log(pdf_lung./length(Ivals)+eps);
Cs2 = look_lung(IND);
Cs2(He_seeds == 1) = big_num;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

varParas = [rows;cols;heights;300; 5e-5; 0.35; 0.10;beta];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[uu1, uu2, erriter,num,tt] = CoSeg_GPU(single(alpha1), single(alpha2), single(Cs1), ...
                  single(Cs2), single(Ct1), single(Ct2), single(varParas));
uu1 = uu1>0.5;
uu2 = uu2>0.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (flag)
    metaImageWrite(uu1,[dir '/H_out_co.mhd'],H_hdr);
else
    metaImageWrite(uu1,[dir '/H_out_s.mhd'],H_hdr);
end  
end
