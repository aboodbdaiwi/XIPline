function CoSeg3D(filename_H,filename_CT,filename_He,filename_CT_H_def,filename_CT_He_def,filename_Def1_Ux,filename_Def1_Uy,filename_Def1_Uz,filename_Def2_Ux,filename_Def2_Uy,filename_Def2_Uz)

g1 = metaImageRead(filename_H);
g1 = double(g1);
img1_orig = g1 ;
g1 = (g1 - min(g1(:)))/(max(g1(:)) - min(g1(:)))*120;%150


g2 = metaImageRead(filename_CT);
g2_hdr = metaImageInfo(filename_CT);
g2 = double(g2);
img2_orig = g2 ;
g2 = (g2 - min(g2(:)))/(max(g2(:)) - min(g2(:)))*90;%150

g3 = metaImageRead(filename_He);
g3 = double(g3);
img3_orig = g3 ;
g3 = (g3 - min(g3(:)))/(max(g3(:)) - min(g3(:)))*80;%80

[rows, cols, heights] = size(g1);

numIter = 200;
steps = 0.10; 
alpha1 = 5;  % 5, weight of R(|\delta(u1(x))|)
alpha2 = 8;  % 8, weight of R(|\delta(u2(x))|)
gamma = 0.5; % weight of E(u1-u2)
weight1 = 0.5; %0.6 weight of E(CT-H)
weight2 = 0.8; %0.8 weight of E(CT-He)

cc = 0.35;

U1x = zeros(2,2,2,'single'); 
U1y = U1x; 
U1z = U1x;

U2x = U1x;
U2y = U1x; 
U2z = U1x;

wi = 0.8;

errbound = 2e-4;
warps = [4 4 2];
levels = [4 2 1];
sigma = 0.7;
offset = 1;

tic
for j = 1:length(levels)
       
    hs = fspecial('gaussian',[10,1],levels(j));
    
    g1_l = volresize(volfilter(g1,hs),ceil(size(g1)./levels(j)));
    g2_l = volresize(volfilter(g2,hs),ceil(size(g2)./levels(j)));
    g3_l = volresize(volfilter(g3,hs),ceil(size(g3)./levels(j)));
    
    % resize volumes for current level
    
    [U1x,U1y,U1z] = resizeFlow(U1x,U1y,U1z,size(g1_l));
    [U2x,U2y,U2z] = resizeFlow(U2x,U2y,U2z,size(g1_l));
    
    
    u1x = zeros(ceil(size(g1)./levels(j)));
    u1y = zeros(ceil(size(g1)./levels(j)));
    u1z = zeros(ceil(size(g1)./levels(j)));
    
    u2x = zeros(ceil(size(g1)./levels(j)));
    u2y = zeros(ceil(size(g1)./levels(j)));
    u2z = zeros(ceil(size(g1)./levels(j)));
    
    [n_rows,n_cols,n_heights]=size(g1_l);
    
    for k = 1:warps(j)
        
        U1x = U1x + u1x*wi;
        U1y = U1y + u1y*wi;
        U1z = U1z + u1z*wi;     
        
        U2x = U2x + u2x*wi;
        U2y = U2y + u2y*wi;
        U2z = U2z + u2z*wi;  
        
        
        g1_w = g1_l;
        g2_w1 = volWarp(g2_l,U1x,U1y,U1z);
        g2_w2 = volWarp(g2_l,U2x,U2y,U2z);
        g3_w = g3_l;
        
        [~, ssc_q1_w] = SSC_descriptor_H(g1_w, sigma, offset);
        [~, ssc_q2_w1] = SSC_descriptor_CT(g2_w1, sigma, offset);
        [~, ssc_q2_w2] = SSC_descriptor_CT(g2_w2, sigma, offset);
        [~, ssc_q3_w] = SSC_descriptor_He(g3_w, sigma, offset);
        
    
        % switch between SAD and SSC calculation of image similarity.
        g1t = hammingDist3D(ssc_q2_w1, ssc_q1_w);
        g2t = hammingDist3D(ssc_q2_w2, ssc_q3_w);

        g1x = (hammingDist3D(volshift(ssc_q2_w1,1,0,0), ssc_q1_w) - hammingDist3D(volshift(ssc_q2_w1,-1,0,0), ssc_q1_w))/2;
        g1y = (hammingDist3D(volshift(ssc_q2_w1,0,1,0), ssc_q1_w) - hammingDist3D(volshift(ssc_q2_w1,0,-1,0), ssc_q1_w))/2;
        g1z = (hammingDist3D(volshift(ssc_q2_w1,0,0,1), ssc_q1_w) - hammingDist3D(volshift(ssc_q2_w1,0,0,-1), ssc_q1_w))/2;
        
        g2x = (hammingDist3D(volshift(ssc_q2_w2,1,0,0), ssc_q3_w) - hammingDist3D(volshift(ssc_q2_w2,-1,0,0), ssc_q3_w))/2;
        g2y = (hammingDist3D(volshift(ssc_q2_w2,0,1,0), ssc_q3_w) - hammingDist3D(volshift(ssc_q2_w2,0,-1,0), ssc_q3_w))/2;
        g2z = (hammingDist3D(volshift(ssc_q2_w2,0,0,1), ssc_q3_w) - hammingDist3D(volshift(ssc_q2_w2,0,0,-1), ssc_q3_w))/2;
        
               
        g1x(:,[1,n_cols],:) = 0;
        g1y([1,n_rows],:,:) = 0;
        g1z(:,:,[1,n_heights]) = 0;
        
        g2x(:,[1,n_cols],:) = 0;
        g2y([1,n_rows],:,:) = 0;
        g2z(:,:,[1,n_heights]) = 0;
        
        g1f = g1x.*g1x + g1y.*g1y + g1z.*g1z;
        g2f = g2x.*g2x + g2y.*g2y + g2z.*g2z;
        
        % - para: a sequence of parameters for the algorithm
        %   para[0,1,2]: rows, cols, heights of the given image
        %   para[3]: the maximum iteration number
        %   para[4]: the error bound for convergence
        %   para[5]: cc for the step-size of augmented Lagrangian method
        %   para[6]: the step-size for the graident-projection step to the
        %   total-variation function.
        varaParas = [n_rows; n_cols; n_heights; levels(j)*numIter; errbound; cc; steps; alpha1; alpha2; gamma; weight1; weight2];
        
        % CPU based flow adpating
%         [u1x, u1y, u1z, u2x, u2y, u2z, err, num, timet] = CoReg_TVL1_Newton_mex(single(varaParas), single(U1x), single(U1y), single(U1z), single(g1x), single(g1y), single(g1z), single(g1t), single(g1f), ...
%                                                                                                        single(U2x), single(U2y), single(U2z), single(g2x), single(g2y), single(g2z), single(g2t), single(g2f)); 
        
        % GPU based flow adpating
        [u1x, u1y, u1z, u2x, u2y, u2z, err, num, timet] = CoReg_TVL1_Newton_GPU_Ind_opt(single(varaParas), single(U1x), single(U1y), single(U1z), single(g1x), single(g1y), single(g1z), single(g1t), single(g1f), ...
                                                                                                       single(U2x), single(U2y), single(U2z), single(g2x), single(g2y), single(g2z), single(g2t), single(g2f));  
        
    end
    %     % testing
    %     out_vol = volWarp(g2q,Ux,Uy,Uz);
    %     slice=uint16(size(out_vol,3)/2);
    %     figure();
    %     subplot(2,3,1); imshow(g1q(:,:,slice),[]); title('g1q')
    %     subplot(2,3,2); imshow(g2q(:,:,slice),[]); title('g2q')
    %     subplot(2,3,3); imshow(out_vol(:,:,slice),[]); title('g2q_w')
    %     subplot(2,3,4); imshow(abs(g1q(:,:,slice)-out_vol(:,:,slice)),[]); title('abs(g1q-g2q_w)')
    %     subplot(2,3,5); quiver(Uy(:,:,slice),Ux(:,:,slice)); axis('equal'); title('Uy,Ux')
    %     tmp = zeros([size(g1q(:,:,slice)),3]);
    %     tmp(:,:,1) = mat2gray(Uy(:,:,slice)); tmp(:,:,2) = mat2gray(Ux(:,:,slice)); tmp(:,:,3) = mat2gray(Uz(:,:,slice));
    %     subplot(2,3,6); imshow(tmp);

    U1x = U1x + u1x*wi;
    U1y = U1y + u1y*wi;
    U1z = U1z + u1z*wi;
    
    U2x = U2x + u2x*wi;
    U2y = U2y + u2y*wi;
    U2z = U2z + u2z*wi;
end

toc
timet = toc

if(levels(end)>1)
    [U1x,U1y,U1z] = resizeFlow(U1x,U1y,U1z, size(g1));
    [U2x,U2y,U2z] = resizeFlow(U2x,U2y,U2z, size(g1));
end

% record the deformation for future reference.
% metaImageWrite(U1x, filename_Def1_Ux, g2_hdr);
% metaImageWrite(U1y, filename_Def1_Uy, g2_hdr);
% metaImageWrite(U1z, filename_Def1_Uz, g2_hdr);

metaImageWrite(U2x, filename_Def2_Ux, g2_hdr);
metaImageWrite(U2y, filename_Def2_Uy, g2_hdr);
metaImageWrite(U2z, filename_Def2_Uz, g2_hdr);

% g2_def1 = volWarp(g2,U1x,U1y,U1z);
% g2_def1 = volWarp(img2_orig,U1x,U1y,U1z);
% metaImageWrite(g2_def1, filename_CT_H_def, g2_hdr);

g2_def2 = volWarp(g2,U2x,U2y,U2z);
g2_def2 = volWarp(img2_orig,U2x,U2y,U2z);
metaImageWrite(g2_def2, filename_CT_He_def, g2_hdr);
end
