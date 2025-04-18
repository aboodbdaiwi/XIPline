
function [Image] = cartesian_2D_recon(MainInput)
    %% Working with an existing ISMRMRD data set
    
    % This is a simple example of how to reconstruct images from data
    % acquired on a fully sampled cartesian grid
    %
    % Capabilities:
    %   2D/3D
    %   multiple slices/slabs
    %   multiple contrasts, repetitions
    %   
    % Limitations:
    %   only works with a single encoded space
    %   fully sampled k-space (no partial fourier or undersampling)
    %   multiple repetitions
    %   doesn't handle averages, phases, segments and sets
    %   ignores noise scans (no pre-whitening)
    % 
    
    % We first create a data set using the example program like this:
    %   ismrmrd_generate_cartesian_shepp_logan -r 5 -C -o shepp-logan.h5
    % This will produce the file shepp-logan.h5 containing an ISMRMRD
    % dataset sampled evenly on the k-space grid -128:0.5:127.5 x -128:127
    % (i.e. oversampled by a factor of 2 in the readout direction)
    % with 8 coils, 5 repetitions and a noise level of 0.5
    % with a noise calibration scan at the beginning
    %
    %
    clc
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    % Loading an existing file %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% 

    cd(MainInput.XeDataLocation)

    if strcmp(MainInput.ReconImageMode,'xenon')
        if exist(MainInput.XeFileName, 'file')
            dset = LoadData.ismrmrd.Dataset(MainInput.XeFileName, 'dataset');
        else
            error(['File ' MainInput.XeFileName ' does not exist.  Please generate it.'])
        end
    elseif strcmp(MainInput.ReconImageMode,'proton')
        if exist(MainInput.HFileName, 'file')
            dset = LoadData.ismrmrd.Dataset(MainInput.HFileName, 'dataset');
        else
            error(['File ' MainInput.HFileName ' does not exist.  Please generate it.'])
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Read some fields from the XML header %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % We need to check if optional fields exists before trying to read them
    
    hdr = LoadData.ismrmrd.xml.deserialize(dset.readxml);
    
    %% Encoding and reconstruction information
    % Matrix size
    enc_Nx = hdr.encoding.encodedSpace.matrixSize.x;
    enc_Ny = hdr.encoding.encodedSpace.matrixSize.y;
    enc_Nz = hdr.encoding.encodedSpace.matrixSize.z;
    rec_Nx = hdr.encoding.reconSpace.matrixSize.x;
    rec_Ny = hdr.encoding.reconSpace.matrixSize.y;
    rec_Nz = hdr.encoding.reconSpace.matrixSize.z;
    
    % Field of View
    enc_FOVx = hdr.encoding.encodedSpace.fieldOfView_mm.x;
    enc_FOVy = hdr.encoding.encodedSpace.fieldOfView_mm.y;
    enc_FOVz = hdr.encoding.encodedSpace.fieldOfView_mm.z;
    rec_FOVx = hdr.encoding.reconSpace.fieldOfView_mm.x;
    rec_FOVy = hdr.encoding.reconSpace.fieldOfView_mm.y;
    rec_FOVz = hdr.encoding.reconSpace.fieldOfView_mm.z;
    
    % Number of slices, coils, repetitions, contrasts etc.
    % We have to wrap the following in a try/catch because a valid xml header may
    % not have an entry for some of the parameters
    
    try
      nSlices = hdr.encoding.encodingLimits.slice.maximum + 1;
    catch
        nSlices = 1;
    end
    
    try 
        nCoils = hdr.acquisitionSystemInformation.receiverChannels;
    catch
        nCoils = 1;
    end
    
    try
        nReps = hdr.encoding.encodingLimits.repetition.maximum + 1;
    catch
        nReps = 1;
    end
    
    try
        nContrasts = hdr.encoding.encodingLimits.contrast.maximum + 1;
    catch
        nContrasts = 1;
    end
    
    % TODO add the other possibilities
    
    %% Read all the data
    clc
    % Reading can be done one acquisition (or chunk) at a time, 
    % but this is much faster for data sets that fit into RAM.
    D = dset.readAcquisition();
    
    % Note: can select a single acquisition or header from the block, e.g.
    % acq = D.select(5);
    % hdr = D.head.select(5);
    % or you can work with them all together
    
    %% Ignore noise scans
    % TODO add a pre-whitening example
    % Find the first non-noise scan
    % This is how to check if a flag is set in the acquisition header
    isNoise = D.head.flagIsSet('ACQ_IS_NOISE_MEASUREMENT');
    firstScan = find(isNoise==0,1,'first');
    if firstScan > 1
        noise = D.select(1:firstScan-1);
    else
        noise = [];
    end
    meas  = D.select(firstScan:D.getNumber);
    clear D;
    
    %% Reconstruct images
    % Since the entire file is in memory we can use random access
    % Loop over repetitions, contrasts, slices
    reconImages = {};
    nimages = 0;
    for rep = 1:nReps
        for contrast = 1:nContrasts
            for slice = 1:nSlices
                % Initialize the K-space storage array
                K = zeros(enc_Nx, enc_Ny, enc_Nz, nCoils);
                % Select the appropriate measurements from the data
                acqs = find(  (meas.head.idx.contrast==(contrast-1)) ...
                            & (meas.head.idx.repetition==(rep-1)) ...
                            & (meas.head.idx.slice==(slice-1)));
                for p = 1:length(acqs)
                    ky = meas.head.idx.kspace_encode_step_1(acqs(p)) + 1;
                    kz = meas.head.idx.kspace_encode_step_2(acqs(p)) + 1;
                    K(:,ky,kz,:) = meas.data{acqs(p)};
                end
                % Reconstruct in x
                K = fftshift(ifft(fftshift(K,1),[],1),1);
                % Chop if needed
                if strcmp(MainInput.Scanner,'Siemens')
                    if strcmp(MainInput.AnalysisType,'Ventilation') 
    %                      centerIndex = size(K, 1) / 2;                    
    %                     ind1 = round(centerIndex - enc_Ny/2) + 1;
    %                     ind2 = round(centerIndex + enc_Ny/2);                    
    %                     im = K(ind1:ind2,:);
                        if strcmp(MainInput.ReconImageMode,'xenon')
                            crop_factor = round(enc_Ny*0.75);
                        elseif strcmp(MainInput.ReconImageMode,'proton')
                            crop_factor = round(enc_Ny*0.5);
                        end
                        ind1 = crop_factor + 1;
                        ind2 = crop_factor + enc_Ny;                    
                        im = K(ind1:ind2,:,:,:);
                        im = flip(im,1);
                        im = flip(im,2);  
                        %disp('true 2')
                    elseif strcmp(MainInput.AnalysisType,'Diffusion') 
                        crop_factor = round(enc_Ny*1.3);
                        ind1 = crop_factor + 1;
                        ind2 = crop_factor + enc_Ny;                    
                        im = K(ind1:ind2,:,:,:);
                        im = flip(im,1);
                        im = flip(im,2);  
                        %disp('true 3')
                    end
                elseif (enc_Nx == rec_Nx)
                    im = K;
                else
                    ind1 = floor((enc_Nx - rec_Nx)/2)+1;
                    ind2 = floor((enc_Nx - rec_Nx)/2)+rec_Nx;
                    im = K(ind1:ind2,:,:,:);
                    %disp('true 4')
                end

                % Reconstruct in y then z
                im = fftshift(ifft(fftshift(im,2),[],2),2);
                if size(im,3)>1
                    im = fftshift(ifft(fftshift(im,3),[],3),3);
                end
                
                % Combine SOS across coils
                im = sqrt(sum(abs(im).^2,4));
                
                % Append
                nimages = nimages + 1;
                reconImages{nimages} = im;
            end
        end
    end
    
    %% Display the first image
%     figure
%     colormap gray
%     imagesc(reconImages{2}); axis image; axis off; colorbar;
%     Img1 = double(reconImages{1,2});
%     Image = zeros(size(Img1,1),size(Img1,2),size(reconImages,1),size(reconImages,2),size(reconImages,3));
    Image = zeros(rec_Nx,rec_Nx,size(reconImages,1),size(reconImages,2),size(reconImages,3));
    for i = 1:size(reconImages,1)
        for j = 1:size(reconImages,2)
            for k = 1:size(reconImages,3)
%                 Image(:,:,i,j,k) = double(reconImages{i,j,k});
                Img = double(reconImages{i,j,k});
                Image(:,:,i,j,k) = imresize(Img, [rec_Nx,rec_Nx]);
            end
        end
    end

    if strcmp(MainInput.Scanner,'Siemens')
        if nContrasts > 1
            Image = squeeze(Image);
            Image = reshape(Image, [rec_Nx,rec_Nx,nSlices,nContrasts]);
        else
            Image = squeeze(Image);
            if length(size(Image)) == 3
                Image = squeeze(Image(:,:,1:nSlices));
            elseif length(size(Image)) == 4
                Image = squeeze(Image(:,:,1:nSlices,:));
            end
        end
    end
    %figure; imslice(Image)

end