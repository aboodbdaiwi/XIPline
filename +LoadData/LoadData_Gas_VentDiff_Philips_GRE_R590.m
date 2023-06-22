
function [Image, parentPath, FileNames] = LoadData_Gas_VentDiff_Philips_GRE_R590(DataLocation)

    DataFiles = dir([DataLocation,'\*.data']);
    FileNames = DataFiles.name;
    parentPath = DataFiles.folder;
    ch_range = [1 1]; % assume one channel for data import by default
    filterim = 'Fermi';
    cd(parentPath)
    filename = FileNames;
    filename(end-4:end)=[];
    [ImgK,NoiK,kx_oversample_factor] = LoadData.load_philips_extr1_2D(filename,ch_range);
    diffK = ImgK(:,:,:,:);
    
    
    if  length(size(diffK)) == 4
        data_size = size(diffK);
        ky = data_size(1);
        kx =data_size(2);
        slices = data_size(3);
        nbs = data_size(4);
        diffimgcmplx = zeros(ky,kx/kx_oversample_factor,slices, nbs);
        if strcmp(filterim,'Fermi') == 1
                filter= LoadData.fermi_filter_2D_AB(ky,kx,10);
        else
                filter = 1;
        end
        kx_size = kx/kx_oversample_factor;
%         Image = zeros(128,128,slices,nbs);
        for sl =1:slices
            for b = 1:nbs
                data_slice = diffK(:,:,sl,b).*filter;
                recon_slice = fftshift(fft2(data_slice),2);
                recon_slice(:,1:(kx/2)-(kx_size/2))=[]; % crop extra pixels from oversample factor
                recon_slice(:, (1+kx_size):end)=[]; 
                diffimgcmplx(:,:,sl,b) = recon_slice;
                diffimg(:,:,sl,b) = abs(diffimgcmplx(:,:,sl,b)); % make magnitude data
                diffimg = rot90(rot90(diffimg));
%                 Image(:,:,sl,b) = imresize(diffimg(:,:,sl,b),[128 , 128]);
                Image(:,:,sl,b) = diffimg(:,:,sl,b) ;
            end
        end

    elseif length(size(diffK)) == 3
            data_size = size(diffK);
            ky = data_size(1);
            kx =data_size(2);
            slices = data_size(3);
            diffimgcmplx = zeros(ky,round(kx/kx_oversample_factor),slices);
            if strcmp(filterim,'Fermi') == 1
                    filter= LoadData.fermi_filter_2D_AB(ky,kx,10);
            else
                    filter = 1;
            end
            kx_size = kx/kx_oversample_factor;
%             Image = zeros(128,128,slices);
            for sl =1:slices
                    data_slice = diffK(:,:,sl).*filter;
                    recon_slice = fftshift(fft2(data_slice),2);
                    recon_slice(:,1:(kx/2)-(kx_size/2))=[]; % crop extra pixels from oversample factor
                    recon_slice(:, (1+kx_size):end)=[]; 
                    diffimgcmplx(:,:,sl) = recon_slice;
                    diffimg(:,:,sl) = abs(diffimgcmplx(:,:,sl)); % make magnitude data
                    diffimg = rot90(rot90(diffimg));
%                     Image(:,:,sl) = imresize(diffimg(:,:,sl),[128 , 128]);
                    Image(:,:,sl) = diffimg(:,:,sl) ;
            end
            
    end
    
end