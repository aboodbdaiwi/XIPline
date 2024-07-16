function [snrMap] = recon_snr_map_2Dgr(d,h,fname,cm)
%RECON_SNR_MAP_2Dgr Reconstruct coils SNR map
%                                                            (default)
%         d  Raw (p-file) data  (or pfile fname)
%         h  Header from p-file (or empty)
%     fname  Print <fname>.png and save reco as <fname>.mat ([]) 
%        cm  colormap (min max)                             (0 max(snrMap))
%
%  9/2018 Galen Reed
% 11/2020 Rolf Schulte
if (nargin<1), help(mfilename); return; end

if ~exist('cm','var'), cm = []; end


%% input parameters
if ~exist('fname','var'),    fname = []; end
if ~isempty(fname)
    if ~islogical(fname)
        if ~isempty(regexpi(fname,'\.7$')),   fname = fname(1:end-2); end
        if ~isempty(regexpi(fname,'\.h5$')),  fname = fname(1:end-3); end
        if ~isempty(regexpi(fname,'\.mat$')), fname = fname(1:end-4); end
        fname = [fname '_snr_map2D'];
    else
        warning('islogical(fname); ignoring');
        fname = [];
    end
end


%% reconstruction parameters
params.integrationWindow = 500; % [Hz] spectra integration width for generating image
params.lineBroadening = 1; % [Hz] line broadening filter width 
params.noiseRegionSize = 8; % [pixels] noise calculated from a square with this edge size
params.noiseStdThresh = 5; % threshold for noise masks
params.reconMode = 0; % 0 for multiple images in SNR units, 1 for B1 mapping. 
params.doPlot = 1;% make a plot of the summed spectra with integration limits
%params.noiseBandwidth = 1200; % [Hz] bandwidth of spectra over hich to determine noise
params.plotFontSize = 15;

RECONSNRMAPS = 0;
RECONB1MAP = 1;


%% reading in data, if pfile name is given
if ~isnumeric(d)
    if exist(d,'file')
        [d,h] = read_p(d);
    else
        warning('strange input d/file not existing');
    end
end


%% scaling factor for different FOV/slice thickness
slthick = h.image.slthick;        % slice thickness [mm]
fov = h.rdb_hdr.fov;              % FOV [mm]
xcsi = h.rdb_hdr.xcsi;
ycsi = h.rdb_hdr.ycsi;
zcsi = h.rdb_hdr.zcsi;
if xcsi~=16, warning('xcsi(=%g)~=16',xcsi); end
if ycsi~=16, warning('ycsi(=%g)~=16',ycsi); end
if zcsi>1,   warning('zcsi(=%g)>1',zcsi); end
res = fov./[xcsi ycsi];
scale = 8000/(slthick*prod(res)); % factor relative to (2cm)^3 voxel size   



%% reconstruction
squeezedData = squeeze(d);
ncoils = size(d,6);

% reconstruct individual coil images
[MRSIImages]  = sub_fftAndZeroPad(squeezedData, params, h);

% do a sum of squares over channels if needed
if ncoils>1
    sosImages = sqrt(sum(MRSIImages.*conj(MRSIImages),4)); 
else
    sosImages = MRSIImages;
end

%% MRSI to image
[integratedData totalSpec] = sub_MRSIToImage(sosImages, params, h);


%% turn magnitude images into SNR maps
if(params.reconMode == RECONSNRMAPS)
    [mask, noiseSTD, noiseMEAN] = sub_createMaskAndCalculateNoise(integratedData, params);
    snrMap = (integratedData- noiseMEAN) / noiseSTD;
else
    snrMap = integratedData;
end
snrMap = snrMap*scale;


%% saving
if ~isempty(fname)
    save(fname, 'snrMap','h','params');
end


%% plotting
figure;
if isempty(cm), cm = [0 max(snrMap(:))]; end
imagesc(snrMap, cm);
colormap jet;
colorbar();
set(gca, 'xtick', [], 'ytick', []);

figstr = sprintf('P%05d Exam%d Series%d',...
    h.image.rawrunnum,h.exam.ex_no,h.series.se_no);
set(gcf,'name',figstr);
if ~isempty(fname)
    print(fname,'-dpng','-r600','-painters');
end


end      % main function recon_snr_map.m


%% sub-functions

%% fftAndZeroPad
function [reconImage] = sub_fftAndZeroPad(inputData, params, header)
% fftAndZeroPad
% written by Galen Reed 09/18
% this is a simple routine to reconstruct MRSI data from fidcsi.e
% frequency domain filtering, and kx and ky zero padding are applied prior to FT
%
%
% input data should be parsed as (freq, phase) or (freq, phase, coils)


%
% read headers, define FID filter and t, f, domains
%
dataSize = size(inputData);
nc = 0;
if(length(dataSize) == 2)
    nc = 1;
else
    nc = dataSize(3);
end

nx = round(sqrt(dataSize(1)));
nf = dataSize(2);
ny = nx; % hopefully these are the same
nex = header.image.nex;

disp(['frequency samples = ', num2str(nf)]);
disp(['X samples = ', num2str(nx)]);
disp(['Y samples = ', num2str(ny)]);
disp(['channels = ', num2str(nc)]);
disp(['nex = ', num2str(nex)]);

sweepWidth = header.image.user0;
FIDTimeDomain = linspace(0, nf/sweepWidth, nf);
fidWindow = exp(-FIDTimeDomain * params.lineBroadening);
%
% step 1: parse, compensate for RF chopping
%

apodizedData = zeros([nf, nx, ny, nc]);
summedFIDs = zeros([1 nf]);
summedApodizedFIDs = zeros([1 nf]);

for kk = 1:nc
    for ii = 1:ny
        for jj = 1:nx
            
            ind = jj + (ii-1)*ny;
            
            %NOTE: FIDcsi adds RF chopping when opnex>1
            %an improvement would be to only perform this correction
            % if the header indicates nex>1
            rfChopping = (-1)^ind;
            if(nex>1)
                rfChopping = 1;
            end
            thisFid = squeeze(inputData(ind, :, kk)) * rfChopping;
            
            % window the FID with a filter.
            apodizedFID = thisFid .* fidWindow;
            apodizedData(:,jj,ii,kk) = apodizedFID;
            
            summedFIDs = summedFIDs + thisFid;
            summedApodizedFIDs = summedApodizedFIDs + apodizedFID;
        end
    end
end


%
% step 2: zero pad spatial dimensions (2X)
%
% zeroPaddedData = padarray(apodizedData, [0 nx/2 ny/2 0], 0, 'both');
% padarray requires image_toolbox; removed to reduce ctf file size
zeroPaddedData = truma(apodizedData,false,[size(apodizedData,1) 2*nx 2*ny nc]);

%
% step 3: FFT each dimension (aside from coils)
%
fftData = zeros(size(zeroPaddedData));
for kk = 1:nc
    singleCoilData = squeeze(zeroPaddedData(:, :, :, kk));
    singleCoilFFTData = 1/sqrt(length(singleCoilData(:)))*...
        fftshift(fftn(ifftshift(singleCoilData)));
    % singleCoilFFTData =  fftnc(singleCoilData);
    fftData(:, :, :, kk) = abs(singleCoilFFTData);
end


reconImage = fftData;

end      % sub-function fftAndZeroPad.m

%% MRSIToImage
function [integratedData totalSpec] = sub_MRSIToImage(inputData, params, header)
% written by Galen Reed

maskType = 1;

dataSize = size(inputData);
nf = dataSize(1);
nx = dataSize(2);
ny = dataSize(3);

sweepWidth = header.image.user0;
FIDTimeDomain = linspace(0, nf/sweepWidth, nf);
fidWindow = exp(-FIDTimeDomain * params.lineBroadening);
freqDomain = linspace(-sweepWidth/2, sweepWidth/2, nf);

integratedData = zeros([nx ny]);

sosSummedOverY = sum(inputData, 3);
totalSpec = sum(sosSummedOverY, 2);

if(maskType == 1) % a hard coded window
    integrationPoints = round((params.integrationWindow /sweepWidth)* nf);
    lb = nf/2 - round(integrationPoints/2);
    ub = nf/2 + round(integrationPoints/2);
    integrationBounds = [lb,ub];
    mask = max(totalSpec)*ones(size(totalSpec));
    mask(1:lb-1) = 0;
    mask(ub+1:end) = 0;
else % use peak detection algorithm
    
    dirtyMask = peakDetect(totalSpec);
    
    % clean up the mask a bit.
    structuringElement = strel('line',5,90);
    erodedMask = imerode(dirtyMask, structuringElement);
    dilatedMask = imdilate(erodedMask, structuringElement);
    dilatedMask = imdilate(dilatedMask, structuringElement);
    
    mask = dilatedMask * max(totalSpec);
end
maskIndices = find(mask>0);

for ii = 1:ny
    for jj = 1:nx
        thisSpec = squeeze(inputData(:, jj, ii));
        integratedData(jj, ii) = sum(thisSpec(maskIndices));
    end
end


if (params.doPlot == 1)
    figure();
    plot(freqDomain, totalSpec,...
        freqDomain, mask);
    set(gca, 'fontsize', params.plotFontSize);
    xlabel('freq [Hz]');
    h1 = legend('summed spec', 'integration mask');
    set(h1, 'fontsize', params.plotFontSize);
    set(gca, 'fontsize', params.plotFontSize);
end

end      % sub-function MRSIToImage

%% createMaskAndCalculateNoise
function [mask, noiseSTD, noiseMEAN] = sub_createMaskAndCalculateNoise(inputImage, params)
% creatMaskAndCalculateNoise
% written by Galen Reed 10/26/18


% take a small patch in the top left corner of the image
topLeftCorner = inputImage(1:params.noiseRegionSize, 1:params.noiseRegionSize);
stdCornerNoise = std(topLeftCorner(:));
meanCornerNoise = mean(topLeftCorner(:));


% now use the mean and std of the corner to generate a mask of the phantom
dirtyMask = inputImage>(meanCornerNoise + params.noiseStdThresh*stdCornerNoise);
structuringElement = strel('line',5,90);

% clean up the mask a bit.
erodedMask = imerode(dirtyMask, structuringElement);
dilatedMask = imdilate(erodedMask, structuringElement);

% use pixels outside mask to calculate noise.
% we run the dilation again just to ensure no signal pixels contaminate measurement
doubleDilatedMask = imdilate(dilatedMask, structuringElement);
noisePixels = inputImage(find(doubleDilatedMask == 0));

% output
noiseSTD = std(noisePixels(:));
noiseMEAN = mean(noisePixels(:));
mask = dilatedMask;

if(params.doPlot ==1)
    fs = params.plotFontSize;
    
    figure();
    subplot(2,2,1);
    imagesc(dirtyMask);
    title('dirtyMask');
    set(gca, 'fontsize', params.plotFontSize);
    
    subplot(2,2,2);
    imagesc(erodedMask);
    title('erodedMask');
    set(gca, 'fontsize', params.plotFontSize);
    
    subplot(2,2,3);
    imagesc(dilatedMask);
    title('dilated mask');
    set(gca, 'fontsize', params.plotFontSize);
    
    subplot(2,2,4);
    imagesc(doubleDilatedMask);
    title('doubleDilatedMask');
    set(gca, 'fontsize', params.plotFontSize);
end

end      % sub-function createMaskAndCalculateNoise

