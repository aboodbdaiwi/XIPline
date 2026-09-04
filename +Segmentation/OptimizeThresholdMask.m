function [testmask,bestThreshold,QC] = OptimizeThresholdMask(Image)
Image = double(squeeze(Image));
thresholds = 0.3:0.1:4.0;
nT = numel(thresholds);
Masks = cell(1,nT);
SNRvals = nan(1,nT);
Volumes = nan(1,nT);
Fragmentation = nan(1,nT);
Stability = nan(1,nT);
sx = size(Image,1);
sy = size(Image,2);
cornerSizeX = max(round(sx*0.10),2);
cornerSizeY = max(round(sy*0.10),2);
bg1 = Image(1:cornerSizeX,1:cornerSizeY,:);
bg2 = Image(end-cornerSizeX+1:end,1:cornerSizeY,:);
bg3 = Image(1:cornerSizeX,end-cornerSizeY+1:end,:);
bg4 = Image(end-cornerSizeX+1:end,end-cornerSizeY+1:end,:);
background = [bg1(:); bg2(:); bg3(:); bg4(:)];
noiseSD = std(background,'omitnan');
if noiseSD <= 0 || isnan(noiseSD)
    noiseSD = eps;
end
for ii = 1:nT
    th = thresholds(ii);
    mask = logical(Segmentation.SegmentLungthresh(Image,th,0.6));
    Masks{ii} = mask;
    Volumes(ii) = nnz(mask);
    if Volumes(ii) == 0
        Fragmentation(ii) = 0;
        continue
    end
    lungSignal = mean(Image(mask),'omitnan');
    SNRvals(ii) = lungSignal/noiseSD;
    CC = bwconncomp(mask,26);
    if CC.NumObjects > 0
        componentSizes = cellfun(@numel,CC.PixelIdxList);
        componentSizes = sort(componentSizes,'descend');
        nKeep = min(2,numel(componentSizes));
        Fragmentation(ii) = sum(componentSizes(1:nKeep))/sum(componentSizes);
    else
        Fragmentation(ii) = 0;
    end
end
for ii = 2:nT-1
    M0 = Masks{ii};
    M1 = Masks{ii-1};
    M2 = Masks{ii+1};
    dice1 = 2*nnz(M0 & M1)/max(nnz(M0)+nnz(M1),1);
    dice2 = 2*nnz(M0 & M2)/max(nnz(M0)+nnz(M2),1);
    Stability(ii) = mean([dice1 dice2]);
end
if nT > 2
    Stability(1) = Stability(2);
    Stability(end) = Stability(end-1);
end
validVolumes = Volumes(Volumes > 0);
if isempty(validVolumes)
    warning('No valid threshold mask could be generated.');
    testmask = false(size(Image));
    bestThreshold = NaN;
    QC = struct;
    return
end
medianVolume = median(validVolumes,'omitnan');
validMask = Volumes > 0.5*medianVolume & ...
            Volumes < 1.5*medianVolume & ...
            Fragmentation > 0.85 & ...
            ~isnan(SNRvals);
if ~any(validMask)
    validMask = Volumes > 0 & ~isnan(SNRvals);
end
snrNorm = zeros(1,nT);
validSNR = SNRvals(validMask);
if ~isempty(validSNR)
    snrMin = min(validSNR);
    snrMax = max(validSNR);
    if snrMax > snrMin
        snrNorm(validMask) = ...
            (SNRvals(validMask)-snrMin)/(snrMax-snrMin);
    else
        snrNorm(validMask) = 1;
    end
end
Score = 0.5*snrNorm + ...
        0.3*Stability + ...
        0.2*Fragmentation;
Score(~validMask) = -Inf;
[~,bestIdx] = max(Score);
bestThreshold = thresholds(bestIdx);
testmask = Masks{bestIdx};
QC.Thresholds = thresholds;
QC.SNR = SNRvals;
QC.Volume = Volumes;
QC.Stability = Stability;
QC.Fragmentation = Fragmentation;
QC.Score = Score;
QC.BestIndex = bestIdx;
QC.BestThreshold = bestThreshold;
QC.BestSNR = SNRvals(bestIdx);
QC.BestStability = Stability(bestIdx);
QC.BestFragmentation = Fragmentation(bestIdx);
fprintf('Optimal threshold = %.2f | SNR = %.2f | Stability = %.2f | Fragmentation = %.2f\n', ...
    bestThreshold,QC.BestSNR,QC.BestStability,QC.BestFragmentation);
end