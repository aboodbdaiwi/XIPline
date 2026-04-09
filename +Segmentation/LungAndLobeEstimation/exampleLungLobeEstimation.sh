dataDir=./
outputDir="${dataDir}/Output/"

mkdir -p ${outputDir}

intensityImages=( `ls ${dataDir}/ProtonImages/*H1.nii.gz` );
lungMaskImages=( `ls ${dataDir}/ProtonLungMasks/*Mask.nii.gz` );
lobeMaskImages=( `ls ${dataDir}/ProtonLobeMasks/*Lobes.nii.gz` );

############
#
#    For both the lung and lobe estimation, we need to normalize to the target
#    image.  Simply as a matter of convenience, we're going to make the target
#    image the first image in the arrays extracted above, i.e. ${intensityImages[0]}
#

targetImage=${intensityImages[0]}
targetLungMask=${lungMaskImages[0]}
targetLobeMask=${lobeMaskImages[0]}

targetLungMaskEstimated="${outputDir}/targetLungMaskEstimated.nii.gz"
targetLobeMaskEstimated="${outputDir}/targetLobeMaskEstimated.nii.gz"

warpedImages=()
warpedLungMasks=()
warpedLobeMasks=()

for (( i = 1; i < ${#intensityImages[@]}; i++ ))
  do
    outputPrefix=${outputDir}/proton${i}_

    normalizeCommand="antsRegistrationSyNQuick.sh -d 3 -f $targetImage -m ${intensityImages[$i]} -o $outputPrefix -t b -s 26"
    $normalizeCommand

    # We transform the lung and lobe labels to the space of the target image
    transformCommand="antsApplyTransforms -d 3 -v 1 -i ${lungMaskImages[$i]} -r ${intensityImages[0]} -o ${outputPrefix}LungsWarped.nii.gz -n NearestNeighbor -t ${outputPrefix}1Warp.nii.gz -t ${outputPrefix}0GenericAffine.mat"
    $transformCommand

    transformCommand="antsApplyTransforms -d 3 -v 1 -i ${lobeMaskImages[$i]} -r ${intensityImages[0]} -o ${outputPrefix}LobesWarped.nii.gz -n NearestNeighbor -t ${outputPrefix}1Warp.nii.gz -t ${outputPrefix}0GenericAffine.mat"
    $transformCommand

    warpedImages[${#warpedImages[@]}]="${outputPrefix}Warped.nii.gz"
    warpedLungMasks[${#warpedLungMasks[@]}]="${outputPrefix}LungsWarped.nii.gz"
    warpedLobeMasks[${#warpedLobeMasks[@]}]="${outputPrefix}LobesWarped.nii.gz"
  done

##############
#
#   Estimate target lung mask
#

jointFusionCommand="antsJointFusion -d 3 -v 1 -t $targetImage -o $targetLungMaskEstimated -p 2 -s 3 -a 0.1"
for (( i = 0; i < ${#warpedImages[@]}; i++ ))
  do
    jointFusionCommand="${jointFusionCommand} -g ${warpedImages[$i]}"
    jointFusionCommand="${jointFusionCommand} -l ${warpedLungMasks[$i]}"
  done
$jointFusionCommand


##############
#
#   Estimate target lobe mask
#

jointFusionCommand="antsJointFusion -d 3 -v 1 -t $targetImage -o $targetLobeMaskEstimated -p 2 -s 3 -a 0.1"
for (( i = 0; i < ${#warpedImages[@]}; i++ ))
  do
    jointFusionCommand="${jointFusionCommand} -g ${warpedImages[$i]}"
    jointFusionCommand="${jointFusionCommand} -l ${warpedLobeMasks[$i]}"
  done
$jointFusionCommand


##############
#
#   Perform label overlap evaluation
#


echo
echo
echo
echo

echo "Evaluation"
echo

echo "Lung mask estimation label overlap measures:"
LabelOverlapMeasures 3 $targetLungMask $targetLungMaskEstimated

echo
echo

echo "Lobe mask estimation label overlap measures:"
LabelOverlapMeasures 3 $targetLobeMask $targetLobeMaskEstimated
