# LungAndLobeEstimationExample

Example code and data for the techniques discussed in
Tustison et al., [Atlas-based estimation of lung and lobar anatomy in proton MRI. _Magn Reson Med_](http://www.ncbi.nlm.nih.gov/pubmed/26222827).

## To do:

* Need to run it through twice as discussed in the paper.  The first run is to simply extract the lungs which gives a better registration for the second round.

## Instructions

* Download and install [ANTs](https://github.com/stnava/ANTs).
* Add the ANTs binary directory to your path environment variable.
* Run ``exampleLungLobeEstimation.sh`` from the ``LungandLobeEstimationExample/`` directory.

## Results

```
Evaluation

Lung mask estimation label overlap measures
                                          ************ All Labels *************
                      Total  Union (jaccard)      Mean (dice)      Volume sim.   False negative   False positive
                   0.978303         0.956787         0.977916      0.000791194        0.0216966        0.0224703
                                       ************ Individual Labels *************
     Label           Target  Union (jaccard)      Mean (dice)      Volume sim.   False negative   False positive
         1         0.975608         0.959328         0.979242      -0.00745017        0.0243923        0.0170967
         2         0.981569         0.953745         0.976325        0.0106853        0.0184308        0.0288635


Lobe mask estimation label overlap measures
                                          ************ All Labels *************
                      Total  Union (jaccard)      Mean (dice)      Volume sim.   False negative   False positive
                   0.956222         0.911176         0.953524       0.00564215        0.0437785        0.0491584
                                       ************ Individual Labels *************
     Label           Target  Union (jaccard)      Mean (dice)      Volume sim.   False negative   False positive
         1         0.963516         0.919172         0.957884        0.0116906        0.0364842        0.0476829
         2         0.970837         0.927132         0.962188        0.0178171        0.0291631        0.0463078
         3         0.919346         0.900792         0.947807        -0.061917        0.0806545        0.0219128
         4         0.966887         0.855415         0.922074        0.0926963        0.0331126         0.118769
         5           0.9711         0.917791         0.957133        0.0287653        0.0288997        0.0564376
```

