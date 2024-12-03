![XIPlineMainLogo_noBackground](https://github.com/aboodbdaiwi/HP129Xe_Analysis_App/assets/36932337/ee751c64-065b-4d8c-94b3-5edc89e14ac1)

# 129Xenon Image Processing Pipeline (XIPline)

XIPline application is developed at the [CPIR](https://www.cincinnatichildrens.org/research/divisions/c/cpir). 

## Table of contents:

1. [Setup](#setup)

2. [Usage](#Usage)

3. [Acknowledgments](#acknowledgements)


## Setup
### Run Standalone Application
To install XIPline, download the [`XIPlineInstaller.exe`](https://drive.google.com/drive/folders/1CIFWiEEHiJ0wWNwVQTAFq40OrirYuE7q?usp=drive_link) file. and fellow the normal instllation steps. This step doesn't require additional steps except downloading the pre-trained deep-learning models for auto-segmentation (see below). 

#### auto-segmentation: Downloading the .h5 models for machine learning
1 - Download the h5 models for auto segmentation from [here](https://drive.google.com/drive/folders/1CIFWiEEHiJ0wWNwVQTAFq40OrirYuE7q?usp=drive_link).

2 - Create this path `C:\XIPline\models` in your local `C:\` drive and add the downloaded .h5 files in the `models` folder. 

3 - Install [Python](https://www.python.org/downloads/) 3.10.2 and make sure to add the Python path to the environment variables.

Please refer to `User_Manual.pdf` document for more detailed steps. 

#### For GE Data: Waveform Files
1 - Create this path `C:\XIPline\GE_WFs` (in the same folder above) in your local `C:\` drive and add text files with the names (cal_wf_freq.txt, cal_wf_tr.txt, vent_wf.txt, diff_wf.txt, diss_wf_freq.txt) in the `GE_WFs` folder. 

2- Add the path to the waveform file in the first line of each file. Note: the waveform files should be .mat for vent_wf and diff_wf and .fdl for the rest. 

Please refer to `User_Manual.pdf` document for more detailed steps. 

### Run Developer Mode 

To customize the application, the user will need MATLAB version R2023b or newer. 

Clone XIPline to your local directory.

```
git clone https://github.com/aboodbdaiwi/XIPline.git
```

Start by adding the local folder in the MATLAB path and execute the `XIPline.mlapp` file. 

When implementing new features, or debugging, we recommend using the debugging MATLAB script `XIPline_Code_Script.m`. This script can run the entire application but without the graphical user interface (GUI). 

## Usage
These are a few analysis demos (please refer to the user's manual for step-by-step guide on performing all analysis)

1 - [Installation demo](https://www.youtube.com/watch?v=mWbWL6vIEUc&t=8s&ab_channel=AbdullahBdaiwi)

2 - [Calibration Analysis demo](https://www.youtube.com/watch?v=x1zQrBrFOZ8&ab_channel=AbdullahBdaiwi)

3 - [Ventilation Analysis demo](https://www.youtube.com/watch?v=qLTG6Hiz-q4&ab_channel=AbdullahBdaiwi)

4 - [Diffusion Analysis demo](https://www.youtube.com/watch?v=kItn_P4dDyw&ab_channel=AbdullahBdaiwi)

5 - [Gas Exchange Analysis demo](https://www.youtube.com/watch?v=_aerEFhWbm0&ab_channel=AbdullahBdaiwi)

6 - [Developer Mode demo](https://www.youtube.com/watch?v=fEjruhWYejA&t=494s&ab_channel=AbdullahBdaiwi)

## Acknowledgments:
Author: Abdullah S. Bdaiwi 

Work email: abdullah.bdaiwi@cchmc.org

Personal email: abdaiwi89@gmail.com

Website: [CPIR](https://www.cincinnatichildrens.org/research/divisions/c/cpir)

[please cite this paper](https://onlinelibrary.wiley.com/doi/full/10.1002/mrm.30347): Bdaiwi AS, Willmering MM, Plummer JW, et al. 129Xe Image Processing Pipeline: An open-source, graphical user interface application for the analysis of hyperpolarized 129 Xe MRI. Magn Reson Med.2024;1-18. doi: 10.1002/mrm.30347

