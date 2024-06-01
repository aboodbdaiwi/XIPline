![XIPlineMainLogo_noBackground](https://github.com/aboodbdaiwi/HP129Xe_Analysis_App/assets/36932337/ee751c64-065b-4d8c-94b3-5edc89e14ac1)

# 129Xenon Image Processing Pipeline (XIPline)

XIPline application is developed at the [CPIR](https://www.cincinnatichildrens.org/research/divisions/c/cpir). 

## Table of contents:

1. [Setup](#setup)

2. [Usage](#Usage)

3. [Acknowledgments](#acknowledgements)


## Setup
### Run Standalone Application
To run the standalone application, download the `XIPline.exe` file. This step doesn't require additional steps except downloading the pre-trained deep-learning models for auto-segmentation (see below) and the MATLAB runtime. 
Download the MATLAB runtime compatible with the application's MATLAB version (R2023b) from [here](https://www.mathworks.com/products/compiler/matlab-runtime.html)

You can also clone the application to your local directory for future releases.
```
git clone https://github.com/aboodbdaiwi/HP129Xe_Analysis_App.git
```

#### auto-segmentation: Downloading the h5 models for machine learning
1 - Download the h5 models for auto segmentation from [here](https://drive.google.com/drive/folders/1gcwT14_6Tl_2zkLZ_MHsm-pAYHXWtVOA?usp=sharing).

2 - Create this path `C:\XIPline\models` in your local `C:\` drive and add the downloaded h5 files in the `models` folder. 

3 - Install [Python](https://www.python.org/downloads/) 3.10 and install these packages: 
```
pip install numpy
pip install keras
pip install tensorflow
pip install nibabel
pip install scipy
```
Please refer to `User_Manual.pdf` document for more detailed steps. 

### Run Developer Mode 
To customize the application, the user will need MATLAB version R2023b or newer. Start by adding the local folder in the MATLAB path and execute the `XIPline.mlapp` file. 
When implementing new features, or debugging, we recommend using the debugging MATLAB script `XIPline.m`. This script can run the entire application but without the graphical user interface (GUI). 

## Usage
These are a few analysis demos (please refer to the user's manual for step-by-step quid on performing all analysis)
1 - Calibration Analysis demo:
video
2 - Ventilation Analysis demo:
video
3 - Diffusion Analysis demo:
video
4 - Gas Exchange Analysis demo:
video

## Acknowledgments:
Author: Abdullah S. Bdaiwi

Work email: abdullah.bdaiwi@cchmc.org

Personal email: abdaiwi89@gmail.com

Website: [CPIR](https://www.cincinnatichildrens.org/research/divisions/c/cpir)

please cite this paper: TBA

