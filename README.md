![XIPlineMainLogo_noBackground](https://github.com/aboodbdaiwi/HP129Xe_Analysis_App/assets/36932337/ee751c64-065b-4d8c-94b3-5edc89e14ac1)

# 129Xenon Image Processing Pipeline (XIPline)

XIPline application is developed at the [CPIR](https://www.cincinnatichildrens.org/research/divisions/c/cpir). 

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.16053654.svg)](https://doi.org/10.5281/zenodo.18124693)

## Table of contents:

1. [Running the Application](#Setup)

2. [Usage](#Usage)

3. [Acknowledgments](#acknowledgements)

3. [Updates](#Updates)


## Setup

### Standalone Mode

Due to ongoing improvements and updates, we encourage users to run the application in **developer mode**. If you prefer to use the application as a standalone executable, please contact the authors to obtain the `XIPlineInstaller.exe` installer.


### Developer Mode (Recommended)

To customize the application, the user will need MATLAB version R2023b or newer. 

Clone XIPline to your local directory.

```
git clone https://github.com/aboodbdaiwi/XIPline.git
```

Start by adding the local folder in the MATLAB path and execute the `XIPline.mlapp` file. 

When implementing new features, or debugging, we recommend using the debugging MATLAB script `XIPline_Code_Script.m`. This script can run the entire application but without the graphical user interface (GUI). 

##### Set up XIPline folder

1. Download the [`XIPline`](https://zenodo.org/records/16053654) folder.
   
2. Unzip the folder and place it directly in the `C:\` drive **without renaming or moving** it.
   
3. This directory includes all required components:
   - AI segmentation models
   - Offline reconstruction scripts

These components are essential for image reconstruction, registration, and segmentation.

##### Set up Python

1. Install [Python](https://www.python.org/downloads/) 3.10.2.

2. Add the Python path to the environment variables.
   
3. Navigate to `C:\XIPline\python_path.txt` and change the Python path to your local location.

Note: If you're having problems with the auto-segmentation, please take a look at [`issue #3`](https://github.com/aboodbdaiwi/XIPline/issues/3) for potential solutions. 

##### For GE Data: Waveform Files
 Add waveform files in the `C:\XIPline\GE\waveforms` folder. Please refer to `User_Manual.pdf` document for more detailed steps. 


## Usage
These are a few analysis demos (please refer to the user's manual for step-by-step guide on performing all analysis)

1 - [Installation demo (Please refer to the instructions above for the latest version)](https://www.youtube.com/watch?v=mWbWL6vIEUc&t=8s&ab_channel=AbdullahBdaiwi)

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

If you use this software, please cite:
> Bdaiwi AS, Willmering MM, Plummer JW, et al. *129Xe Image Processing Pipeline: An open-source, graphical user interface application for the analysis of hyperpolarized 129Xe MRI*. Magn Reson Med. 2024;1–18. https://doi.org/10.1002/mrm.30347
> 
 

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.16053654.svg)](https://doi.org/10.5281/zenodo.18124693)

## Updates:

| **Task**                                                                 | **Date**        |
|--------------------------------------------------------------------------|-----------------|
| Update Vent AI segmentation model                                        | Jan 01, 2026    |
| Update Diffusion Analysis (SNR, healthy ref, report)                     | Dec 11, 2025    |
| Update AI segmentation model for Diffusion Analysis                      | Oct 20, 2025    |
| Add new ventilation analysis report                                      | July 17, 2025   |
| Add offline reconstruction for 2D spiral Philips ventilation images      | July 17, 2025   |
| Add ANTs and manual registration for ventilation analysis                | July 17, 2025   |
| Remove Python requirement for AI segmentation models                     | July 17, 2025   |
| Update AI segmentation models                                            | May 20, 2025    |

## Analysis Report Examples:
![Vent_report](https://github.com/aboodbdaiwi/XIPline/blob/main/XIPline_resources/Vent_report_example.png?raw=true)

![Diff_report](https://github.com/aboodbdaiwi/XIPline/blob/main/XIPline_resources/Diff_report_example.png?raw=true)

![Gx_report](https://github.com/aboodbdaiwi/XIPline/blob/main/XIPline_resources/Gx_report_example.png?raw=true)

