# MATLAB Development with Orchestra SDK

This document provides instructions for installing and using the MATLAB extension provided by the Orchestra SDK. The reconstruction content is provided in a single MATLAB MEX file, making installation as easy as adding the file to your MATLAB path. The MEX utility has an inventory of reconstruction content that can be called with string commands. Within the package, there are example m-file reconstructions to help you get started with the utility.


## Requirements

The MATLAB toolbox is supported on Windows, Mac, and Linux platforms. It has been built and tested with the following MATLAB versions:

- **Linux**: Built with MATLAB 2013a, tested with MATLAB 2022a.
- **Windows**: Built MATLAB 2018a, tested with MATLAB 2024a.
- **Mac**: Buit with MATLAB 2019b, tested with MATLAB 2024a.

MATLAB does not guarantee compatibility of MEX functions across versions, but they typically work with newer versions.

#### Recommended Versions

- **Linux**: MATLAB 2018a or later.
    - OpenSUSE/SLES 15.3 or later.
    - Ubuntu 20.04 or later.
- **MacOS**: MacOS 11+ with MATLAB 2019b or later.
- **Windows**: Windows 10/11 with MATLAB 2018a or later.

## Installation

The MATLAB extension is distributed as a ZIP file provided by GE HealthCare. Follow these steps to install:

1. Obtain the ZIP file: `orchestra-sdk_3.0.3_matlab.zip`.
2. Extract the ZIP file to a directory of your choice. This directory will be referred to as `MATLAB_SDK_ROOT`.
3. Add the `MATLAB_SDK_ROOT` directory (including subdirectories) to your MATLAB path. This allows you to access the example reconstructions and utility functions.

### Verifying Installation

1. Open MATLAB and type the following command:
    ```matlab
    >> GERecon
    ```
2. If you see a usage message with a list of available functions, the installation is successful.

### Troubleshooting on Windows

If you encounter an "Invalid MEX file" error, you may need to install the required Visual Studio Runtime libraries:

- For SDK versions 1.8-1 and later, download and install the **Microsoft Visual C++ Redistributable for Visual Studio 2017** (`vc_redist.x64.exe`).
- For older SDK versions, download and install the **Microsoft Visual C++ Redistributable for Visual Studio 2012**.

Refer to the following links for downloads:
- Visual Studio 2017: [Microsoft Support](https://support.microsoft.com/en-in/help/2977003/the-latest-supported-visual-cdownloads)
- Visual Studio 2012: [Microsoft Download Center](http://www.microsoft.com/en-pk/download/details.aspx?id=30679)

## Getting Started

The `MATLAB_SDK_ROOT` directory contains the following resources to help you get started:

- **MEX Utility**: The main utility (`GERecon.mex*`) is platform-specific:
    - `w64` for Windows
    - `maci64` for Macs with an Intel processor
    - `maca64` for Macs with an Apple processor
    - `a64` for Linux
- **Sub-directories**:
    - `Doc`: HTML-based documentation.
    - `Examples`: Example MATLAB scripts (`.m` files) that can be run directly.
        - The example data is available as a separate download: `orchestra-sdk_3.0.3_matlab_examples.zip`. After downloading, unzip the contents and place them in the `Examples/Data` directory.
    - `Scripts`: Additional MATLAB scripts for reconstruction, serving as a starting point for development.

Ensure the `GERecon` utility is on your MATLAB path to call it from the command line.
