# GE Recon Orchestra

- This repo has been initialized with the Windows Matlab version of the 3.0.3 GE Orchestra code from weconnect
    - Windows was chosen due to the majority of CPIR using Windows computers
    - Matlab was chosen since GE states the most up-to-date code and most features are in the c++ files, followed by Matlab, then Python
    - Future commits may add the others versions to this repo
 - Below is the Readme for the Orchestra 3.0.3 SDK
 - Make sure to add the GE_Recon_Orchestra and fidall folders are added to the path

---

# Orchestra Software Development Kit

The Orchestra SDK provides tools for GE HealthCare research customers to develop reconstruction software from MR system raw data, supporting C++, MATLAB, and Python development across Linux, Mac, and Windows platforms.

For detailed instructions, refer to the language-specific README files available for C++, MATLAB, and Python in their respective packages.

## Table of Contents
- [Feedback](#feedback)
- [Overview](#overview)
    - [C++ Development](#c-development)
    - [MATLAB Development](#matlab-development)
    - [Python Development](#python-development)
    - [Development Environment](#development-environment)
    - [Hardware Requirements](#hardware-requirements)
    - [Installation Requirements](#installation-requirements)
    - [Support](#support)
- [Version 3.0.3 Release Notes – May 10, 2025](#version-303-release-notes--may-10-2025)
    - [New Content](#new-content)
    - [Compatibility](#compatibility)
- [Frequently Asked Questions (FAQ)](#frequently-asked-questions-faq)
    - [Can I uninstall the SDK?](#can-i-uninstall-the-sdk)
    - [How do I include content from the SDK into my code?](#how-do-i-include-content-from-the-sdk-into-my-code)
    - [How do I generate private fields in the DICOM images generated from my reconstruction?](#how-do-i-generate-private-fields-in-the-dicom-images-generated-from-my-reconstruction)
    - [Registering Private Tags](#registering-private-tags)
- [Legal](#legal)
    - [GE HealthCare Proprietary Information](#ge-healthcare-proprietary-information)
    - [Trademark Legends](#trademark-legends)

## Feedback

We value your input! If you have suggestions, questions, or feedback about the Orchestra SDK, please don't hesitate to reach out. Your insights help us improve and better serve your needs. Join the discussion in the Orchestra subgroup of the GE HealthCare MR Collaboration Community on [WeConnect](https://weconnect.gehealthcare.com/).

## Overview

The Orchestra Software Development Kit (SDK) is a software package that is available to research-focused GE HealthCare customers. The tools in the package facilitate the development of reconstruction software from raw data captured by a variety of MR systems. The software allows for development of **C++**, **MATLAB**, and **Python** programs.

### C++ Development

For C++ development, the software comes in the form of C++ libraries and headers that can be used in the compilation of a custom reconstruction. Example source files are also provided. The C++ development environment is available for **Linux** and **Mac**.

### MATLAB Development

The MATLAB functionality comes in the form of a **MATLAB Extension (MEX)** file that allows for compiled C++ routines to be called like ordinary MATLAB functions. A set of **m-files** are provided as example reconstructions, along with sample data to run a variety of routines out of the box. The MEX utility is available in **Windows**, **Linux**, and **Mac**.

### Python Development

The Python functionality comes in the form of a `GERecon.so` file for **Linux** and **Mac** development or a `GERecon.pyd` file for **Windows** development that allows for compiled C++ routines to be called like ordinary Python functions. A set of **.py-files** are provided as example reconstructions, along with sample data to run a variety of routines out of the box. The Python utility is available in **Windows**, **Linux**, and **Mac**.

### Development Environment

Development of the software is performed offline on a **Linux**, **Mac**, or **Windows** workstation or virtual machine. The MR system may not be used for development or for building reconstruction programs.

### Hardware Requirements

Hardware required for development work may be purchased from any 64-bit **Windows**, **Mac**, or **Linux**-compatible personal computer (PC) distributor or a third-party vendor thereof. GE HealthCare does not sell and will not provide support for the **Windows**, **Mac**, or **Linux** operating systems (OS), PCs, or software. Customers may purchase a service contract from the OS provider or the PC distributor for support and maintenance of hardware and software obtained from these manufacturers.

### Installation Requirements

It is expected that the individual installing the software has the necessary background in setting up and maintaining PC equipment and should have some understanding of application software, including **Linux**, **Mac OS X**, or **Unix-based systems**. The installer must have the necessary super-user privileges and must be authorized to perform software installations on the machines that they are responsible for maintaining.

### Support

Technical questions on the internals of the software provided in this package may be directed to the GEHC SDK support team through membership in the Orchestra subgroup of the GE HealthCare MR Collaboration Community on [WeConnect](https://weconnect.gehealthcare.com/).

## Version 3.0.3 Release Notes – May 10, 2025

### New Content

#### SDK Package Restructuring

The Orchestra SDK has been restructured into distinct packages, each tailored for a specific operating system and language pair. Each package contains:

- The language and OS-specific SDK.
- This README, a CHANGELOG, and a LICENSE file.
- A language-specific README with more detailed install and usage instructions.
- Language-specific examples and starter scripts.

Example data have been placed in their own language-specific downloads to keep the main packages smaller.

#### Full Native Support for Apple Silicon

The Orchestra SDK now includes native ARM64 support across C++, MATLAB, and Python. Both C++ and Python SDKs are provided as universal builds, ensuring compatibility across architectures. The MATLAB SDK, however, is available in separate x86_64 Mac and ARM64 Mac versions due to MEX file constraints.

#### Portability to Ubuntu Environments

Compatibility issues with Ubuntu releases have been resolved. The Orchestra SDK now supports all languages on Ubuntu. Compatibility has been tested with Ubuntu 20.04 (Focal), 22.04 (Jammy), and 24.04 (Noble).

#### Python Module
You can now install the Python version of the Orchestra SDK using `pip`:

```bash
pip install /path/to/orchestra-sdk_3.0.3_python/GERecon
```

This command integrates the Orchestra SDK into your Python environment, automatically resolving and installing all necessary dependencies.

If you prefer, you can still circumvent `pip` and use the Python library directly by adding `/path/to/orchestra-sdk_3.0.3_python/GERecon/GERecon/GERecon.so` to your `$PYTHONPATH`.

#### Other Features & Enhancements

- Added support for Propeller archives and control packets in Python.
- Added raw data output for raw control packets in MATLAB and Python.

### Compatibility

#### C++

- Support for reading raw data from systems with header revision 30.100
- Recommended Versions:
    - **Linux**: CMake 3.21 or later, glibc 2.31 or later.
    - **MacOS**: CMake 3.21 or later, Xcode 

#### MATLAB

- Recommended Versions:
    - **Linux**: MATLAB 2018a or later.
        - OpenSUSE/SLES 15.3 or later.
        - Ubuntu 20.04 or later.
    - **MacOS**: MacOS 11+ with MATLAB 2019b or later.
    - **Windows**: Windows 10/11 with MATLAB 2018a or later.

#### Python

- Python 3.6 and 3.7 are **deprecated** and no longer supported.
- Recommended Versions:
    - **Linux**: Python 3.8 or later with OS:
        - OpenSUSE/SLES 15.3 or later.
        - Ubuntu 20.04 or later.
    - **MacOS**: MacOS 11+ with Python 3.10 or later.
    - **Windows**: Windows 10/11 with Python 3.10. 

## Frequently Asked Questions (FAQ)

### Can I uninstall the SDK?

Yes, uninstalling the SDK is straightforward:

- **C++**: Simply remove the installation directory.
- **MATLAB**: Delete the `GERecon.mex*` file.
- **Python**: Run `pip uninstall GERecon`, or manually delete the `GERecon.so` or `GERecon.pyd` file if `pip` was not used for installation.

### How do I include content from the SDK into my code?

When modifying or creating new C++ header or source files, you can include a header file by specifying its relative path starting from `$SDKTOP/recon`. For example, to include the file located at `$SDKTOP/recon/Orchestra/Cartesian2D/Homodyne.h`, use the following syntax:

```cpp
#include <Orchestra/Cartesian2D/Homodyne.h>
```

The library containing the implementation should already be included in the build. If not, you can modify the `CMakeLists.txt` file to explicitly include the library. This can also be useful for including header file paths or libraries that are not part of the SDK.

### How do I generate private fields in the DICOM images generated from my reconstruction?

It is often desirable to add custom or user-defined information to a DICOM image. This allows any arbitrary information to be stored with the image and retrieved for later use. In the [DICOM Standard](https://www.dicomstandard.org/), these are known as Private Tags.

The concept is fairly straightforward, and for the most part, Private Tags behave like the standard ones. However, a few rules must be followed to make sure any generated DICOMs stay compliant and do not interfere with any GE HealthCare or external devices.

The first requirement is that the group number (gggg) in the standard notation of a group/element tag (gggg,eeee) **MUST** be odd. From [Part 5, Section 7.8.1](https://dicom.nema.org/medical/dicom/current/output/html/part05.html#sect_7.8) of the DICOM Standard, the rest of the requirements state:

> ------
> 1. Private Creator Data Elements numbered (gggg,0010-00FF) (gggg is odd) shall be used to reserve a block of Elements with Group Number gggg for use by an individual implementer. The implementer shall insert an identification code in the first unused (unassigned) Element in this series to reserve a block of Private Elements. The VR of the private identification code shall be LO (Long String) and the VM shall be equal to 1.
> 
> 2. Private Creator Data Element (gggg,0010), is a Type 1 Data Element that identifies the implementer reserving element (gggg,1000-10FF), Private Creator Data Element (gggg,0011) identifies the implementer reserving elements (gggg,1100-11FF), and so on, until Private Creator Data Element (gggg,00FF) identifies the implementer reserving elements (gggg,FF00-FFFF).
> 
> 3. Encoders of Private Data Elements shall be able to dynamically assign private data to any available (unreserved) block(s) within the Private group, and specify this assignment through the blocks corresponding Private Creator Data Element(s). Decoders of Private Data shall be able to accept reserved blocks with a given Private Creator identification code at any position within the Private group specified by the blocks corresponding Private Creator Data Element.
>
>     - The versions of this standard prior to V3.0 described shadow groups. These were groups with a group number one greater than the standard groups. Elimination of conflicts in Private Data Element Tags have made this distinction obsolete and this terminology has been retired.
>
>     - The versions of this standard prior to V3.0 specified private group element numbers (gggg,10FF-7FFF) reserved for manufacturers and private group element numbers (gggg, 8100-FFFF) reserved for users. Elimination of conflicts in Private Data Element Tags has made this distinction obsolete and this specification has been retired.
> 
>     - The requirements of this section do not allow any use of elements in the ranges (gggg,0001-000F) and (gggg,0100-0FFF) where gggg is odd.
> 
> 4. Elements with Tags (0001,xxxx), (0003,xxxx), (0005,xxxx), and (0007,xxxx) shall not be used.
>
> ------

So what does this mean? It comes down to registering your private tags (with a Private Creator ID) in an arbitrary odd group. The group (gggg) needs to be an odd number greater than 0008. The Private Creator is encoded as such:

```
(gggg,00bb)
```

Where `bb` is an open private block in the DICOM object and has a value in the range of `10-FF`. This private block is where conflicts between vendors are dealt with. The actual private elements must then be assigned in this block using the following convention:

```
(gggg,bbxx)
```

But how do you tell if a block is open? All GE HealthCare MR DICOMs are created using the block number `10`. Therefore, any odd group number greater than or equal to `0009` and any block number greater than or equal to `11` can be used to define private tags with the toolkit. If using DICOM tags from another vendor, use a different block number from your own.

Normally in DICOM, block numbers can be changed when rewriting an image, but in this toolkit, GE HealthCare always uses block `10`. It is suggested that you use your institution name or abbreviation in the Private Identification code (LO) string to ensure uniqueness. Do not use strings that could interfere with GE HealthCare’s, such as `GEMS`, `GEHC`, etc.

## Legal

### GE HealthCare Proprietary Information

Copyright © 2025 GE HealthCare. All rights reserved.

The Orchestra SDK software, documentation, and all copies of it are the proprietary property of GE HealthCare, and title to it remains with GE HealthCare. This software and documentation may be used only for the purpose of development of Magnetic Resonance research experiments. This software and documentation may not be used for commercial purposes. The contents of this package contain proprietary, confidential information and trade secrets of GE HealthCare. No use may be made of this package or any software developed with its use, except for non-commercial purposes at the institution to which GE HealthCare ships the package.

Without limiting the rights under copyright, no part of this document may be reproduced, stored in or introduced into a retrieval system, or transmitted in any form or by any means (electronic, mechanical, photocopying, recording, or otherwise), or for any purpose, without the express written permission of GE HealthCare.

GE is a trademark of General Electric Company used under trademark license.

### Trademark Legends

Linux® is the registered trademark of Linus Torvalds in the United States and other countries.

MATLAB is the registered trademark of The MathWorks, Inc. in the United States and other countries.

Windows is a registered trademark of Microsoft Corporation in the United States and other countries.

Mac is a trademark of Apple Inc., registered in the U.S. and other countries.

CMake is a trademark of Kitware Inc, registered in the U.S. and other countries.

All other trademarks and copyrights referred to are the property of their respective owners.
