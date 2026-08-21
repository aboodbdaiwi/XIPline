# Release Notes

## 3.0.3 – May 10, 2025

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

--------------------------------------------------------------------------------

## 2.1-1

### New Content

- Support for reading raw data from systems with header revision 30.000
- Support for reading PRICE raw data in Python and MATLAB SDK
- Native ARM64 build for Mac Python SDK
- Non-sequential data access to Scan Archives (acquired in 3D) in Python SDK using the `Cartesian3DAcquiredData`.

#### Accessing Non-Sequential Data from Cartesian3D ScanArchive in Python

1. The function calls which return the k-space information for a volume and scan from an acquired data object:
    - `volumeData = acquiredData.GetVolume(acquisition, echo)`
    - `scanData = acquiredData.GetScan()`

    **Example:**
    To read data for acquisition count 1 and echo count 2, invoke:

    ```python
    acquiredData = Cartesian3DAcquiredData(archive)
    volume1 = acquiredData.GetVolume(0, 1)  # 1st acquisition, 2nd echo
    ```

--------------------------------------------------------------------------------

## 2.0-1

### New Content

- Support for reading raw data from systems with header revision 28.003
- Non-sequential data access to Scan Archives in Python SDK using the `Cartesian2DAcquiredData`.
- Non-sequential data access to Scan Archives in C++ SDK (Rehearsals)
- Compiler upgrade to gcc6 to enable C++ 11 features for Linux C++ SDK. (All the libraries compiled for Linux C++ SDK were compiled on gcc6)

#### Compiler Upgrade for Linux C++ SDK

All the details of the compiler upgrade have been discussed in Chapter 2 and Chapter 3.

#### Accessing Non-Sequential Data from Cartesian2D ScanArchive in Python

1. The function calls which return the k-space information for a slice, baseline, channel, and scan from an acquired data object:
    - `channelData = acquiredData.GetChannel(acquisition, slice, echo, channel)`
    - `sliceData = acquiredData.GetSlice(acquisition, slice, echo)`
    - `baselineData = acquiredData.GetBaseline(acquisition, slice, echo)`
    - `scanData = acquiredData.GetScan()`

    **Example:**
    To read data for acquisition count 1, slice count 3, echo count 1, and channel count 4, invoke:

    ```python
    acquiredData = Cartesian2DAcquiredData(archive)
    channel4 = acquiredData.GetChannel(0, 2, 0, 3)  # 1st acquisition, 3rd slice, 1st echo, 4th channel
    ```

#### Accessing Non-Sequential Data from Cartesian2D ScanArchive in C++

1. The function calls which return the k-space information for a slice, baseline, channel, and scan from an acquired data object:
    - `MDArray::ComplexFloatCube sliceData; acquiredData.GetSlice(sliceData, acquisition, echo, slice);`
    - `MDArray::ComplexFloatCube baselineData; acquiredData.GetBaseline(baseline, acquisition, echo, slice);`
    - `MDArray::ComplexFloatMatrix channelData; acquiredData.GetChannel(channelData, acquisition, echo, slice);`
    - `MDArray::Array<std::complex<float>, 6> scan; acquiredData.GetScan(scan);`

    **Example:**
    To read data for acquisition count 1, slice count 3, echo count 1, and channel count 4, invoke:

    ```cpp
    const Cartesian2DAcquiredData acquiredData(scanArchive);
    MDArray::ComplexFloatMatrix channelData;
    acquiredData.GetChannel(channelData, 0, 2, 0, 3);  // 1st acquisition, 3rd slice, 1st echo, 4th channel
    ```

--------------------------------------------------------------------------------

## 1.10-1

### New Content

- Support for reading raw data from systems with header revision 28.002
- Python support for reading slice-information object from Pfile and ScanArchive.
- From 1.10-1 and later releases, the Windows MEX utility will be supported only on Matlab release 2019 and above.

#### Reading Slice Information Object in Python

1. The function call which returns the slice-information object from an input Pfile or ScanArchive object:
    - `slice_info = pfile.Info(slice)`
    - `slice_info = pfile.Info(slice, pass)`
    - `slice_info = archive.Info(slice)`
    - `slice_info = archive.Info(slice, pass)`

    where `slice_info` is a Python dictionary.

#### Examples

To read and print the slice-information for slice count 1 and phase count 1, invoke:

```python
pfile = Pfile(filename)
info = pfile.Info(0, 0)
print(info)
```

Output:

```python
{'freqShift': 0.0, 'slthickScale': 1.0, 'sliceShift': 0.0, 'freqScale': 1.0, 'phaseScale': 1.0, 'sliceNumber': 0, 'phaseShift': 0.0}
```

--------------------------------------------------------------------------------

## 1.9-1

### New Content

- Python support for reading raw data headers (ScanArchive & Pfile)
- Python PURE support
- Mac OS X & Windows support for Python development

### Bug Fix

- EPI Phase Correction bug fix in Python and MATLAB for body data collected on Signa Premier

**Note:** No new enhancements for the Linux & Mac C++ development environments are added in the 1.9-1 release. Users are requested to download the following for the C++ development:
- `orchestra-sdk-1.8-1.x86_64.tgz` (for Linux)
- `orchestra-sdk-1.8-1.mac64.tgz` (for Mac)

--------------------------------------------------------------------------------

## 1.8-1

### New Content

- Support for reading raw data from systems with header revision 28.000
- From 1.8-1 and later releases, the Windows MEX utility will be supported only on Visual Studio 2017 Runtime libraries. All older versions before 1.8-1 of the Windows MEX utility will be supported only on Visual Studio 2012 Runtime libraries.

Refer to Chapter 6 (MATLAB Development -> Installation) for instructions to download and install the Visual Studio 2012 and 2017 Runtime libraries.

#### Python Development in Linux

Python access has been added for a number of reconstruction pipelines in this release. The ability to read raw data from either a Pfile or ScanArchive is available for all reconstructions listed below, with the exception of EPI relying solely on ScanArchive input. DICOM read and write functionality is also made available in Python for this release.

- 2D Cartesian Reconstruction
- Calibration Reconstruction
- ARC Reconstruction
- ASSET Reconstruction
- EPI Reconstruction (including EPI Diffusion)

--------------------------------------------------------------------------------

## 1.7-1

### New Content

- Support for reading raw data from systems with header revision 27.001
- Support for reading and processing ScanArchive files is incorporated into the Diffusion C++ Rehearsal
- Multiband reconstruction logic is incorporated into MATLAB and C++ examples for fMRI.
- MATLAB & C++ support for DW-EPI distortion correction, motion correction, and Eddy current correction, which are based on Reversed Phase Gradient (RPG) acquisition, rigid body registration, and polynomial registration, respectively.
- MATLAB & C++ support for High Order Phase Correction for EPI.

#### Distortion Correction Usage in MATLAB

1. The function call which computes the B0 displacement map for the given forward and reverse polarity image:

    ```matlab
    displacementMap = GERecon('DistortionCorrection.ComputeMap', forwardImage, reverseImage);
    ```

2. The function which applies the displacement map and/or does motion correction on the distorted image to give corrected images:

    ```matlab
    correctedImage = GERecon('DistortionCorection.Apply', uncorrectedImage, displacementMap, referenceImage);
    ```

3. The function to set the parameters required for the RPG algorithm:

    ```matlab
    GERecon('DistortionCorrection.SetParameters', kernelWidth, zBlurringFactor, kernelStep, regularizationMultiplier, residualCompareFactor, maxIter, cubicInterpolation);
    ```

    **Note:** This function is optional and can be used to override default algorithm parameters.

#### HOPC Usage in MATLAB

1. The function call which computes phase correction coefficients and returns the computed coefficients:

    ```matlab
    coefficients = GERecon('Epi.ComputeCoefficients', kSpace2D, doHighOrderPhaseCorrection);
    ```

    - `doHighOrderPhaseCorrection`: 0 - linear phase correction, 1 - high order phase correction

2. The function which will use the above computed coefficients to apply on a k-space data:

    ```matlab
    phaseCorrected2DKSpace = GERecon('Epi.ApplyPhaseCorrection', coefficients, kSpace2D);
    ```

    **Note:** K-space data is a 2D data per channel per slice.

--------------------------------------------------------------------------------

## 1.6-1

### New Content

- **CMake C++ Build Configuration**
- **C++ libraries for reading ScanArchive files**  
    (See 1.5-2 release notes for a description of the ScanArchive file format.)
- **Embedded file extraction from ScanArchive files**  
    - In MATLAB, embedded files are extracted to the directory that the `GERecon` function was first invoked from.  
    - In the C++ examples, embedded files are extracted to the ScanArchive directory.  
    - Configurability for this extract location may be added in a future release.
- **Multiband reconstruction logic**  
    Incorporated into the C++ MultiPhaseEpi Example.
- **MATLAB & C++ Pfile examples**  
    For 2D spiral reconstruction.
- **MATLAB & C++ ScanArchive examples**  
    For Spiral2D, 3DASL, & Epi reconstruction.

### New CMake C++ Build Configuration

The C++ build has been moved to a CMake-based build system. This allows the C++ build to be portable across platforms and IDEs. The `reconmake` command has been retired. To configure a C++ build using CMake, see Chapter 3 (Linux) or Chapter 4 (Mac) below.

### Modified ScanArchive API

This release also contains a modification to the structure returned by the `Archive.Next` command in MATLAB. The raw data associated with a control packet is moved up one level in this structure. The raw data is no longer in the `Frame` substructure; it is now accessible from the `Data` field in the top level of the structure. Furthermore, the `Frame` substructure is renamed to `FrameInfo`.

The raw data in the `Data` field has dimensions of `[xRes x numChannels x numFrames]`. Many control packets have a single frame (readout) of data associated with them. For these cases, the third dimension of the raw data matrix will have a size of 1. For cases in which multiple frames are associated with a single control packet, the third dimension of the raw data matrix will be equal to the number of frames associated with the control packet.

#### Example Output from the `Archive.Next` Command

```matlab
>> control = GERecon('Archive.Next', archive)

control =

    Type: 'HyperFrameControlPacket'
    opcode: 6 
    echoNum: 0
    operation: 0
    cardRepeat: 1
    kRead: 1 
    sliceNum: 0
    viewNum: 40
    viewSkip: -1
    numViews: 40
    FrameCount: 40
    FrameInfo: [1x40 struct]
    Data: [172x18x40 single] 
```

--------------------------------------------------------------------------------

## 1.5-2

This release contains an update to the MATLAB package to add functionality to read a new type of raw file generated by systems with DV26 software (and beyond). These raw files are known as ScanArchives. They do not replace Pfiles, but rather provide an alternative and unsorted form of acquisition data which may better serve some use cases.

Below is a quick introduction to ScanArchive: where to find them, what they contain, and how to read them. As with any new functionality, there is still a lot to enhance and tools to develop. Please give them a try and let us know what you think. User feedback will play a big part in the evolution of ScanArchives.

For a short ScanArchive reconstruction, see the example provided in the package under `<TOP>/Examples/RunCartesianScanArchive.m`.

### Where to Find ScanArchives

Almost every scan that is run will generate a unique ScanArchive. The system does NOT need to be in research mode nor is a CV setting required to save the raw data. Unlike Pfiles, ScanArchives are saved directly on the VRE to a partition that is allocated only for raw data. This disk space is managed automatically, assuring the latest scans are always available on disk. The generation of hardware the scanner has will determine how much and for how long the data will remain on the system. In general, this can range from a few hours to a few days’ worth of scans.

The archived data files will be saved on the VRE under the directory: `/data/arc/`. To log in to the VRE from the host, simply type: `rsh vre`. Data can be copied to a mounted drive on the host and exported in the same mechanism used for Pfiles.

Once a series is finished, the ScanArchive file will be sorted into sub-directories based on the corresponding Exam and Series numbers. This should make locating files much simpler. ScanArchives that have completed or been closed will be located under:

```
/data/arc/Closed/Exam##/Series##/ScanArchive_<ID>_<TIMESTAMP>.h5
```

The ID is usually the system name, and TIMESTAMP is typically shown in `YYYYMMDD_HHMMSS###`. For example, a full archive may look like the following:

```
/data/arc/Closed/Exam10/Series1/ScanArchive_MR450_20170112_115629319.h5
```

When sorting through the directories for the ScanArchives, it’s likely that there will be multiple files in the same `Series#` directory. This would be the case if Prescan is run or if a series was duplicated on the UI. For a quick check of the ScanArchive of interest, you can parse the header using a tool on the scanner. Example:

```bash
# /usr/g/bin/ArchiveTool --input-file ScanArchive.h5 --list-info
```

Example output from an Autoshim and 3-Plane Localizer may look like the following:

```
Exam: 455, Series: 1, Seq: 1, Chans: 1, Coil: C-GE_HNS Head, File Type: ScanArchive, Scan Type: Autoshim, Rev: 0, Desc: Autoshim
Exam: 455, Series: 1, Seq: 2, Chans: 12, Coil: C-GE_HNS Head, File Type: ScanArchive, Scan Type: Scan, Rev: 27.000, Desc: 3-Plane Localizer
```

The same `ArchiveTool` can be used to anonymize the data before pulling it off the scanner, as well as some other handy utilities for inspecting and manipulating the files. For a full list of commands supported by the tool, type:

```bash
# /usr/g/bin/ArchiveTool --help
```

### ScanArchive Format

At their core, ScanArchives are HDF5 files. The open-source hierarchical data format provides a lot of functionality and tools that make creating, compressing, modifying, and configuring the files very easy. The files can be inspected with any HDF5 viewer, but the raw data is serialized in binary format and most of the headers are serialized as XML to better support versioning and portability. Therefore, ScanArchives will be most useful with the tools and APIs provided by the Orchestra SDK.

The data stored in the file is the raw and unsorted stream of acquisition data received by the hardware. The data is in time order as it was acquired, instead of traditional sorted K-Space order like most Pfiles. For those familiar with EPIC sequence development, the SSP/DAB packets filled and sent by the PSD are the controls that show up in the ScanArchive. The corresponding frames of data (where applicable) are matched to these PSD controls so the user has access to the raw frame and its header for everything played out in a scan.

The ScanArchive also contains all secondary inputs needed for retrospective reconstruction or processing of the raw data. Traditional secondary inputs include the gradient non-linearity coefficients file (`gw_coils.dat`), calibration data for ASSET/PURE, EPI parameter files, and so on. This, combined with all acquisition and timing information (including real-time updates), makes a complete retrospective reconstruction possible from a single ScanArchive file.

To further understand the format of the data and how controls and frames of data are stored, see the next section on reading ScanArchives.

### Reading ScanArchives

Currently, the method for reading ScanArchives is very basic. There is one call to load the file from disk and populate the header information and a second call to iterate over each control/frame of data (in acquired time order). The current iterator is only forward-looking but could be enhanced to allow for bi-directional or random access iterating.

#### Loading

The call to load a ScanArchive in MATLAB is:

```matlab
archiveHandle = GERecon('Archive.Load', <FullPathToArchive>)
```

For example, loading a simple Lx scan leads to the following:

```matlab
>> archive = GERecon('Archive.Load', 'ScanArchive.h5')

archive =

    struct with fields:

        HandleType: 'ScanArchive'
        ID: 1 
        ArchiveType: 'ScanArchive'
        ControlCount: 1569
        FrameCount: 1568
        UpdateCount: 0
        StudyUID: ' '
        ScanType: 'Scan'
        ReconName: ' '
        ExamNumber: 17
        SeriesNumber: 5 
        SequenceNumber: 9 
        DateTime: '2014-Sep-03 12:41:44.033231' 
        DateTimeLocalOffset: '05:30:00'
        ReconPath: ' '
        SystemID: [1x1 struct]
        EntryPoints: {'scan'}
        DownloadData: [1x1 struct]
        IsLxScan: 1
        Slices: 4
        Passes: 1
        SlicesPerPass: 4
```

All header information is provided directly on load. The `DownloadData` field in the structure contains the detailed “raw header” traditionally associated with an Lx scan (sometimes referred to as the Pool Header). The structure is a direct mapping of the C++ structure and its fields found in `<TOP>/lx/include/rdbm.h` and `imagedb.h`.

```matlab
>> archive.DownloadData

ans = 

    struct with fields:

        rdb_hdr_rec: [1x1 struct]
        rdb_hdr_per_pass: [1x2048 struct]
        rdb_hdr_unlock_raw: [1x2048 struct]
        rdb_hdr_data_acq_tab: [1x2048 struct]
        rdb_hdr_nex_tab: [1x1 struct]
        rdb_hdr_nex_abort_tab: [1x1 struct]
        rdb_hdr_tool: [1x1 struct]
        rdb_hdr_ps: [1x1 struct]
        rdb_hdr_exam: [1x1 struct]
        rdb_hdr_series: [1x1 struct]
        rdb_hdr_image: [1x1 struct]
        rdb_hdr_grad_data: [1x1 struct]
        rdb_hdr_cttentry: [1x1 struct]

>> archive.DownloadData.rdb_hdr_rec

ans =

    struct with fields:

        rdb_hdr_rdbm_rev: 27
        rdb_hdr_off_data: 219828 
        rdb_hdr_off_per_pass: 4096
        rdb_hdr_off_unlock_raw: 20480 
        rdb_hdr_off_data_acq_tab: 36864
        rdb_hdr_off_nex_tab: 184320 
        rdb_hdr_off_nex_abort_tab: 186372
        rdb_hdr_off_tool: 188424
        rdb_hdr_off_exam: 193660 
        rdb_hdr_off_series: 195620
        rdb_hdr_off_image: 198180
        rdb_hdr_off_ps: 190472 
        rdb_hdr_off_grad_data: 200628
        rdb_hdr_off_CTT_data: 214260
```

It should be noted that different scan types (Prescan, Autoshim, Spectro Prescan) will contain a `DownloadData` structure that matches the C++ structure used for its scan, not the Lx Pool Header. However, most work will be done with traditional Lx scans and match the format above.

#### Iterating

Because the iterator only supports the forward direction, there is one call to receive the “next” control/frame in the ScanArchive.

```matlab
control = GERecon('Archive.Next', archiveHandle)
```

For example, iterating the same Lx ScanArchive produces the following:

```matlab
>> control = GERecon( 'Archive.Next', archive)

control =

    struct with fields:

        Type: 'ProgrammableControlPacket' 
        opcode: 1
        echoTrainIndex: 0
        unused: [0 0 0 l]
        sliceNum: 0
        echoNum: 0
        operation: 0
        viewNum: 0
        FrameCount: 1
        Frame: [1x1 struct]

>> control.Frame

ans =

    struct with fields:

        SequenceNumber: 0
        SampleType: 'Extended'
        SampleSize: 8
        NumberSamples: 256
        ChannelCount: 8
        Header: [1x1 struct]
        Footer: [18 struct]
        Data: [256x8 int32]
```

By default, the control type is filled and tagged by the standard definitions used by PSD/Recon. In this case, the standard SSP programmable control packet is used to describe the first frame of data received in the scan. The Frame structure itself provides a full description of the acquisition data as well as the raw buffer in the `Data` field.

For advanced use cases or when the standard control definitions are not accurate (i.e., the data is not slice/echo/view), the control packets can be viewed in a raw mode. This allows the 31-bytes (plus 1 for opcode) to be viewed and parsed in any user-specified fashion. Simply add `raw` to the call above to see the data in an alternative format.

```matlab
>> control = GERecon( 'Archive.Next', archive, 'raw')

control =

    struct with fields:
        Type: 'RawPacket' 
        opcode: 1
        data: [0 0 0 0 0 1 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0] 
        sequenceNumber: 0
        FrameCount: 1
        Frame: [1x1 struct]
```

While iterating over the controls, it’s important to note that not every control will have frames of data associated with it. Traditional examples may include Plane 3D controls, Pass Done/Scan Done controls, and many others. Every Lx scan should end with a Scan Done control indicating the end of the scan. This can be used to end any program loops. If the iterator exceeds the number of controls, a `null` string is returned. An example of the typical last control packet would look like the following:

```matlab
>> control

control = 

    struct with fields:

        Type: 'ScanControlPacket' 
        opcode: 0
        scanControl: 12
        pathControl: 0
        opcode2: 0
        isAcqDone: 1
        isScanDone: 1
        FrameCount: 0
```

#### Closing

It is recommended to explicitly close the ScanArchive using `GERecon` when it is no longer being used. This will release it from memory and properly close the file. Other methods to close the file (e.g., `clear mex`) may result in undesired behavior. To close an archive, use:

```matlab
GERecon('Archive.Close', archiveHandle);
```

--------------------------------------------------------------------------------

## 1.5-1

### Changes to Pfile Content

Previously, only Pfiles with ARC enabled were Z-encoded. Starting with DV26 scanners, all 3D datasets are Z-encoded and must be Z-transformed and 3D scaled. For EPI scans with reference views, the reference views are now flipped and will have to be “unflipped.” There is an API for unflipping the reference views.

To account for the changes to DV26 Pfiles, new APIs are provided in both C++ and MATLAB. This will help describe the data in both the new and older formats so code can be written without knowledge of the software revision from which a given Pfile was obtained.

### New Convenience APIs

The Pfile object has some new APIs available for both C++ and MATLAB:

- **Convenience APIs to check if ARC is enabled**  
    - C++ API: `bool IsArc() const;`  
    - MATLAB API: the field in the MATLAB Pfile struct `isArc`

- **Convenience API to check if the Pfile is Z-Encoded**  
    - C++ API: `bool IsZEncoded() const;`  
    - MATLAB API: the field in the MATLAB Pfile struct `isZEncoded`

- **Convenience API to get the appropriate 3D Scaling Factor to apply after Z-transform**  
    - C++ API: `float ScaleFactor3D();`  
    - MATLAB API: the field in the MATLAB Pfile struct `scaleFactor3d`

- **Convenience API to check if the EPI Reference views are flipped**  
    - C++ API: `bool EpiReferenceViewsFlipped() const;`  
    - MATLAB API: the field in the MATLAB Pfile struct `areReferenceViewsFlipped`

- **Convenience API to return the number of slices actually reconstructed**  
    - C++ API: `int ReconstructedSlicesPerAcq() const;`  
    - MATLAB API: the field in the MATLAB Pfile struct `reconstructedSlicesPerPass`

### Modified APIs

The Pfile object modified/renamed an existing API to differentiate between acquired slices and reconstructed slices for both C++ and MATLAB:

- **Convenience API to return the number of slices acquired in K-space**  
    - C++ API: `int AcquiredSlicesPerAcq() const;`  
        _Note – this previously was `int SlicesPerAcq() const;`_  
    - MATLAB API: the field in the MATLAB Pfile struct `slicesPerPass`  
        _Note – this field has the same name in MATLAB as previous revisions, but might return a different value now._

### Rehearsals Updated (C++ and MATLAB)

The following rehearsals have been updated to account for the changes in Pfile content:

- **Cartesian**  
    - For C++, added a new Cartesian3D Rehearsal Application. The new example can be built and copied from:  
        `<SDKTOP>/recon/Orchestra/Cartesian3D/Rehearsal`  
        _Note that the former Cartesian 2D Rehearsal has been removed, since the new example handles both 3D and 2D Pfiles._  
        Removed: `<SDKTOP>/recon/Orchestra/Cartesian2D/Rehearsal`  
    - For MATLAB, split existing `RunCartesian` into three separate files: `RunCartesian2D`, `RunCartesian3D`, and `RunCartesian`. The latter determines which of the two former scripts to run and runs them.

- **Flex**  
    - Added the Z-transform step if Pfile is Z-encoded.

- **Calibration**  
    - For C++, added necessary Z-transforms for 3D calibration. No changes to MATLAB.

- **EPI**  
    - For C++ and MATLAB, added necessary checks for `ReferenceViewsFlipped`, and implemented the `RowFlip` plugin to flip accordingly.

### Linux Package Deployed as Tar File

The Linux C++ development environment is no longer deployed as an RPM file, but rather a simple Tar file. This allows any user to simply unzip the file to a desired location, update environment settings, and compile against Orchestra. This will also enable having multiple versions installed next to each other. The installation instructions below describe the new changes.

--------------------------------------------------------------------------------

## 1.4-722

Matlab only release to support Matlab 2015 and newer versions on the Mac platform.

--------------------------------------------------------------------------------

## 1.3-24

### Fat/Water (Flex) Reconstruction Functionality

APIs for doing Fat/Water processing were added in both MATLAB and C++. Complete Flex reconstruction examples can be found within the packages at the following locations:

**MATLAB:**
- `<MATLAB_TOP>/Examples/RunFlex.m`
- `<MATLAB_TOP>/Scripts/FlexRecon.m`

**C++:**
- `<INSTALL_ROOT>/recon/Orchestra/Flex/Rehearsal/RunFlex.cpp`

### 3D ARC Reconstruction Enhancements

Updates were made to the 3D ARC reconstruction script in MATLAB to better handle functionality such as slice zipping, reverse slice ordering, and when no kacq file is available. The main change involves pulling the data out of the Pfile in acquisition/time order instead of final geometric order. Because 3D ARC Pfiles are a special case, this is the only example where such logic is required.

**Updated file:**
- `<MATLAB_TOP>/Scripts/Arc3DRecon.m`

--------------------------------------------------------------------------------

## 1.2-115

### Mac OS X support for C++ and MATLAB development

### True Color DICOM support in C++ and MATLAB

A True Color image explicitly specifies a Red, Green, and Blue (RGB) value for each pixel. The values are 8-bit unsigned values and can thus range from 0-255. Each of the three color channels must be specified for each pixel. For example, a purely Red pixel would be `[255, 0, 0]` and a purely Blue pixel would be `[0, 0, 255]`.

In MATLAB, True Color images are stored as `[Rows, Cols, 3]` where the last dimension represents the color channel. In C++, the type `MDArray::Array<unsigned char, 3>` is used to represent a True Color image. This is typedef’d to `GEDicom::RGBImage` and used in many APIs. Either the first or last dimension can represent the color channel; however, only one of these dimensions should have an extent of 3.

**Example:** Read a grayscale DICOM image in MATLAB and use a built-in color map to write it as a color image.

```matlab
% Read reference image
img = GERecon('Dicom.Read', 'gray_ref.dcm');
% Convert to Color
gray = mat2gray(img); % Scale from 0 to 1
scaled = uint16(65536 * gray); % Use dynamic 16-bit range
cmp = jet(65536); % Create 'JET' colormap
img_color = uint8(255 * ind2rgb(scaled, cmp)); % Scale to 8-bit range
% Write color image using original image as reference
GERecon('Dicom.Write', 'color.dcm', img_color, 'gray_ref.dcm', 100);
```

### DICOM API Enhancements

The DICOM classes now have a more complete set of APIs for both reading and writing images. All classes have getter/setter functions for standard and GE private tags. See `<INSTALL_ROOT>/recon/Dicom/MR/Image.h` and the corresponding module classes for available functions.

Example: Use the `MR::ImageModule` class to specify parameters for each image. Code is inserted into the Cartesian2D example Rehearsal. Don’t forget to include the `ImageModule` header (`#include <Dicom/MR/ImageModule.h>`).

```cpp
// Create DICOM image
const GEDicom::MR::ImagePointer dicom = dicomSeries.NewImage(finalImage, imageNumber, imageCorners);

// Specify temporal position (phase) information and current echo number.
dicom->ImageModule()->TemporalPositions(GEDicom::ToString(numPhases));
dicom->ImageModule()->TemporalPosIdentifier(GEDicom::ToString(currentPhase));
dicom->ImageModule()->EchoNumber(GEDicom::ToString(currentEcho));
```

### Orientation and Corner Point API Improvement

In the previous release, a very subtle (and unfortunately common) programming error led to DICOM images having incorrect geometry annotation. New C++ APIs and classes were introduced to prevent the programming error from occurring. Also, one API was removed. This will force any code developed from previous examples to be updated.

Below is a snippet of the **OLD** key code used to orient an image matrix and create a DICOM image:

```cpp
// Get information for current slice
const SliceOrientation& sliceOrientation = pfile->Orientation(currentSlice);
const SliceCorners& sliceCorners = pfile->Corners(currentSlice);
// Other steps...
// Rotate/Transpose image and corner points accordingly
SliceCorners imageCorners = sliceCorners;
FloatMatrix finalImage = RotateTranspose::Apply<float>(imageCorners,
    magnitudeImage, sliceOrientation.RotationType(), sliceOrientation.TransposeType());
// Other steps...
// Create DICOM image
const GEDicom::MR::ImagePointer dicom = dicomSeries.NewImage(finalImage,
    imageNumber, sliceOrientation, imageCorners);
```

In this old example, a `SliceCorners` object extracted from the Pfile needed to be rotated/transposed or oriented in the same fashion as the final image data. This is done to get the image and the corner point in the desired presentation display. However, if the `SliceCorners` were not oriented, the DICOM image creation had no way of knowing. If this subtle, but very important step was missed, the annotation was likely to be very wrong!

To solve the issue, the Rotate/Transpose of the corners was replaced by a new class, `ImageCorners`. This object takes the original `SliceCorners` and a `SliceOrientation` object to internally calculate the desired presentation corner points for the output image. The DICOM creation functionality now requires an `ImageCorners` object, which is guaranteed to have the information required. This is a simple example of leveraging the type-safe nature of the language, as opposed to using a `SliceCorners` object that at any given time may represent very different things.

The explicit Rotate/Transpose of the `SliceCorners` was removed, and the DICOM image creation API was updated to use the new `ImageCorners` class. **The `RotateTranspose::Apply()` function with `SliceCorners` has been completely removed.** Existing scripts will need to be updated. Below is the **NEW** code and the correct way to orient images and create DICOMs:

```cpp
// Get information for current slice
const SliceOrientation& sliceOrientation = pfile->Orientation(currentSlice);
const SliceCorners& sliceCorners = pfile->Corners(currentSlice);
// Other steps...
// Rotate/Transpose output image
FloatMatrix finalImage = RotateTranspose::Apply<float>(magnitudeImage,
    sliceOrientation.RotationType(), sliceOrientation.TransposeType());
// Other steps...
// Create DICOM image
const ImageCorners imageCorners(sliceCorners, sliceOrientation);
const GEDicom::MR::ImagePointer dicom = dicomSeries.NewImage(finalImage,
    imageNumber, imageCorners);
```

Note that the MATLAB APIs did not change, but the underlying implementation of the structures did. The `Orient` function with Pfile corner points will still run, but will return a wrapper of the `ImageCorners` class. The MATLAB DICOM function will accept either the new `ImageCorner` structure or the existing `SliceCorners`. This was provided for backward compatibility and should not affect any existing MATLAB scripts.

--------------------------------------------------------------------------------

## 1.1-88

- Initial public release
- C++ development in Linux
- MATLAB development in Linux and Windows

--------------------------------------------------------------------------------
