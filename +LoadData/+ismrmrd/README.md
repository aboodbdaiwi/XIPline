# Generic_Data_Format_Tools
 Tools to write hyperpolarized 129Xe MRI data acquired according to 129Xe MRI Clinical Trials Consortium Recommendations to MRD format. Currently, Siemens and GE data are supported.

## Usage
All functions are called from the master function xemri_mrd.m. Once run, this function will prompt the user to enter the participant ID and locate the raw data files that need to be converted. Converted MRD files will be written to the same path where the xenon MRI data is stored. 

## Dependencies
In order to write to MRD files, you'll need the ismrmrd repository: https://github.com/ismrmrd/ismrmrd
Add the ismrmrd/matlab folder to your matlab path.

## File Contents
Detailed specifications for the contents of each generated MRD file are contained in MRD_Spec.md

## Site Specific Edits
### In +GE/dissolved_to_ismrmrd.m
- Line 32: File containing frequency offset
- Lines 96-98: Scanner make, model, and institution
### In +GE/calibration_to_ismrmrd.m
- Line 20: File containing frequency offset
- Lines 66-68: Scanner make, model, and institution 
