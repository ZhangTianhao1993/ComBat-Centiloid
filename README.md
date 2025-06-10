# ComBat-Centiloid

A MATLAB-based implementation for converting SUVR values from amyloid PET tracers into Centiloid (CL) format using a modified ComBat harmonization approach.

## Overview

This project provides a method to standardize amyloid PET imaging data across different tracers by transforming SUVR values into the Centiloid scale. It employs a modified version of the `ComBat` algorithm adapted for PET harmonization.

## Input Description

- **`TRACERSUVR`**  
  A matrix of size `n × 84`, where `n` is the number of amyloid PET images and each column represents a cortical brain region.  
  The 84 cortical regions correspond to those defined in the article and can be extracted using the atlas provided in the `ROI_FOR_SUVR` folder.

- **`TRACERmod`**  
  A matrix of size `n × 4` representing demographic and diagnostic information for each subject:  
  - Column 1: Age (in years)  
  - Column 2: Gender (`1` for male, `0` for female)  
  - Column 3: Diagnostic status (`1` if AD, `0` otherwise)  
  - Column 4: Diagnostic status (`1` if MCI, `0` otherwise)  

> ⚠️ Note: For stable performance of the ComBat algorithm, it is recommended that `n > 30` and that the dataset includes a balanced mix of cognitively normal (CN), mild cognitive impairment (MCI), and Alzheimer’s disease (AD) subjects.
> 
 You can find sample input data for testing in the TestData folder.

## Output Description

- **`HCL`**  
  A vector of harmonized Centiloid values for each subject.

- **`HRCL`**  
  A matrix of region-wise Centiloid values for each subject.

## Key Function

- `combat_centiloid.m`  
  Main function that takes in `TRACERSUVR` and `TRACERmod`, and outputs `HCL` and `HRCL`.

- `combat_to_center.m`  
  A modified version of `combat.m` from the [ComBatHarmonization MATLAB project](https://github.com/Jfortin1/ComBatHarmonization/tree/master/Matlab), adjusted to allow alignment with a specific center.

## Example Usage

```matlab
[HCL, HRCL] = combat_centiloid(AV45_ADNI_SUVR, AV45_ADNI_Demo);
```

## Citation
If you use this code in your research, please cite the following paper:

Zhang Tianhao, Nie Binbin, Shan Baoci.
ComBat-Centiloid: A Calibration-Free Method for Quantifying Centiloid Values in Amyloid PET Imaging.
Brain Science Advances.
## Author
Zhang Tianhao  
December 2, 2024  
thzhang@ihep.ac.cn
