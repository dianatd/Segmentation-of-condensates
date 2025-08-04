# Condensate Analysis: Image Segmentation and Contact Angle Measurement

This repository contains a MATLAB script designed to analyze two-channel fluorescence microscopy images of phase-separated condensates. The code automates the process of identifying condensates, quantifying their properties, and measuring the contact angles between phases in de-mixed systems.

---

##  Overview

The main analysis script performs the following steps:

### 1. User Input
Prompts the user to:
- Select a data folder
- Enter a basename string (common to all image filenames)
- Accept or reject contact angle fitting

### 2. Image Segmentation
Automatically detects condensate boundaries in each two-channel `.tif` image using thresholding and Gaussian filtering.

### 3. Metrics Calculation
For each identified condensate, it calculates:
- **Correlation coefficients** (Pearson’s) between the two channels
- **Phase area** for each channel (in µm²)
- **Fluorescence partitioning coefficients** (ρ, rho) to quantify mixing 

### 4. Contact Angle Analysis
For de-mixed condensates:
- Fits circles to phase boundaries and their interface
- Calculates the two contact angles
- Prompts the user to validate fits for accuracy

### 5. Data Export
Saves all quantitative results to `.csv` files and exports visualization images of selected condensates as `.png`.

---

## Requirements

### Software
- MATLAB (tested on R2021a or later)
- Image Processing Toolbox

### Dependencies
Ensure the following helper functions are present in the same folder or in your MATLAB path:
- `ContactAnglesDT.m`: Computes contact angles from segmented images (included)
- `CircleFitByPratt.m`: Pratt-style circle fitting (included)

---

##  How to Use

### Setup
1. Place all required `.m` files in one folder:
   - `Image_segmentation_for_condensates.m`
   - `ContactAnglesDT.m`
   - `CircleFitByPratt.m`
2. Organize your two-channel `.tif` images in a single directory.

### Run the Analysis
1. Open MATLAB
2. Run: Image_segmentation_for_condensates.m

Licensed under [CC BY-NC 4.0](https://creativecommons.org/licenses/by-nc/4.0/) – for academic use only.