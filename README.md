# QPM Cell Analysis Pipeline

This repository contains a MATLAB-based pipeline for processing and analyzing **Quantitative Phase Microscopy (QPM)** images of biological cells. It performs phase unwrapping, cell cropping, and quantitative feature extraction from live and dead cells.

---

## 📋 Pipeline Overview

### Step 1: Phase Unwrapping (`step1_phase_unwrapping.m`)
- Applies phase-shifting interferometry on `.tif` image stacks.
- Uses Goldstein's phase unwrapping algorithm to retrieve quantitative phase maps.
- Calculates Optical Path Difference (OPD) and cell thickness.
- Outputs:
  - Normalized grayscale images
  - Thickness maps
  - `.mat` files for OPD and unwrapped phases

### Step 2: Cell Cropping (`step2_cell_cropping.m`)
- Loads QPM phase results and metadata from `.csv`.
- Allows the user to select a reference point on the image.
- Crops regions of interest (ROIs) for both live and dead cells.
- Outputs:
  - Individual cropped grayscale `.tif` files
  - `.mat` files with ROI and distance information

### Step 3: Feature Extraction (`step3_feature_extraction.m`)
- Processes cropped cell images and masks.
- Extracts multiple features such as:
  - Area, volume, surface area
  - Circularity, roundness, flatness
  - Dry mass and related shape descriptors
- Performs statistical comparisons and generates plots.
- Outputs:
  - `.mat` feature files
  - `.bmp` plots
  - Sign-rank test results in the console

---

## 📦 Output Files

| File                          | Description                         |
|-------------------------------|-------------------------------------|
| `resultsXXX.mat`              | Unwrapped phase and thickness       |
| `OPDXXX.mat`, `phaseXXX.mat`  | Optical path and raw phase data     |
| `ImcropmatrixXXX.mat`         | Bounding boxes for live cells       |
| `deadcropmatrixXXX.mat`       | Bounding boxes for dead cells       |
| `finalmaskXXX.mat`            | Final binarized masks               |
| `.bmp`                        | Boxplots per feature                |

---

## 📖 Citation and Acknowledgments

**Goldstein Phase Unwrapping Algorithm:**

This project uses the Goldstein 2D phase unwrapping implementation by C. Smith:

> Smith, C. (2023). *Goldsteinunwrap2d_r1* [MATLAB Central File Exchange].  
> https://www.mathworks.com/matlabcentral/fileexchange/29497-goldsteinunwrap2d_r1

Please cite this tool if you use this pipeline for publications involving phase unwrapping.


