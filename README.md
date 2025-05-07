QPM Cell Analysis Pipeline
==========================

This repository provides a MATLAB and Python-based pipeline for processing and analyzing 
Quantitative Phase Microscopy (QPM) images of biological cells. The pipeline consists of 
three main steps: phase unwrapping, cell cropping, and feature extraction.

Contents
--------

- step1_phase_unwrapping.m         : MATLAB script for phase unwrapping and thickness calculation
- step2_cell_cropping.m            : MATLAB script for cropping regions of interest (ROIs)
- step3_feature_extraction.m       : MATLAB script for feature extraction from cropped images
- step3_feature_extraction_py.py   : Python script (partial reimplementation of step 3)
- LIScent.tiff                     : Reference image for user input
- ResultsXXX.csv / deadXXX.csv     : Bounding box metadata for live/dead cells
- resultsXXX.mat                   : QPM phase map data
- README.txt                       : This file

Overview of Steps
-----------------

Step 1: Phase Unwrapping (MATLAB)
  - Loads QPM `.tif` image stacks
  - Applies phase-shifting and unwrapping using Goldstein's algorithm
  - Outputs OPD maps, thickness maps, and normalized grayscale images

Step 2: Cell Cropping (MATLAB)
  - Loads center image (`LIScent.tiff`) and waits for a user click
  - Uses bounding boxes from `.csv` to extract regions for each cell
  - Saves cropped grayscale images for live and dead cells

Step 3: Feature Extraction (MATLAB / Python)
  - Loads cropped images and applies image segmentation
  - Extracts features like volume, surface area, roundness, dry mass
  - Performs statistical analysis and saves plots

Python Version
--------------
- The Python script (`step3_feature_extraction_py.py`) reimplements a subset of MATLAB Step 3
- Includes polynomial background subtraction, adaptive thresholding, and morphology operations
- Uses OpenCV, scikit-image, and matplotlib for visualization
- Feature extraction includes volume, area, and dry mass (others are outlined but not calculated)


Acknowledgments
---------------
Goldstein phase unwrapping algorithm:
C. Smith. Goldsteinunwrap2d_r1 (2023). MATLAB Central File Exchange.
https://www.mathworks.com/matlabcentral/fileexchange/29497-goldsteinunwrap2d_r1

Please cite this resource if using the unwrapping method in your publication.

License
-------
This project is provided for academic and non-commercial use only.
