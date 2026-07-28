# Bates2026
A repressive regulatory cascade shapes temporal patterning of activity-regulated gene expression in a defined sensory neuron type
Code repository

This repository contains MATLAB image analysis code used to quantify GFP reporter expression and sub-cellular localization. Code has only been verified to work 
on MATLAB2023b and will not work on older versions (pre addition of volshow 2017a) or more recent versions which have updated UI figure drawing/viewer syntax. 
Functions included in this repository are as follows:

  measureCTCF: Intensity threshold-based ROI detection and selection from a directory of .ims 3D volumes, returns per ROI CTCF calculation 
  (mean(ROI_intensity) x ROI_area) - (Background_intensity x ROI_area) and uncorrected mean as Nx2 array.
  
  segmentCells3D: Intensity threshold-based 3D ROI detection and selection on 3D volumes. Returns selected 3D ROIs in individual 3D arrays.
  
  segmentCells3Dbatch: run segmentCells3D on a directory on 3D .ims volumes. Returns a cell array of selected 3D ROIs.
  
  annotate_cells: Rotate, project, and annotate nuclear centers of volumes output by segmentCells3D. Outputs struct array with fields:
                  .R            3×3 rotation matrix
                  .proj_axis    1=XY, 2=XZ, 3=YZ
                  .nuc_center   [row col] in chosen projection
                  .projection   2D double mean projection
                  
  average_cell_intensity: Min/max normalize, draw nuclear ROIs, and calculate a variety of signal distribution metrics on projections output from annotate_cells. Also,
  aligns all input projections over their annotated nuclei. Outputs a struct with fields containing the many different metrics.
  
  replay_figure: plotting function to visualize results contained in the average_cell_intensity output structure and produce mean projections of aligned cells.
  Produces the following plots

  1. Cumulative histogram
  <img width="560" height="420" alt="image" src="https://github.com/user-attachments/assets/af9fc8b5-a21a-43db-a340-51425d491e8c" />
  
  
  2. Radial sum
  <img width="560" height="420" alt="image" src="https://github.com/user-attachments/assets/739fffe9-8861-4ca6-a632-20caf9d0584e" />



  3.Radial mean
  <img width="560" height="420" alt="image" src="https://github.com/user-attachments/assets/0508a7a1-b55f-49e5-b683-7a5377a5efab" />


  4. Scalar metrics#1
  <img width="560" height="420" alt="image" src="https://github.com/user-attachments/assets/f5702347-f43e-4216-924e-f89ef4051a79" />


  5+6. Line scan and column mean profiles
  <img width="560" height="420" alt="image" src="https://github.com/user-attachments/assets/958b434f-caaa-4bbc-bd2a-aaf6ee8abbd7" />


  7-9. Metric distributions
  <img width="560" height="420" alt="image" src="https://github.com/user-attachments/assets/905be30f-fb4c-48de-98d5-200aeb47d347" />


  10. Scalar metrics#2
  <img width="560" height="420" alt="image" src="https://github.com/user-attachments/assets/4677f525-9a5b-4ce4-924d-acc972b99e8c" />






