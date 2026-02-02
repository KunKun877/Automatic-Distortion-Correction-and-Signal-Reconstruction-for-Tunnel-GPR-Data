# Automatic Distortion Correction and Signal Reconstruction for Tunnel GPR Data

[![MATLAB](https://img.shields.io/badge/Made_with-MATLAB-orange.svg)](https://www.mathworks.com/products/matlab.html)
[![License](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

This repository contains the source code and demonstration data for the paper:  
**"Automatic Distortion Correction and Signal Reconstruction in GPR Data in Tunnel Lining Surveys"**.

## 📖 Introduction

Ground-Penetrating Radar (GPR) data acquired in time-measurement mode often suffers from two major quality issues:
1.  **Geometric Distortion:** Caused by non-uniform scanning speeds.
2.  **Random Data Loss:** Caused by wireless transmission instability.

This project implements an integrated processing workflow to automatically address these challenges. It combines a statistical velocity analysis module for distortion correction with an **Improved Adaptive Rank-Reduction Interpolation** algorithm (based on DMSSA) to recover missing data with high fidelity.

## ✨ Key Features

* **Automatic Distortion Correction:**
    * Estimates the standard scanning velocity via Probability Density Function (PDF) mode analysis.
    * Performs trace position remapping using "normalized velocity" as the resampling factor.
* **Adaptive Signal Reconstruction:**
    * **DMSSA Framework:** Utilizes Damped Multi-channel Singular Spectrum Analysis for simultaneous denoising and interpolation.
    * **Adaptive Rank Selection:** Selects the optimal damping coefficient ($K \approx 4$) based on Random Matrix Theory (4$\sigma$ rule) to distinguish signal subspaces from the noise floor.
    * **POCS Solver:** Incorporates a Projection Onto Convex Sets (POCS) iterative loop to strictly preserve valid observed data.

