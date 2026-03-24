# Bioimpedance Signal Processing and TEER Validation Suite

This repository contains the MATLAB framework developed for my Master's Thesis in Biomedical Engineering at the **University of Florence** (Graduated 110/110 Cum Laude with Academic Mention).

## 🎯 Project Scope
The project involves the design, characterization, and validation of an in-vitro monitoring system for biological tissues. Using the **MAX30009 Analog Front-End (AFE)**, the system performs precise **TEER (Trans-Epithelial/Endothelial Electrical Resistance)** measurements via low-frequency microcurrent stimulation (MCS).

## 🚀 Key Technical Features
- **Hardware Optimization**: Automated logic to configure PLL dividers and DAC parameters for the MAX30009, ensuring frequency precision and stimulus stability.
- **Metrological Validation**: A complete pipeline to convert raw ADC counts into physical Impedance (Ohm), validated against passive benchmarks and Gold Standard instruments (Millicell ERS-2).
- **Physical Modeling**: Simulation of current density distribution ($J$) based on actual 24-well plate and STX01 electrode geometries to ensure stimulus uniformity.
- **Signal Processing**: (Incoming) Advanced filtering and outlier removal using Median Absolute Deviation (MAD) for biological signal integrity.

## 📂 Repository Structure
- **`Stimulus_Frequency_Calculator.m`**: Computes hardware-level parameters (PLL, KDIV, OSR) for the desired frequency.
- **`Target_Current_Calculator.m`**: Optimizes voltage and resistance stages to achieve the target RMS current.
- **`Current_Settings.m`**: Core function for calculating effective voltage and resistance combinations (Dependency).
- **`Current_Density_Distribution.m`**: Geometric field simulation and visualization of the electrical stimulus.
- **`ADC_Counts_to_Ohm_Conversion.m`**: Automated processing of large CSV datasets for calibration and impedance conversion.

## 🎓 About the Author
**Gabriele Rellini** *Biomedical Engineer specialized in Biomechanics, Biomaterials, and Tissue Engineering.* University of Florence | Master's Degree with Honors.

Expertise: MATLAB, Bioimpedance, Signal Processing, R&D.

[![LinkedIn](https://img.shields.io/badge/LinkedIn-Connect-blue)](https://www.linkedin.com/in/gabrielerellini/)
