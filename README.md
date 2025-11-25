# SiO₂–Au Nanorod LSPR Analysis

This repository contains a single Python script that analyses the localised surface plasmon resonance (LSPR) of silica‐coated gold nanorods (AuNR@SiO₂) in water/glycerol mixtures.

The script:

- Reads UV–Vis spectra from Excel files for a series of refractive indices (water/glycerol mixtures).
- Finds the LSPR peak positions using a Gaussian fit.
- Optionally loads previously fitted peak positions from a separate Excel file.
- Fits a simple dielectric model to extract:
  - Nanorod aspect ratio (`AR`)
  - Effective AuNR volume fraction (`f`)
- Plots the experimental data and model for different CTAB concentrations.

---

## 1. File and data structure

The script expects **five** Excel files in the same directory:

```text
porosity -no CTAB.xlsx
porosity -5 mM CTAB.xlsx
porosity -10 mM CTAB.xlsx
porosity -15 mM CTAB.xlsx
All LSPR peak positions 21-7-22.xlsx
You can change these names in the __main__ block if needed.

1.1 UV–Vis files
Files:

porosity -no CTAB.xlsx

porosity -5 mM CTAB.xlsx

porosity -10 mM CTAB.xlsx

porosity -15 mM CTAB.xlsx

Expected format for each file:

Column 0: Wavelength (nm)
Columns 1–7: UV–Vis spectra (absorbance) data for different water/glycerol mixtures
