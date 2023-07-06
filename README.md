# Vallat, Shah and Walker (2023).

This repository contains Python code to reproduce the results and figures of:

- Vallat, Shah and Walker(2023). Coordinated human sleeping brainwaves map peripheral body glucose homeostasis. _Cell Reports Medicine_.

**Data and Code Availability**: The CFS dataset can be obtained, with the appropriate permissions, and for non-commercial use, at https://sleepdata.org/datasets/cfs. The MESA dataset can be obtained, with the appropriate permissions, through the BioLINCC repository at: https://biolincc.nhlbi.nih.gov/studies/mesa/. In addition to the public access repository, interested investigators may also access the data through the MESA Coordinating Center at the University of Washington. Use of the data via this mechanism is overseen by standard MESA policies and procedures, which assure that participant consent is honored. Additional information can be found at: https://www.mesa-nhlbi.org/.

**Lead Contact**: Further information and requests for resources should be directed to the lead contact, Dr. Raphael Vallat (raphaelvallat9@gmail.com)

*********

## Code structure

### 01_coupling_*

Calculate the SO-spindles (sigma) coupling from the PSG recording.

### 01_extract_HRV.ipynb

Calculate overnight HRV metrics from the PSG recording.

### 01_extract_sleepstats.ipynb

Calculates the sleep statistics from the hypnogram.

### 01_extract_spectral_power.ipynb

Calculates the EEG spectral power from the PSG recording.

### 02_data_correlation_*

These notebooks combine the different data streams into a single dataframe. Specifically, the following steps are applied:

1. Concatenate the demographics & health, coupling, spectral power and HRV data
2. Remove participants <15 yrs old (not applicable in MESA)
3. Preprocessing: outlier removal and data transformation
4. Calculate the HOMA variables
5. Calculate and plot the correlation between coupling and glucose, adjusted for age (partial correlation)
6. Export `df_concat_R_1sec_[CFS/MESA].csv`, the main spreadsheet that will be used for subsequent analyses

### 03_regression_*

Main statistical analysis pipeline, written in R.

### Table_Demographics.ipynb

Generate tables for CFS and MESA datasets with descriptive statistics.

*********

## Figures

The Python/R code to generate the figures of the paper can be found in:

- Figure 1A: `Fig_ERGCPAC.ipynb`
- Figure 1B: `02_data_correlation_CFS.ipynb`
- Figure 1C-D: `Fig_SO-coupling.ipynb`
- Figure 2: `02_data_correlation_CFS.ipynb`
- Figure 3: `03_regression_MESA.R`
- Figure 4: `02_data_correlation_CFS.ipynb`
- Figure 5: `Fig_TopSleepAdjusted.ipynb`