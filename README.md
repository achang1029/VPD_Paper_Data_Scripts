# VPD_Paper_Data_Scripts

This repository contains all code and scripts used for the Manuscript titled "Rapid Increases in Heat-Driven Extreme Moisture1
Demand in the Southern Amazon".

Descriptions on each of the files are listed below:

**Data Regridding.ipynb**

This notebook was used to regrid any CMIP6 model data downloaded into a common resolution of 1 degree by 1 degree for intercomparison purposes.

**ERA5_land_sea_mask.nc**

This file is a land/sea mask to filter out land and sea areas in ERA5 data.

**Flow Analog Data Plotting.ipynb**

This Jupyter notebook is used to plot all figures (main and supplementary) in the manuscript.

**Flow Analog Training and Data Generation.ipynb**

This Jupyter notebook is used to create the multivariate constructed flow analog data from the ERA5 dataset for VPD, es, and ea.

**Variable Creation.ipynb**

This Jupyter notebook is used solely to concatenate datasets containing multiple files with different years for the same variable into single files containing all of the years for the associated variables.

**bias_correction_daily.m**

This Matlab file is the script used to bias-correct any CMIP6 projections data using Quantile Delta Mapping (QDM).

**d2m_to_hurs.m**

This Matlab file is a script that converts dew point to relative humidity.

**d2m_to_svp.m**

This Matlab file is a script that converts dew point to actual vapor pressure (ea).

**hurs_to_d2m.m**

This Matlab file is a script that converts relative humidity to dew point.

**qdm_z.m**

This Matlab file is a script that contains the Quantime Delta Mapping (QDM) algorithm.
