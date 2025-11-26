# An improved model of particle attenuation reduces estimates of Southern Ocean carbon transfer efficiency
(under revision)
Code for the analysis and figure production

## Python code for FIgure production

**Datasets:**
- AVISO sea_surface_height_above_geoid
"Climatological Absolute Dynamic Topography 1993 to 2021/01"

- BGC_argo_processed.nc (https://doi.org/10.25959/4t8m-a057)
- ds_cbpm_TE.nc (https://doi.org/10.25959/4t8m-a057)
- stats_df.pkl (in data folder)



## Data processing (see paper for detail)

**BGC-Argo floats**
- Data download using the Matlab OneArgo toolbox (https://github.com/NOAA-PMEL/OneArgo-Mat)
- Data QC and interpolation on 1 dbar grid, smoothing for bbp and chla, SO specific correctionf actor for chl

**Remote Sensing**
- AQUA-MODIS PAR monthly climatology
- surface carbon from the carbon productivity website based on AQUA-MODIS backscatter


## Computation of bgc parameter
- dynamic height and sea surface height for zones
- POC from Johnson and Koestner
- euphotic zone
- NPP from Arteaga

## Broken power law
- least-square fit
- monte carlo simulation for uncertainty estimate
- AIC for statistical test
- Spatial modelt preference, regridding on map
- Seasonal cycle of parameter, surface  averaged, monthly bins
- Linear regression for correlation test with temperature and NPP
- R2, p-values,...

