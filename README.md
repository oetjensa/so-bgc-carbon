# An improved model of particle attenuation reduces estimates of Southern Ocean carbon transfer efficiency
Code for the analysis and figure production of the paper:
(submitted to \textit{communications earth & envrionment})

## 1: Data processing

**BGC-Argo floats**
- Data download using the Matlab OneArgo toolbox (https://github.com/NOAA-PMEL/OneArgo-Mat)
- Data QC and interpolation on 1 dbar grid, smoothing for bbp and chla, SO specific correctionf actor for chl

**Remote Sensing**
- AQUA-MODIS PAR monthly climatology
- surface carbon from the carbon productivity website based on AQUA-MODIS backscatter

  

## 2: Computation of bgc parameter
- dynamic height and sea surface height for zones
- POC from Johnson and Koestner
- euphotic zone
- NPP from Arteaga

## 3: Broken power law
- least-square fit
- monte carlo simulation for uncertainty estimate
- AIC for statistical test
- Spatial modelt preference, regridding on map
- Seasonal cycle of parameter, surface  averaged, monthly bins
- Linear regression for correlation test with temperature and NPP
- R2, p-values,...

  ## 4: Figures
  - Figure 1
  - Figure 2
  - Figure 4
  - Figure 5
  - Supplementary figures

