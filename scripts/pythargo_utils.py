#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: aoetjens
"""


# Import modules
import glob
import xarray as xr
import numpy as np
import pandas as pd
from pandas.plotting import register_matplotlib_converters
import seaborn as sns
import geopandas as gpd
import cmocean
import gsw
# MATPLOTLIB
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import cm, ticker
from matplotlib import colors as clrs 
from matplotlib.ticker import LogFormatter 
import matplotlib.path as mpath
from matplotlib.patches import Patch, Rectangle
from matplotlib.colors import ListedColormap


# CARTOPY
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature

from shapely.ops import unary_union
from shapely import geometry
#SCIPY
from scipy.interpolate import griddata
import scipy.io
from scipy.io import loadmat
import scipy.signal as signal
from scipy import interpolate
from scipy.optimize import curve_fit
import scipy.stats as stats
from scipy.stats import linregress, ks_2samp
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score

# ADDONS
import copy
import dask
import calendar
from datetime import datetime
import time
from netCDF4 import Dataset as netcdf_dataset
import h5py
# DASK
import dask
import dask.array as da
import dask.dataframe as dd
from dask import delayed
from dask.diagnostics import ProgressBar  # Import ProgressBar


# other python scripts
from pythargo_CbPMargo import  daylength, AustinPetzold_1986, cbpm_bgcfloats
import sys 
import os
path_main = '/home/581/ao2582/SONIF/data/ARGO/'  # Location of this file and where the subdir live
sys.path.append(os.path.abspath(path_main))
path_g = '/../../../g/data/jk72/ao2582/'

# Import files
#aviso_file= glob.glob(path_g + 'seasonal_clim/madt_u*/*.nc')[0]
#ugos = xr.open_dataset(aviso_file)['ugos'].data[0,::3,::3]

#aviso_file= glob.glob(pat_g + 'AVISO/climatology/global/delayed-time/seasonal_clim/madt_v*/*.nc')[0]
#vgos = xr.open_dataset(aviso_file)['vgos'].data[0,::3,::3]

#lat_geo, lon_geo = xr.open_dataset(aviso_file)['latitude'].data[::3],xr.open_dataset(aviso_file)['longitude'].data[::3]

#%%
region_dict = {'SOUTHGEORGIA' : {'region' : [-60,0,-65,-40], 'contour' : 0.45, 'geopoint' : (-37,-54), 'file' : 'SONIF_SOUTHGEORGIA_bo_qcv3.nc'},
                'SOAZ' : {'region' :[-120,30,-80,-40], 'contour' : 0, 'geopoint' : (-60,-55)},
                'KERGUELEN' : {'region' : [65, 100,-60,-40], 'contour' : 0.2, 'geopoint' : (70,-49), 'file' : 'SONIF_KERGUELEN_bo_qcv3.nc'},
                'KERGUELEN_c' : {'region' : [0,180,-80,-40], 'contour' : 0, 'geopoint' : (-37,-54)},
                'BALLENY' :{'region' : [140, 175, -70, -55], 'contour' :0.3 ,'geopoint' : (163,-67), 'file' : 'SONIF_BALLENY_bo_qcv3.nc'},
               'SOUTHERNOCEAN' :{'region' : [-180, 180, -90, -40], 'contour' :0 ,'geopoint' : (163,-67)},
               'SOPACIFIC' :{'region' : [-180, -60, -90, -40], 'contour' :0 ,'geopoint' : (-120,-50)},
               'SOATLANTIC' :{'region' : [-60, 40, -90, -40], 'contour' :0 ,'geopoint' : (0,-50)},
               'SOINDIAN' :{'region' : [40, 180, -90, -40], 'contour' :0 ,'geopoint' : (120,-50)},
               'test' :{'region' : [0, 2, -50, -40], 'contour' :0 ,'geopoint' : (5,-45)},
               'BOUVET' :{'region' : [-5,30,-70,-50], 'contour' :0.3 ,'geopoint' : (3.3,-54.4)},
               'SICE' :{'region' : [-180,180,-80,-45], 'contour' :0.3 ,'geopoint' : (3.3,-60.4)}
                }

#print(path_g)
#print(glob.glob(path_g + 'SONIF_*'))
aviso_file= path_g + 'AVISO/dt_global_allsat_eke_y1993_2021_m04_06.nc'
eke_clim = xr.open_dataset(aviso_file)['eke']

#%% 1. Read in data

zonal_dict = {
    'saz' : {'label' : 'Sub Antarctic Zone', 'color' : 'blue'},
    'pfz' : {'label' : 'Polar Frontal Zone', 'color' : 'orange'},
    'az' : {'label' : 'Antarctic Zone', 'color' : 'cyan'},
    'sz' : {'label' : 'Southern Zone', 'color' : 'purple'},
    'nz' : {'label' : 'No Zone', 'color' : 'grey'}
    }


# basic functions to work with Argo data or satellite products

def read_sat(fosor):
    print('read sat')
    print(region_dict[fosor])
    sat_file = path_g +"AQUA_MODIS/CHL/season_clim/AQUA_MODIS.20021221_20230320.L3m.SCWI.CHL.x_chlor_a.nc"
    region = region_dict[fosor]['region']
  #  sat_file = glob.glob(path_sat + '*' + str(year) + str(month) + '*.nc')[0]
    dq = xr.open_dataset(sat_file)
    ds = dq.where((dq.lon > region[0])&(dq.lon < region[1])&(dq.lat > region[2])&(dq.lat < region[3]), drop=True)
    #Check satellite data
    chla_oc = ds.chlor_a.coarsen(lat=5, lon=5).mean()
    dq.close()
    ds.close()
    return chla_oc

def chl_box(fosor, chl_contour = True):
    from shapely.geometry import Point
    print(fosor)
    region = region_dict[fosor]['region']
    if chl_contour:
        contour = region_dict[fosor]['contour']
        loc = geometry.Point(region_dict[fosor]['geopoint']) #geometry.Point(70,-49)
        chla_oc = read_sat(fosor)
        ax = plt.axes(projection=ccrs.Mercator(), frameon=False)
        #ax.set_extent(extent, ccrs.PlateCarree())
        contours = ax.contour(chla_oc.lon.data,chla_oc.lat.data,chla_oc.data,[contour], transform= ccrs.PlateCarree())
       # ax.plot(-37,-54,'ro', transform=ccrs.PlateCarree())
        for k in range(len(contours.allsegs[0])):
            #print(k)
               #         contour level , first polygon from that level
            CS = contours.allsegs[0][k].T 
        #poly = Polygon([(i[0], i[1]) for i in CS])
            if len(CS[1,:]) > 4:
                poly = geometry.Polygon([(CS[0, i], CS[1, i]) for i in range(len(CS[1,:]))])
                boundingbox = gpd.GeoSeries(poly)
             #   boundingbox.plot()
                if boundingbox.contains(loc)[0]: 
                    
                    print('contour polygon found, k = ',k)
                 #   print(k)/'
                    return boundingbox, poly
    else:
        poly = geometry.Polygon([Point(region[0],region[2]),Point(region[0],region[3]),Point(region[1],region[3]),Point(region[1],region[2])])
        boundingbox = gpd.GeoSeries(poly)
     #   boundingbox.plot()
        return boundingbox, poly
print('chl_box  --- imported')
def argo(wmo, path):
    file = glob.glob(path +str(wmo) +'_Sprof.nc')[0]
    dq = xr.open_dataset(file)
    return dq

def create_boundingbox(region): 
    extent = region
    p1, p2, p3, p4 = (extent[0], extent[2]), (extent[1], extent[2]),(extent[1], extent[3]),(extent[0], extent[3])
    '''provide four corner points as (lon,lat), 
       the order is bottom-left(p1), bottom-right(p2), top-right(p3),top-left(p4) on a map'''
    p1 = geometry.Point(p1)
    p2 = geometry.Point(p2)
    p3 = geometry.Point(p3)
    p4 = geometry.Point(p4)
    pointList = [p1, p2, p3, p4, p1]
    boundingbox = geometry.Polygon([[p.x, p.y] for p in pointList])
    boundingbox = gpd.GeoSeries(unary_union(boundingbox))
    return boundingbox

def nc_get_nearest(lon, lat, nc_file):
    # extract parent netcdf and lat lon variables
    dataset = netcdf_dataset(nc_file)
    nc_var = dataset.variables['Band1']
    dataset.close()
    nc = nc_var.group()
    dims = nc_var.dimensions
    for i, name in enumerate(dims):
        if 'la' in name:
            latdi, latname = i, name
        if 'lo' in name:
            londi, lonname = i, name
    lat_var, lon_var = nc.variables[latname], nc.variables[lonname]

    # convert lon to +- 180
    if lon_var[:].max() > 180.:
        lon_var = lon_var[:]  # extract values
        lon_var[lon_var > 180] = lon_var[lon_var > 180] - 360  # correct >180
        index = abs(lon_var[1:] - lon_var[:-1]).argmax()  # determine crossover point
        lon_var = np.roll(lon_var, index + 1)
        ax = dims.index(lonname)  # dimension index
        nc_var = nc_var[:]  # extract variable data
        nc_var = np.roll(nc_var, index + 1, axis=ax)


    # find nearest lon lat indices
    lati = abs(lat - lat_var[:]).argmin()
    loni = abs(lon - lon_var[:]).argmin()

    slc = [slice(None) for d in dims]
    slc[latdi] = slice(lati, lati + 1)
    slc[londi] = slice(loni, loni + 1)
    slc = tuple(slc)
    

    return nc_var[slc].squeeze()








