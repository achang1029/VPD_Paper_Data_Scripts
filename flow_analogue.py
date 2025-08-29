import os, sys
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import scipy.stats as stats
import pandas as pd
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import datetime as dt
import healpy as hp
import astropy as ap
from netCDF4 import Dataset, MFDataset
import h5py
import xarray as xr
from mpl_toolkits.mplot3d import axes3d
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import matplotlib.colors as mcolors
import cftime
import matplotlib.lines as mlines
from sklearn.decomposition import PCA
from sklearn import preprocessing
from sklearn.linear_model import LinearRegression
from scipy.spatial.distance import cdist
import glob
import csv

os.environ['HDF5_USE_FILE_LOCKING']='FALSE'

prefix = '/data/rong4/Data/'
prefix2 = '/data/rong3/achang1029/'

z500_fname = 'ERA5/z500_1979-2020.nc'
z500file = Dataset(prefix + z500_fname) #1979-2020

z500 = z500file.variables['z500'][:]

era5_lon_array = z500file.variables['longitude'][:]
era5_lat_array = z500file.variables['latitude'][:]

vpd_fname = 'gridMET/vpd_*.nc'
vpdfile = xr.open_mfdataset(prefix2 + vpd_fname)

gridmet_lon_array = vpdfile['lon'][:]
gridmet_lat_array = vpdfile['lat'][:]

vpd = vpdfile['mean_vapor_pressure_deficit'][:]

# for f in files:
#     with Dataset(f) as data:
#         var1 = data.variables['variable_name1'][:]
#         var2 = data.variables['variable_name2'][:]
#         # combine the data as needed
#         result = var1 + var2
        
# for i in range(1980, 2023):
#     print(i)
#     vpd_fname = 'gridMET/vpd_' + str(i) + '.nc'
#     vpdfile = Dataset(prefix2 + vpd_fname) #1979-2022
#     vpd = np.concatenate((vpd, vpdfile.variables['mean_vapor_pressure_deficit'][:]), axis=0)

def pairwise_distances(x, y):
    return cdist(x, y)
    
def check_feb_29(x):
    date = dt.datetime(1979,1,1) + dt.timedelta(x)
    return (date.month == 2 and date.day == 29)
    
#remove all feb 29 days to make time calculations easier
feb_29 = []
for i in range(z500.shape[0]):
    if check_feb_29(i):
        feb_29.append(i)
z500_feb29 = np.delete(z500, feb_29, axis=0)
print(feb_29)

feb_29 = []
for i in range(vpd.shape[0]):
    if check_feb_29(i):
        feb_29.append(i)
vpd_feb29 = vpd[np.arange(vpd.shape[0])[~np.isin(np.arange(vpd.shape[0]), feb_29)],:,:]
print(feb_29)

#cos(lat) weighted mean
z500_cos = z500_feb29 * np.expand_dims(np.cos(era5_lat_array * np.pi / 180) ** 2, axis=(0,2))
vpd_cos = vpd_feb29 * np.expand_dims(np.cos(gridmet_lat_array * np.pi / 180) ** 2, axis=(0,2))

#nMov-day moving avg 
nMov = 5
z500_xr = xr.DataArray(z500_cos, dims=['time', 'longitude', 'latitude'])
z500_nMov_avg = z500_xr.rolling(time=nMov, min_periods=1, center=True).mean().to_numpy()
vpd_nMov_avg = vpd_cos.rolling(day=nMov, min_periods=1, center=True).mean()

#Subtract the global spatial mean z500 from each day. i.e. subtract Jan 1 1979 spatial mean z500 from Jan 1 1979 values
z500_global_mean_removed = (z500_nMov_avg.T - np.mean(z500_nMov_avg, axis=(1, 2))).T
vpd_global_mean_removed = (vpd_nMov_avg.T - np.mean(vpd_nMov_avg, axis=(1, 2))).T

z500_reshaped = z500_global_mean_removed.reshape(42, 365, era5_lat_array.shape[0], era5_lon_array.shape[0])

z500_mean = np.mean(z500_reshaped, axis=0)
z500_std = np.std(z500_reshaped, axis=0)
z500_std_anom_reshaped = (z500_reshaped - z500_mean) / z500_std

z500_std_anom = z500_std_anom_reshaped.reshape(42*365, era5_lat_array.shape[0], era5_lon_array.shape[0])

#from https://stackoverflow.com/questions/43015638/xarray-reshape-data-split-dimension
multi_index = pd.MultiIndex.from_product((range(44), range(365)), names=("year", "day1"))
vpd_reshaped = vpd_global_mean_removed.assign_coords(day=multi_index).unstack("day").transpose("year", "day1", "lat", "lon")

vpd_mean = np.mean(vpd_reshaped, axis=0)
vpd_std = np.std(vpd_reshaped, axis=0)
vpd_std_anom_reshaped = (vpd_reshaped - vpd_mean) / vpd_std

vpd_std_anom = vpd_std_anom_reshaped.stack(day=("year","day1")).transpose("day", "lat", "lon").rename({"day":"time", "lat":"latitude", "lon":"longitude"})

# flow_analogue flow analogue search
# [DIST, IDOY, IYR, ID] = flowanalogue(X, Lon, Lat, indYrRef, nFA) 
# returns analogue day index and distance between analogue and observation.
# Lon and Lat are the longitudes and latitudes in a vector; X is the input
# circulation pattern variable array (e.g. geopotential) with the size of
# [N(lon), N(lat), N(time)]; indYrRef is the reference
# year index for searching analogues, e.g. 1:32; nFA is the number of flow
# analogue returned to output. Ouput analogue indices include ID, IDOY, 
# and IYR: IDOY and IYR are index of DOY (day of year) and year, and ID  
# equals IDOY + (IYR-1)*N(day-of-year). All output has the same size of 
# [N(day-of-year), N(year), nFA, N(ensemble)].
#  
# [DIST, IDOY, IYR, ID] = flowanalogue(..., OPTIONAL_ARGUMENT, ARG_VALUE)
# Optional arguments include:
#     nd: half window size for analogue search, window size is nd*2+1.
#       Default is 30;
#
#    distFunc: method for calculating distance, include "Pearson", 
#      "Spearman", and"Euclidean" (default)
#
#    optMov: option to perform moving average, true (default) or false. 
#    nMov: moving average size. Default is 5.
#
#    nEns: flow analogue ensemble size; for each ensemble, analogue is 
#      search every nEns days. Default is 5.
#
#    optCrop: option to crop certain region, true or false (default).
#    lonb: bounding box for longitude, ie. [-160,-80]
#    latb: bounding range for latitude, ie. [20, 60]
#
#    optEof: option to perform EOF and then use reconstructed data as 
#      input, true (default) or false.
#    thEof: cutoff explained variance (%) for EOF. Default is 90.

#    indYrRef: index of years to search for analogues, default is 1:N(year)

#    optNormalize: option to normalize data, include "none", "anomaly", 
#      "standardized anomaly", "smooth anomaly", and "smooth standardized 
#      anomaly" (default).
#    indClm: index of years for calculating climatology. Default is the 
#      same as the reference period (indYrRef)
#    nSmClm: window size for calculating smooth climatology. Default is 31.


#  Yizhou Zhuang (zhuangyz@atmos.ucla.edu)
#  References:
#    Zhuang et al. 2021 (PNAS)
#    Yiou et al. 2007 (GRL)
#

def flowanalogue_vpd(X, lon, lat, nFA, distfunc='euclidean'): 
    
    #check input size
    assert X.ndim == 3, 'Input X must be a 3-dimensional array' #originally 4-dimensional to split by day and year

    nlon = X.shape[2]
    nlat = X.shape[1]
    assert nlon == len(lon), 'Size mismatch between input X and lon'
    assert nlat == len(lat), 'Size mismatch between input X and lat'

#     nDoy = len(X.shape[2] / 365)  #number of days in year, must be <= 365
#     assert nDoy <= 365, 'Days of a year exceed 365 days'
    nDoy = 365

    nYr = X.shape[0] // 365   #number of years
    
    #default parameters
    nd = 30
    optMov = False
    nMov = 5
    nEns = 5
    optCrop = False
    lonb = [-160, -80] 
    latb = [20, 60]
    optEof = False
    thEof = 90 
    optNormalize = 'smooth standard anomaly'
    nSmClm = 31
    indYrRef = np.arange(nYr)
    indClm = indYrRef
    
#     #parameters with default value
#     nd = 30  #calendar window for analogue search is 2*nd+1 days centered on current day

#     optMov = True #option to perform moving average on input data X
#     nMov = 5 #moving average length

#     nEns = 5 #search analogue every nEns days, creating nEns ensembles

#     optCrop = False #option to use certain region for flow analogue search
#     lonb = bounds(lon) #bounding range for longitude, ie. [-160, -80]
#     latb = bounds(lat) #bounding range for latitude, ie. [20, 60]

#     optEof = True  #perform EOF on data X and then use reconstructed X as input
#     thEof = 90 #cutoff explained variance for EOF

#     #indYrRef = 1 : 32   #use 1979-2010 (input is 1979-2020)
#     indYrRef = 1 : nYr  #use all years

#     optNormalize = "smooth standardized anomaly"
#     indClm = indYrRef  #year index for calculating climatology; default is the same as reference period
#     nSmClm = 31 #window size for calculating smooth climatology

#     distFunc = "Euclidean", "Pearson", "Spearman"

    #Optional pre-processing
    #perform moving average for input X
    if optMov:
        X = X.rolling(time=nMov, min_periods=1, center=True).mean() #may still want to set endpoints to NAN
    #crop certain region of X for searching flow analogue
    if optCrop:
        X = X.where(X['longitude'] >= lonb[0]).dropna('longitude')
        X = X.where(X['longitude'] <= lonb[1]).dropna('longitude')
        X = X.where(X['latitude'] >= latb[0]).dropna('latitude')
        X = X.where(X['latitude'] <= latb[1]).dropna('latitude')
        nlon = X.shape[1]
        nlat = X.shape[2]
    #Option to normalize X data
    if optNormalize != "none":
        #from http://nicolasfauchereau.github.io/climatecode/posts/eof-analysis-with-scikit-learn/
#         init_scaler = preprocessing.StandardScaler()
#         X = X.stack(space=("latitude", "longitude"))
#         scaler = init_scaler.fit(X)
#         X = scaler.transform(X)
        X = (X - X.mean()) / X.std()
#         assert np.round(X.mean()) == 0.0
#         assert np.round(X.std()) == 1.0

#         X = xr.DataArray(X, dims=['time', 'longitude', 'latitude'])
    #perform EOF reconstruction
    ### DOUBLE CHECK THIS SECTION LATER
    if optEof: #ignore EOF parts for now
        yq, tmp = np.meshgrid(lat, lon)
        w = np.sqrt(np.cos(yq * np.pi/180)).flatten()  #grid weighting
        
        skpca = PCA()
        skpca.fit((X * w))
        ipc = np.where(skpca.explained_variance_ratio_.cumsum() >= thEof / 100)[0][0]
        print(ipc)
        P = skpca.transform((X * w))
        P = P[:,:ipc]
        E = skpca.components_
        E = E[:ipc,:].T
        print(X.shape, E.shape, P.shape, w.shape)
        X = np.matmul(E, P.T).T / w

        nCCA = P.shape[1] #use EOF number to repsent CCA number
        if nFA < nCCA:
            print('Warning! Analogue number is less than EOF number.')
    else:
        nCCA = nFA
    #Data preparation and output initialization
    #prepare input array X for flow analogue searching
    X0 = X.stack(space=("latitude", "longitude"))
    X = X.unstack("time").stack(space=("latitude", "longitude")).drop_vars(["year", "day1", "space", "latitude", "longitude"])
    
    nan_xr = xr.DataArray(np.full((nd, nlon*nlat), np.nan), dims=["day1", "space"], name="mean_vapor_pressure_deficit")
    X2 = xr.concat([nan_xr, X[0,:,:], X[1,:nd,:]], dim="day1") #add nd missing points at beginning
    ti = np.arange(nDoy*nYr) #create a temporary time index array
    ti2 = np.concatenate((np.empty(nd)*np.nan, ti[:nDoy+nd]))
    for j in range(1, nYr-1):
        X2 = xr.concat([X2, X[j-1,-nd:,:], X[j,:,:], X[j+1,:nd,:]], dim="day1") #add nd days on both ends to match window size
        ti2 = np.concatenate((ti2, ti[-nd+nDoy*j:nDoy*(j+1)+nd]))
    X2 = xr.concat([X2, X[-2,-nd:,:], X[-1,:,:], nan_xr], dim="day1")
    ti2 = np.concatenate((ti2, ti[-nd+nDoy*(nYr-1):], np.empty(nd)*np.nan))
    
    multi_index2 = pd.MultiIndex.from_product((range(44), range(425)), names=("year", "day"))
    X = X2.assign_coords(day1=multi_index2).unstack("day1").transpose("year", "day", "space")
    ti = np.reshape(ti2, (nYr,(nDoy+2*nd)))
#     X = np.concat(np.concat(np.empty(nlon*nlat,nd,1), X[:,nDoy-nd+1:,1:end-1], axis=2), 
#         X, np.concat(X[:,1:nd,2:], np.empty(nlon*nlat,nd,1), axis=2), axis=1)
#     if optEof:
#         nEof = size(P, 2)
#         P = reshape(P.T, nEof, nDoy, nYr); P0 = reshape(P, nEof, []);
#         P = np.concat(2, np.concat(3, nan(nEof,nd,1), P(:,nDoy-nd+1:end,1:end-1)), ...
#             P, np.concat(3, P[:,1:nd,2:end], nan(nEof,nd,1)))
    #create a temporary time index array
#     ti = np.reshape(np.arange(nDoy*nYr), (nYr, nDoy))
#     if nDoy == 365:
#         ti = [[np.empty(nd), ti[nDoy-nd+1:,:ti.shape[1]-1]], ti, [ti[:nd,1:], np.empty(nd)]]
#     else:
#         ti = [np.empty((nd,nYr)), ti, np.empty((nd,nYr))]
#     ti = np.concatenate((np.empty((nYr,nd)), ti, np.empty((nYr,nd))))
    #initialize distance array
    dist = np.empty((nFA, nYr, nDoy, nEns))
    i_d = dist 
    iDoy = dist
    iYr = dist
    b_cca = np.zeros((nFA, nYr, nDoy, nEns))

    #search flow analogue for each calendar day
    for i in range(nDoy):
        #temporary time index array for analogue ID assignment
        ti1 = ti[:,i:i+2*nd+1].flatten()
        #calculate distance between current day and all possible analogue days
        print(X[:,i:i+2*nd+1,:].stack(time=("year", "day")).transpose("time", "space"), '\n')
        print(X[:,i+nd+1,:])
        if distfunc == "pearson":
            X1 = X[:,i:i+2*nd+1,:].stack(time=("year", "day")).transpose("time", "space")
            dist0 = 1 - np.corrcoef(X1.T, X[:,i+nd+1,:].T, rowvar=False)[X1.shape[0]:,:X1.shape[0]].T
        if distfunc == "spearman":
            X1 = X[:,i:i+2*nd+1,:].stack(time=("year", "day")).transpose("time", "space")
            dist0 = 1 - stats.spearmanr(X1.T, X[:,i+nd+1,:].T)[0][X1.shape[0]:,:X1.shape[0]].T
        if distfunc == "euclidean":
            if optEof:
                dist0 = cdist(P[:,i:i+2*nd+1,:].stack(time=("year", "day")).transpose("time", "space"), P[:,i+nd+1,:])  
            else:
                X1 = X[:,i:i+2*nd+1,:].stack(time=("year", "day")).transpose("time", "space")
#                 dist0 = cdist(X1.to_numpy(), X[:,i+nd+1,:].to_numpy())
                dist0 = xr.apply_ufunc(pairwise_distances, X1.drop_vars(["time", "year", "day"]),
                                       X[:,i+nd+1,:].drop_vars(["year", "day"]),
                                       input_core_dims=[["time", "space"], ["space"]],
                                       output_core_dims=[["time", "time"]], vectorize=True, dask='allowed')
        #sort distance to current calendar day of all years to search analogue
        dist0_sorted = np.sort(dist0, axis=0)
        dist0_ind = np.argsort(dist0, axis=0)
        #get index of year
        dist0_indYr = np.floor(dist0_ind / (2*nd+1))
        #current year and years not in the reference period will not be considered
        FAindYr = (dist0_indYr != np.tile(np.arange(nYr),((2*nd+1)*nYr,1))) & (np.isin(dist0_indYr, indYrRef)) 
        #loop for all years
        for j in range(nYr):
            dist0 = dist0_sorted[FAindYr[:,j],j]
            id0 = ti1[dist0_ind[FAindYr[:,j],j]]
            if np.all(np.isnan(dist0)):
                dist[:,j,i,:] = np.nan
                i_d[:,j,i,:] = np.nan
                b_cca[:,j,i,:] = np.nan
            else:
                for k in range(nEns):
                    ind = np.nonzero((id0 % nEns == k).T.flatten())[0]
                    dist[:,j,i,k] = dist0[ind[:nFA]]
                    i_d[:,j,i,k] = id0[ind[:nFA]]
                    #CCA
                    nCCA1 = nCCA
                    if optEof:
                        while np.linalg.matrix_rank(X0[:,i_d[:nCCA1,j,i,k]]) < nCCA:
                            nCCA1 = nCCA1 + 1
                    B = LinearRegression(fit_intercept=False).fit(X0[i_d[:nCCA1,j,i,k].astype(int),:].T, X0[(j)*nDoy+i,:]).coef_
                    b_cca[:nCCA1,j,i,k] = B
        iDoy = i_d % nDoy
        iYr = np.floor(i_d/nDoy) 
    return dist, iDoy, iYr, i_d, b_cca

_, _, _, vpd_index_time_pearson, vpd_b_cca_pearson = flowanalogue_vpd(vpd_std_anom, gridmet_lon_array, gridmet_lat_array, 10, 'pearson')
