# === GRACE Spatial Analysis Script ===
# written by Danielle Tadych
# The purpose of this script is to analyze GRACE Data for Arizona by points and shapes
#  - Importing packages: Line 12
#  - Reading in files: Line 36
#  - EASYMORE Remapping using a shapefile to computer zonal statistics: Line 56
#     *Note: in order for this package to work
#               > the projections need to be in espg4326, exported, re-read in
#               > the time variable for nc needs to be datetime
#               > Value error -> run the "fix geometries" tool in qgis
#  - Calculating the average based off a mask (not weighted): Line 212
# %%
from itertools import count
import os
from tkinter import Label
from typing import Mapping
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from matplotlib.colors import ListedColormap
import datetime as dt
import seaborn as sns
import numpy as np
import pandas as pd
from shapely.geometry import geo
import geopandas as gp
import xarray as xr
import rioxarray as rxr
import netCDF4
import rasterio
from scipy.stats import kendalltau, pearsonr, spearmanr
import easymore
import glob
import scipy.stats as sp

# %% Read in the file
filename = 'CSR_GRACE_GRACE-FO_RL06_Mascons_all-corrections_v02.nc'
datapath = '../Data/Input_files/GRACE/'
outputpath = '../Data/Output_files/'
shapepath = '../Data/Shapefiles/'

grace_dataset = xr.open_dataset(datapath+'/'+filename)
grace_dataset

# %% Read in the mask shapefile
filename = "AZ_counties.shp"
filepath = os.path.join(shapepath, filename)
counties = gp.read_file(filepath)

filename_georeg = 'georeg_reproject_fixed.shp'
filepath = os.path.join(shapepath, filename_georeg)
georeg = gp.read_file(filepath)
# %% Look at that sweet sweet data
metadata = grace_dataset.attrs
metadata

# %% View first 5 values
grace_dataset["lwe_thickness"]["lat"].values[:5]
print("The min and max latitude values in the data is:", 
      grace_dataset["lwe_thickness"]["lat"].values.min(), 
      grace_dataset["lwe_thickness"]["lat"].values.max())
print("The min and max longitude values in the data is:", 
      grace_dataset["lwe_thickness"]["lon"].values.min(), 
      grace_dataset["lwe_thickness"]["lon"].values.max())

print("The earliest date in the data is:", 
    grace_dataset["lwe_thickness"]["time"].values.min())
print("The latest date in the data is:", 
    grace_dataset["lwe_thickness"]["time"].values.max())

# %%
print("Number of Datapoints")
grace_dataset["lwe_thickness"]['time'].values.shape   

# %% Slicing data to get variables
lat = grace_dataset.variables['lat'][:]
lon = grace_dataset.variables['lon'][:]
time = grace_dataset.variables['time'][:]
lwe = grace_dataset['lwe_thickness']
lwe

# %% Now I need to assign a coordinate system to lwe
lwe.coords['lon'] = (lwe.coords['lon'] + 180) % 360 - 180
lwe2 = lwe.sortby(lwe.lon)
lwe2 = lwe2.rio.set_spatial_dims('lon', 'lat')
lwe2 = lwe2.rio.set_crs("epsg:4269")
lwe2.rio.crs

# %% Convert time to datetime format
time = grace_dataset.variables['time'] # do not cast to numpy array yet 
time_convert = netCDF4.num2date(time[:], "days since 2002-01-01T00:00:00Z", calendar='standard')
lwe2['time'] = time_convert
datetimeindex = lwe2.indexes['time'].to_datetimeindex()
lwe2['time'] = datetimeindex


# %% ---- Remapping using EASYMORE Package ----
# Fixing GRACE to be datetime suitable for easymore.
grace2 = grace_dataset
grace2
# %%
grace2['time'] = datetimeindex
grace2
# %% Write the full nc file
fn = "testGRACE_time.nc"
grace2.to_netcdf(fn)
# %% Now remapping following this tutorial
# https://github.com/ShervanGharari/EASYMORE/blob/main/examples/Chapter1_E1.ipynb
# # loading EASYMORE
from easymore.easymore import easymore

# initializing EASYMORE object
esmr = easymore()

# specifying EASYMORE objects
# name of the case; the temporary, remapping and remapped file names include case name
esmr.case_name                = 'easymore_GRACE_georeg'              
# temporary path that the EASYMORE generated GIS files and remapped file will be saved
esmr.temp_dir                 = '../temporary/'

# name of target shapefile that the source netcdf files should be remapped to
esmr.target_shp = '../Data/Shapefiles/georeg_reproject_fixed.shp'

# name of netCDF file(s); multiple files can be specified with *
esmr.source_nc                = 'testGRACE_time.nc'

# name of variables from source netCDF file(s) to be remapped
esmr.var_names                = ["lwe_thickness"]
# rename the variables from source netCDF file(s) in the remapped files;
# it will be the same as source if not provided
esmr.var_names_remapped       = ["lwe_thickness"]
# name of variable longitude in source netCDF files
esmr.var_lon                  = 'lon'
# name of variable latitude in source netCDF files
esmr.var_lat                  = 'lat'
# name of variable time in source netCDF file; should be always time
esmr.var_time                 = 'time'
# location where the remapped netCDF file will be saved
esmr.output_dir               = outputpath
# format of the variables to be saved in remapped files,
# if one format provided it will be expanded to other variables
esmr.format_list              = ['f4']
# fill values of the variables to be saved in remapped files,
# if one value provided it will be expanded to other variables
esmr.fill_value_list          = ['-9999.00']
# if required that the remapped values to be saved as csv as well
esmr.save_csv                 = True
esmr.complevel                 =  9
# if uncommented EASYMORE will use this and skip GIS tasks
#esmr.remap_csv                = '../temporary/ERA5_Medicine_Hat_remapping.csv'
# %% This code can take a while to run, for our machines it was ~18 minutes
# execute EASYMORE
esmr.nc_remapper()

# %% Setting up new remapped variables
# visualize the remapped netCDF for the first file, first time step
# target nc file
nc_names = sorted(glob.glob (esmr.output_dir + esmr.case_name + '*.nc'))
ds       = xr.open_dataset(nc_names[0]) # the first netcdf file
values   = ds.lwe_thickness [0,:] # the first time frame of the first 
IDs      = ds.ID [:] # get the ID
# create a data frame for the model simulation
df = pd.DataFrame()
df ['value'] = values
df ['ID_t']    = IDs  # .astype(int)
df = df.sort_values (by = 'ID_t')
# load the shape file target that is generated by EASYMORE (with consistent IDs)
shp_target = gp.read_file(esmr.temp_dir+ esmr.case_name + '_target_shapefile.shp') # load the target shapefile
shp_target ['ID_t'] = shp_target ['ID_t'].astype(float)
shp_target = shp_target.sort_values(by='ID_t')# sort on values
shp_target = pd.merge_asof(shp_target, df, on='ID_t', direction='nearest')
shp_target = shp_target.set_geometry('geometry') #bring back the geometry filed; pd to gpd
# plotting to make sure it looks okay
f, axes = plt.subplots(1,1,figsize=(15,15))
shp_target.plot(column= 'value', edgecolor='k',linewidth = 1, ax = axes , legend=True)


# %% Now, read in the remapped csv
filename = 'easymore_GRACE_georeg_remapped_lwe_thickness__2002-04-18-00-00-00.csv'
filepath = os.path.join(outputpath, filename)
grace_remapped = pd.read_csv(filepath)
grace_remapped.head()
# %% Assign category Headers
ID_key = shp_target[['GEO_Region', 'ID_t']]
ID_key['ID'] = 'ID_' + ID_key['ID_t'].astype(str)
grace_remapped = grace_remapped.set_index('time')
del grace_remapped['Unnamed: 0'] # tbh not really sure why this column is here but gotta delete it
grace_remapped
georeg_list = ID_key['GEO_Region'].values.tolist()
grace_remapped.columns = georeg_list

# %% Fixing the time element
grace_remapped.index = pd.to_datetime(grace_remapped.index)
grace_remapped.plot()

# %% Averaging based on year
grace_yearly = grace_remapped
grace_yearly['year'] = pd.DatetimeIndex(grace_yearly.index).year
grace_yearly = grace_yearly.reset_index()
grace_yearlyavg = pd.pivot_table(grace_yearly, index=["year"], dropna=False, aggfunc=np.mean)
grace_yearlyavg = grace_yearlyavg.reset_index()
grace_yearlyavg['year'] = pd.to_numeric(grace_yearlyavg['year'])
grace_yearlyavg['year'] = grace_yearlyavg['year'].astype(int)
grace_yearlyavg = grace_yearlyavg.set_index('year')
print("Annual values calculated.")
grace_yearlyavg.plot() #to see if it completed

# %% Write a .csv for now for graphing
grace_remapped.to_csv(outputpath+'grace_remapped.csv')
grace_yearlyavg.to_csv(outputpath+'grace_remapped_yearly.csv')

# %%
# ---- Creating Averages Based off Shape File Mask ----
# Check the cooridnate systems
mask = counties
print("mask crs:", counties.crs)
print("data crs:", lwe2.rio.crs)

# %% Clipping based off the mask (not weighted)
clipped = lwe2.rio.clip(mask.geometry, mask.crs)
print("Clipping finished.")
# %% Checking to see if it clipped correctly
fig, ax = plt.subplots(figsize=(6, 6))
clipped[0,:,:].plot()
mask.boundary.plot(ax=ax)
ax.set_title("Shapefile Clip Extent",
             fontsize=16)
plt.show()

# %%
clipped['time'] = datetimeindex
clipped_mean = clipped.mean(("lon","lat"))
clipped_mean
cm_df = pd.DataFrame(clipped_mean)
cm_df = cm_df.reset_index()
cm_df['index'] = datetimeindex
cm_df.set_index('index', inplace=True)
cm_df

# %%
# Extract the year from the date column and create a new column year
cm_df['year'] = pd.DatetimeIndex(cm_df.index).year
cm_df_year = pd.pivot_table(cm_df, index=["year"], values=[0], dropna=False, aggfunc=np.mean)
cm_df_year

# %% Write the csv
cm_df_year.to_csv(outputpath+'grace_stateavg_yearly.csv')


# %% Creating Arizona Specific Anomalies (Changes from Az Average)
remap_anom = grace_yearlyavg.copy()
for i in remap_anom.columns:
    remap_anom[i] = remap_anom[i]-cm_df_year[0]
remap_anom
remap_anom.plot()

# %% Write the csv
remap_anom.to_csv(outputpath+'grace_remappedanomalies_yearly.csv')

# %% Linear Regression
# For Depth to Water by SW Access
ds = remap_anom
data_type = "LWE (cm)"
min_yr = 2002
mx_yr = 2020
betterlabels = ['CAP','Regulated Groundwater','Surface Water','Unregulated Groundwater','Mixed GW/SW'] 
Name = str(min_yr) + " to " + str(mx_yr) + " Linear Regression for " + data_type
print(Name)

f = ds[(ds.index >= min_yr) & (ds.index <= mx_yr)]
columns = ds.columns
column_list = ds.columns.tolist()
columns = ds.columns
column_list = ds.columns.tolist()
# ------------------------

stats = pd.DataFrame()
for i in column_list:
        df = f[i]
        # df = f[i].pct_change()
        #print(df)
        y=np.array(df.values, dtype=float)
        x=np.array(pd.to_datetime(df).index.values, dtype=float)
        slope, intercept, r_value, p_value, std_err =sp.linregress(x,y)
        stats = stats.append({'slope': slope, 
                              'int':intercept, 
                              'rsq':r_value*r_value, 
                              'p_val':p_value, 
                              'std_err':std_err, 
                              'mean': np.mean(y),
                              'var': np.var(y),
                              'sum': np.sum(y)
                              },
                              ignore_index=True)


stats.index = column_list
stats1 = stats.transpose()
stats1.to_csv(outputpath+'Stats/'+Name+'.csv')
# %%
print('Code completed. Now move to "Graphs.py"')
