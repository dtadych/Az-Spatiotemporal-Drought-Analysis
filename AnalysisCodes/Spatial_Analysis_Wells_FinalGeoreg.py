# ~~    Spatial Analysis Code   ~~
# Written by Danielle Tadych

# The purpose of this script is to create a code to spatially analyze all the wells in 
# the combined database based on management. 

# WORKFLOW 1
# 1. Read in the data from Cyverse
# 2. Run Linear regression from custom functions

# WORKFLOW 2
# 1. Read in the master ADWR database static database, water level database, and 
#       georegions shapefile created in QGIS
# 2. Overlay region shapefile on static well database shapefile
# 3. Export a dataframe (registry list) of combined ID's with the columns we want 
#       (regulation, etc.)
# 4. Join the registry list with the timeseries database so every well has water 
#       levels and is tagged with a category we want
# 5. Create pivot tables averaging water levels based on categories (e.g. regulation, 
#       access to SW, or georegion (finer scale))
# 6. Export pivot tables into .csv's for easy analyzing later
#       * Note: after reading in packages, skip to line 233 to avoid redoing steps 1-5
#         or to read in the web .csv's
# 8. Run Linear regression from custom functions

# %%
from optparse import Values
import os
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from matplotlib.colors import ListedColormap
import datetime
from matplotlib.transforms import Bbox
import numpy as np
import pandas as pd
from shapely.geometry import box
import geopandas as gp
#import earthpy as et
import scipy.stats as sp
from scipy.stats import kendalltau, pearsonr, spearmanr
import pymannkendall as mk
import Custom_functions as cf

# Data paths
datapath_local = '../Data'
inputpath_local = '../Data/Input_files'
outputpath_local = '../Data/Output_files/'
shapepath_local = '../Data/Shapefiles/'

# %%  ==== WORFKLOW 1 ====
#          Default is to read from the web, change as appropriate
# For regulation
filepath = outputpath_local+'/Waterlevels_Regulation_updated.csv'
cat_wl2_reg = pd.read_csv(filepath, index_col=0)

# For Access to SW
filepath = outputpath_local+'/Waterlevels_AccesstoSW_updated.csv'
cat_wl2_SW = pd.read_csv(filepath, index_col=0)

# %% -- Linear regression --
# For Depth to Water by SW Access
ds = cat_wl2_SW
dt = "Access to Surface Water, Depth to Water (Updated Data)"
min = 1975
mx = 2020
betterlabels = ['Recieves CAP (Regulated)'
                ,'GW Dominated (Regulated)'
                ,'Surface Water Dominated'
                ,'GW Dominated'
                ,'Mixed Source']
cf.linearregress(ds,dt,min,mx,betterlabels)
#%%
ds = cat_wl2_reg
data_type = "Regulation, Depth to Water (Updated Data)"
min = 1975
mx = 2020
betterlabels = ['Regulated','Unregulated'] 

cf.linearregress(ds,data_type,min,mx,betterlabels)

# %% === WORKFLOW 2 ===

# Load in the master databases
filename_mdb_nd = 'Master_ADWR_database_noduplicates.shp'
filepath = os.path.join(inputpath_local, filename_mdb_nd)
print(filepath)

masterdb = gp.read_file(filepath)
pd.options.display.float_format = '{:.2f}'.format
print(masterdb.info())

# %%
# Reading in the shapefile
filename_georeg = 'georeg_reproject_fixed.shp'
filepath = os.path.join(shapepath_local, filename_georeg)
georeg = gp.read_file(filepath)
georeg.plot(cmap='viridis')

#%%
georeg['GEOREGI_NU'] = georeg['GEOREGI_NU'].astype('int64')
georeg.info()

# %% Overlay georegions onto the static database
# Going to use sjoin based off this website: https://geopandas.org/docs/user_guide/mergingdata.html
print("Non-cancelled: ", masterdb.crs, "Georegions: ", georeg.crs)

# %%
georeg = georeg.to_crs(epsg=26912)
masterdb2 = masterdb.set_crs(epsg=26912)
# %%
static_geo = gp.sjoin(masterdb2, georeg, how="inner", op='intersects')
static_geo.head()
print(str(filename_mdb_nd) + " and " + str(filename_georeg) + " join complete.")

# %% Exporting or reading in the static geodatabase instead of rerunning
# static_geo.to_csv(outputpath_local+'/Final_Static_geodatabase_allwells.csv')

# %%
# filename = "Final_Static_geodatabase_allwells.csv"
# filepath = os.path.join(outputpath_local, filename)
# static_geo = pd.read_csv(filepath)
# static_geo

# %% Create a dataframe of Final_Region and Well ID's
reg_list = static_geo[['Combo_ID', 'GEO_Region', 'GEOREGI_NU','Water_CAT', 'Loc','Regulation','WELL_DEPTH','WELL_TYPE_']]
reg_list

# %% Converting Combo_ID to int
reg_list['Combo_ID'] = reg_list['Combo_ID'].astype(np.int64, errors = 'raise')

#%%
# Read in the annual time series database
filename_ts = 'Wells55_GWSI_WLTS_DB_annual_updated.csv'
# filename_ts = 'Wells55_GWSI_MAX_WLTS_DB_annual_updated.csv'
# filename_ts = 'Wells55_GWSI_MIN_WLTS_DB_annual_updated.csv'
filepath = os.path.join(outputpath_local, filename_ts)
print(filepath)
annual_db = pd.read_csv(filepath, header=1, index_col=0)
annual_db

# %%
annual_db = annual_db[1:168102]
annual_db
# %%
annual_db.index = annual_db.index.astype('int64')
# %%
annual_db2 = annual_db.reset_index(inplace=True)
annual_db2 = annual_db.rename(columns = {'year':'Combo_ID'})
annual_db2.head()

# %% Add list to the annual database
combo = annual_db2.merge(reg_list, how="outer")
combo.info()

# %% set index to Combo_ID
combo.set_index('Combo_ID', inplace=True)

# %% Sort the values
combo = combo.sort_values(by=['GEOREGI_NU'])
combo

# %% Exporting the combo table
# combo.to_csv(outputpath_local+'Final_WaterLevels_adjusted_updated.csv')

# %% Reading in so we don't have to redo the combining, comment as appropriate
# filepath = outputpath_local+'Final_WaterLevels_adjusted.csv'
# filepath = outputpath_local+'Final_WaterLevels_adjusted.csv'
# combo = pd.read_csv(filepath, index_col=0)
# combo.head()

# %%
combo_new = combo.copy()
combo_new['Old_Water_CAT'] = combo_new['Water_CAT']
combo_new

# %%
combo_new.loc[combo_new['Water_CAT']=='No_CAP'] = 'GW'
combo_new['Water_CAT'].unique()
# %% Now for aggregating by category for the timeseries (Spatial Average)
cat_wl_reg = combo.groupby(['Regulation']).mean()
cat_wl_SW = combo.groupby(['Water_CAT']).mean()
cat_wl_SW
# %% 
cat_wl2_reg = cat_wl_reg.copy()
cat_wl2_SW = cat_wl_SW.copy()

cat_wl2_SW = cat_wl2_SW.sort_values(by=['GEOREGI_NU'])

# Clean up the dataframe for graphing
        
i = cat_wl2_reg
del i['GEOREGI_NU']
del i['WELL_DEPTH']
f = i.transpose()
f.reset_index(inplace=True)
f['index'] = pd.to_numeric(f['index'])
f['index'] = f['index'].astype(int)
f.set_index('index', inplace=True)
f.info()
cat_wl2_reg = f

i = cat_wl2_SW
del i['GEOREGI_NU']
del i['WELL_DEPTH']
f = i.transpose()
f.reset_index(inplace=True)
f['index'] = pd.to_numeric(f['index'])
f['index'] = f['index'].astype(int)
f.set_index('index', inplace=True)
f.info()
cat_wl2_SW = f

# %% Going to export all these as CSV's
# cat_wl2_reg.to_csv(outputpath_local+'Waterlevels_Regulation_updated.csv')
# cat_wl2_SW.to_csv(outputpath_local+'Waterlevels_AccesstoSW_updated.csv')

# cat_wl2_reg.to_csv(outputpath_local+'Waterlevels_Regulation_updated_MAX.csv')
# cat_wl2_SW.to_csv(outputpath_local+'Waterlevels_AccesstoSW_updated_MAX.csv')

# cat_wl2_reg.to_csv(outputpath_local+'Waterlevels_Regulation_updated_MIN.csv')
# cat_wl2_SW.to_csv(outputpath_local+'Waterlevels_AccesstoSW_updated_MIN.csv')
# %%  ==== Reading in the data we created above ====
#          Default is to read from the web, change as appropriate
# For regulation
filepath = outputpath_local+'Waterlevels_Regulation_updated.csv'
cat_wl2_reg = pd.read_csv(filepath, index_col=0)

# For Access to SW
filepath = outputpath_local+'Waterlevels_AccesstoSW_updated.csv'
cat_wl2_SW = pd.read_csv(filepath, index_col=0)

# %% This might be necessary
del cat_wl2_reg['Res'], cat_wl2_SW['Res']

# %% -- Linear regression --
# For Depth to Water by SW Access
ds = cat_wl2_SW
# dt = "Access to Surface Water, Depth to Water"
# dt = "Access to Surface Water, MAX Depth to Water"
dt = "Access to Surface Water, MIN Depth to Water"
min = 1975
mx = 2022
betterlabels = ['Recieves CAP (Regulated)'
                ,'GW Dominated (Regulated)'
                ,'Surface Water Dominated'
                ,'GW Dominated'
                ,'Mixed Source']
cf.linearregress(ds,dt,min,mx,betterlabels)
#%%
ds = cat_wl2_reg
data_type = "Access to Surface Water, Depth to Water"
# dt = "Access to Surface Water, MAX Depth to Water"
# dt = "Access to Surface Water, MIN Depth to Water"
min = 1975
mx = 2022
betterlabels = ['Regulated','Unregulated'] 

cf.linearregress(ds,data_type,min,mx,betterlabels)

# %%
