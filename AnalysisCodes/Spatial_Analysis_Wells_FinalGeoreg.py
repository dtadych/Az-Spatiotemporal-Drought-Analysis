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
datapath_web = 'https://data.cyverse.org/dav-anon/iplant/home/dtadych/AZ_Spatial_Analysis/Data/'
outputpath_web = 'https://data.cyverse.org/dav-anon/iplant/home/dtadych/AZ_Spatial_Analysis/Data/Output_files/'
shapepath_web = 'https://data.cyverse.org/dav-anon/iplant/home/dtadych/AZ_Spatial_Analysis/Data/Shapefiles/'

datapath_local = '../Data'
outputpath_local = '../Data/Output_files/'
shapepath_local = '../Data/Shapefiles/'

# %%  ==== Reading in the data we created above ====
#          Default is to read from the web, change as appropriate
# For regulation
filepath = outputpath_web+'/Waterlevels_Regulation.csv'
# filepath = outputpath_local+'Waterlevels_Regulation.csv'
cat_wl2_reg = pd.read_csv(filepath, index_col=0)

# For Access to SW
filepath = outputpath_web+'/Waterlevels_AccesstoSW.csv'
# filepath = outputpath_local+'Waterlevels_AccesstoSW.csv'
cat_wl2_SW = pd.read_csv(filepath, index_col=0)

# %% -- Linear regression --
# For Depth to Water by SW Access
ds = cat_wl2_SW
dt = "Access to Surface Water, Depth to Water"
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
data_type = "Regulation, Depth to Water"
min = 1975
mx = 2020
betterlabels = ['Regulated','Unregulated'] 

cf.linearregress(ds,data_type,min,mx,betterlabels)

# %% === WORKFLOW 2 ===

# Load in the master databases
filename_mdb_nd = 'Master_ADWR_database_noduplicates.shp'
filepath = os.path.join(outputpath_web, filename_mdb_nd)
print(filepath)

masterdb = gp.read_file(filepath)
pd.options.display.float_format = '{:.2f}'.format
print(masterdb.info())

# %%
filename_mdb_w = 'Master_ADWR_database_water.shp'
filepath = os.path.join(outputpath_web, filename_mdb_w)
print(filepath)

masterdb_water = gp.read_file(filepath)
pd.options.display.float_format = '{:.2f}'.format
print(masterdb_water.info())
# %%
# Reading in the shapefile
filename_georeg = 'georeg_reproject_fixed.shp'
filepath = os.path.join(shapepath_web, filename_georeg)
georeg = gp.read_file(filepath)
georeg.plot(cmap='viridis')

#%%
georeg['GEOREGI_NU'] = georeg['GEOREGI_NU'].astype('int64')
georeg.info()
#%%
# Read in the annual time series database
filename_ts = 'Wells55_GWSI_WLTS_DB_annual.csv'
filepath = os.path.join(outputpath_web, filename_ts)
print(filepath)
annual_db = pd.read_csv(filepath, header=1, index_col=0)
annual_db

# %%
annual_db.index = annual_db.index.astype('int64')

#%% 
only_special = masterdb[masterdb['WELL_TYPE_']=='OTHER']
only_special.info()
#%%
monitoring = masterdb[masterdb['WELL_TYPE_']=='MONITOR']
monitoring.info()

# %%
exempt = masterdb[masterdb['WELL_TYPE_']=='EXEMPT']
exempt.info()
#%%
nonexempt = masterdb[masterdb['WELL_TYPE_']=='NON-EXEMPT']
nonexempt.info()
# %% Overlay georegions onto the static database
# Going to use sjoin based off this website: https://geopandas.org/docs/user_guide/mergingdata.html
print("Non-cancelled: ", masterdb.crs, "Water Wells: ", masterdb_water.crs, "Georegions: ", georeg.crs)

# %%
georeg = georeg.to_crs(epsg=26912)
masterdb = masterdb.set_crs(epsg=26912)
masterdb_water = masterdb_water.set_crs(epsg=26912)
# %%
static_geo = gp.sjoin(masterdb, georeg, how="inner", op='intersects')
static_geo.head()
print(str(filename_mdb_nd) + " and " + str(filename_georeg) + " join complete.")

# %% Exporting or reading in the static geodatabase instead of rerunning
static_geo.to_csv(outputpath_local+'/Final_Static_geodatabase_allwells.csv')

# %% Rerunning this but for the water wells
static_geo2 = gp.sjoin(masterdb_water, georeg, how="inner", op='intersects')
static_geo2.head()
print(str(filename_mdb_nd) + " and " + str(filename_georeg) + " join complete.")

#%%
static_geo2.to_csv(outputpath_local+'/Final_Static_geodatabase_waterwells.csv')

# %%
filename = "Final_Static_geodatabase_allwells.csv"
filepath = os.path.join(outputpath_local, filename)
static_geo = pd.read_csv(filepath)
static_geo

# %% Create a dataframe of Final_Region and Well ID's
reg_list = static_geo[['Combo_ID', 'GEO_Region', 'GEOREGI_NU','Water_CAT', 'Loc','Regulation','WELL_DEPTH','WELL_TYPE_']]
reg_list

# %% Converting Combo_ID to int
reg_list['Combo_ID'] = reg_list['Combo_ID'].astype(int, errors = 'raise')

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
combo.to_csv(outputpath_local+'Final_WaterLevels_adjusted.csv')

# %% Reading in so we don't have to redo the combining, comment as appropriate
filepath = outputpath_web+'Final_WaterLevels_adjusted.csv'
# filepath = outputpath_local+'Final_WaterLevels_adjusted.csv'
combo = pd.read_csv(filepath, index_col=0)
combo.head()

# %% in order to filter deep/mid/shallow wells
shallow = 200
deep = 500

wd1 = combo[(combo["WELL_DEPTH"] > deep)]
wd2 = combo[(combo["WELL_DEPTH"] <= deep) & (combo["WELL_DEPTH"] >= shallow)]
wd3 = combo[(combo["WELL_DEPTH"] < shallow)]

# %% in order to make it where we can actually group these bitches
del combo['WELL_DEPTH']

# %%
combo_new = combo
combo_new['Old_Water_CAT'] = combo_new['Water_CAT']
combo_new

# %%
combo_new.loc[combo_new['Water_CAT']=='No_CAP'] = 'GW'
combo_new['Water_CAT'].unique()
# %% Now for aggregating by category for the timeseries
cat_wl_georeg = combo.groupby(['GEOREGI_NU']).mean()
cat_wl_reg = combo.groupby(['Regulation']).mean()
cat_wl_SW = combo.groupby(['Water_CAT']).mean()
cat_wl_SW
# %% 
cat_wl2_georeg = cat_wl_georeg.copy()
cat_wl2_reg = cat_wl_reg.copy()
cat_wl2_SW = cat_wl_SW.copy()

cat_wl2_georeg = cat_wl2_georeg.sort_values(by=['GEOREGI_NU'])
cat_wl2_SW = cat_wl2_SW.sort_values(by=['GEOREGI_NU'])

# Clean up the dataframe for graphing

i = cat_wl2_georeg
f = i.transpose()
f.reset_index(inplace=True)
f['index'] = pd.to_numeric(f['index'])
f['index'] = f['index'].astype(int)
f.set_index('index', inplace=True)
f.info()
cat_wl2_georeg = f
        
i = cat_wl2_reg
del i['GEOREGI_NU']
f = i.transpose()
f.reset_index(inplace=True)
f['index'] = pd.to_numeric(f['index'])
f['index'] = f['index'].astype(int)
f.set_index('index', inplace=True)
f.info()
cat_wl2_reg = f

i = cat_wl2_SW
del i['GEOREGI_NU']
f = i.transpose()
f.reset_index(inplace=True)
f['index'] = pd.to_numeric(f['index'])
f['index'] = f['index'].astype(int)
f.set_index('index', inplace=True)
f.info()
cat_wl2_SW = f

# %% Going to export all these as CSV's
cat_wl2_georeg.to_csv(outputpath_local+'Waterlevels_georegions.csv')
cat_wl2_reg.to_csv(outputpath_local+'Waterlevels_Regulation.csv')
cat_wl2_SW.to_csv(outputpath_local+'Waterlevels_AccesstoSW.csv')

# %%  ==== Reading in the data we created above ====
#          Default is to read from the web, change as appropriate
# For regulation
filepath = outputpath_local+'Waterlevels_Regulation.csv'
cat_wl2_reg = pd.read_csv(filepath, index_col=0)

# For Access to SW
filepath = outputpath_local+'Waterlevels_AccesstoSW.csv'
cat_wl2_SW = pd.read_csv(filepath, index_col=0)

# %% -- Linear regression --
# For Depth to Water by SW Access
ds = cat_wl2_SW
dt = "Access to Surface Water, Depth to Water"
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
data_type = "Regulation, Depth to Water"
min = 1975
mx = 2020
betterlabels = ['Regulated','Unregulated'] 

cf.linearregress(ds,data_type,min,mx,betterlabels)
