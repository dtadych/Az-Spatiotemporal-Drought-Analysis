# The purpose of this script is to make multiple timeseries databases using data from the GWSI and Wells55 databases
# Written by Danielle Tadych & Matt Ford
# Goals:
# - Create columns in each respective database specifying its origin
# - Find column they have in common
# - Merge based on that column
# - Make GWSI wells the overriding database and fill in the gaps with Wells55

# - Make 3 Timeseries databases: WL, Water Elevation, and Pumping
#     * Note: Pumping might not be available in Wells55
#             Would need to potentially multiply pumping amounts with how long it has been installed
#             Check with Laura
# - Make the columns well ID's and the rows dates
# - for Wells55 water level and elevation, need to make rows the install date and columns REGISTRY_I
# - Merge based on well ID's

# Status as of 7/5/2021:
# - Made water level (depth to water below land surface in feet) timeseries databases
# %%
import os
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import datetime
from pandas.tseries.offsets import BYearBegin
import seaborn as sns
import geopandas as gp

datapath = '../Data/Input_files/'
outputpath = '../Data/Output_files/'

<<<<<<< HEAD
=======
# %% ---- First Creating the GWSI Water level ----
# Skip to line 72 if downloading the data from cyverse
# Read in the the water level file
GWSI_folder = '../Data/Input_files/GWSI/Data_Tables' #GWSI folder name
file_name = 'GWSI_WW_LEVELS.xlsx'
filepath=os.path.join(GWSI_folder, file_name)
print(filepath)

wl_data = pd.read_excel(filepath, parse_dates=['WLWA_MEASUREMENT_DATE'])
print(wl_data.info())

# Rename the columns in wl file to something shorter
wl_data=wl_data.rename(columns={"WLWA_MEASUREMENT_DATE": "date",
                   "WLWA_SITE_WELL_SITE_ID": "wellid",
                   "WLWA_DEPTH_TO_WATER": "depth"}, errors="raise")
# %%
# Read in the file containing basin codes
file_name = 'GWSI_SITES.xlsx'
filepath=os.path.join(GWSI_folder, file_name)
print(filepath)

basin_data = pd.read_excel(filepath)
print(basin_data.info())

# Rename the columns in basin file to something shorter
basin_data=basin_data.rename(columns={"SITE_WELL_SITE_ID": "wellid",
                   "SITE_ADWBAS_CODE_ENTRY": "basinid"}, errors="raise")

# print number of columns
basin_data['basinid'].nunique()            
# %%
# Merge basin codes and water levels by wellid into new datatable called wl_data2
wl_data2 = wl_data.merge(basin_data, left_on='wellid', right_on='wellid')
print(wl_data2.info())

#%%
# Output wl_data2 to a csv in the specified directory
wl_data2.to_csv(outputpath+'wl_data3.csv')

>>>>>>> 00e8cedb276036a5996ed62c2a5cdef2568d0460
# %% 
# ----- Import the Data and Shapefiles with Geometries -----
# Read in Wells 55 Data
# This is a file with water levels from ADWR which has been joined with another ADWR file with variables
filename = 'Well_Registry_05032023.csv'
filepath = os.path.join(datapath, filename)
print(filepath)

wells55 = pd.read_csv(filepath)
pd.options.display.float_format = '{:.2f}'.format
print(wells55.info())

# Read in GWSI collated water level data
filename = 'wl_data3.csv'
<<<<<<< HEAD
=======
# filename = 'wl_data2.csv'
>>>>>>> 00e8cedb276036a5996ed62c2a5cdef2568d0460
filepath = os.path.join(datapath, filename)
print(filepath)

wl_data2 = pd.read_csv(filepath)
pd.options.display.float_format = '{:.2f}'.format
print(wl_data2.info())

#%% ---- Making dataframes with Date as the index ---
# Make dataframes with Columns for GWSI and Wells55, respectively
# Following this method: https://stackoverflow.com/questions/32215024/merging-time-series-data-by-timestamp-using-numpy-pandas
# Confirmed through the variables list that "depth" in wldata2 and "WATER_LEVE" are both depth to water below land surface in feet

gwsi_wl = wl_data2[["date","wellid","SITE_WELL_REG_ID","depth"]].copy()
gwsi_wl.info()

# gwsi_wl = wl_data2[["date","SITE_WELL_REG_ID","depth"]].copy()
# gwsi_wl.info()

wells55_wl = wells55[["INSTALLED", "REGISTRY_ID", "WATER_LEVEL"]].copy()
wells55_wl.info()

# Need to add an original database column
wells55_wl["Original_DB"] = 'Wells55'
gwsi_wl["Original_DB"] = 'GWSI'
wells55_wl.head()

# Renaming columns
wells55_wl.rename(columns = {'INSTALLED':'date','WATER_LEVEL':'depth'}, inplace=True)
<<<<<<< HEAD
# Need to create a combo ID column
gwsi_wl['Combo_ID'] = gwsi_wl.SITE_WELL_REG_ID.combine_first(gwsi_wl.wellid)
gwsi_wl.rename(columns={'Combo_ID':'REGISTRY_ID'}, inplace=True)

# Need to make sure the columns are the same Dtype
=======
# %% Need to create a combo ID column
gwsi_wl['Combo_ID'] = gwsi_wl.SITE_WELL_REG_ID.combine_first(gwsi_wl.wellid)
gwsi_wl.info()

# %%
gwsi_wl.rename(columns={'Combo_ID':'REGISTRY_ID'}, inplace=True)
# gwsi_wl.rename(columns={'SITE_WELL_REG_ID':'REGISTRY_ID'}, inplace=True)
gwsi_wl.info()

# %% Need to make sure the columns are the same Dtype
>>>>>>> 00e8cedb276036a5996ed62c2a5cdef2568d0460
gwsi_wl.REGISTRY_ID = gwsi_wl.REGISTRY_ID.astype('int64')
gwsi_wl.info()

# Setting the date column to datetime
wells55_wl['date'] = pd.to_datetime(wells55_wl.date)
wells55_wl['date'] = wells55_wl['date'].dt.tz_localize(None)
wells55_wl.info()

# 
gwsi_wl.date = pd.to_datetime(gwsi_wl.date)
gwsi_wl.info()

# Merging 
combo = wells55_wl.merge(gwsi_wl, suffixes=['_wells55','_gwsi'], how="outer" 
                                          ,on=["REGISTRY_ID", 'date', 'Original_DB', 'depth']
                                          )
combo.info()

<<<<<<< HEAD
=======
# %% This block of code takes a really long time to run so I would suggest skipping
# WL_TS_DB = pd.pivot_table(combo, index=["REGISTRY_ID"], columns="date", values="depth")
# WL_TS_DB.head()
# # Export data into a csv
# WL_TS_DB.to_csv(outputpath + 'Wells55_GWSI_WLTS_DB.csv')

>>>>>>> 00e8cedb276036a5996ed62c2a5cdef2568d0460
# %% --- Summarizing the data by date ---
# Extract the year from the date column and create a new column year
combo['year'] = pd.DatetimeIndex(combo.date).year
combo.head()

# %% Creating pivot tables
WL_TS_DB_year = pd.pivot_table(combo, index=["REGISTRY_ID"], columns=["year"], values=["depth"], dropna=False, aggfunc=np.mean)
LEN_TS_DB_year = pd.pivot_table(combo, index=["REGISTRY_ID"], columns=["year"], values=["depth"], dropna=False, aggfunc=len)
max_TS_DB_year = pd.pivot_table(combo, index=["REGISTRY_ID"], columns=["year"], values=["depth"], dropna=False, aggfunc=np.max)
<<<<<<< HEAD
min_TS_DB_year = pd.pivot_table(combo, index=["REGISTRY_ID"], columns=["year"], values=["depth"], dropna=False, aggfunc=np.min)
=======
# %% Testing 1980 versus 2020 to see if there's a difference
print(WL_TS_DB_year.iloc[:,115])
# %%
print(WL_TS_DB_year.iloc[:,155])

# %%
WL_TS_DB_1980 = pd.DataFrame(WL_TS_DB_year.iloc[:,115])
WL_TS_DB_2020 = pd.DataFrame(WL_TS_DB_year.iloc[:,155])
WL_TS_DB_1980.index.name = None
WL_TS_DB_2020.index.name = None
WL_TS_DB_1980.columns = ["depth"]
WL_TS_DB_2020.columns = ["depth"]
>>>>>>> 00e8cedb276036a5996ed62c2a5cdef2568d0460

#%%
WL_TS_DB_year.index.name = None
WL_TS_DB_year.head()
# %%
test = WL_TS_DB_year.mean()
test.plot()

# %% 
# stats = WL_TS_DB_year.describe()
stats = stats.transpose()
stats2 = stats[['mean','25%','50%','75%']]
stats2[119:157].plot()

# %%
narrowedstats = stats[110:158]
narrowedstats
# %%
narrowedstats.to_csv(outputpath+"state_average_WL_updated.csv")
<<<<<<< HEAD

# %%
# Export both yearly summary data and monthly into csv
WL_TS_DB_year.to_csv(outputpath + 'Wells55_GWSI_WLTS_DB_annual_updated.csv')
=======
# %%
max = stats['max']
max[119:157].plot()
# %% Exporting data
# WL_TS_DB_1980.to_csv(outputpath + 'comboDB_WL_1980.csv')
# WL_TS_DB_2020.to_csv(outputpath + 'comboDB_WL_2020.csv')
>>>>>>> 00e8cedb276036a5996ed62c2a5cdef2568d0460
# %%
LEN_TS_DB_year.to_csv(outputpath + 'Wells55_GWSI_LEN_WLTS_DB_annual_updated.csv')

# %%
max_TS_DB_year.to_csv(outputpath + 'Wells55_GWSI_MAX_WLTS_DB_annual_updated.csv')

# %%
min_TS_DB_year.to_csv(outputpath + 'Wells55_GWSI_MIN_WLTS_DB_annual_updated.csv')

# %% Creating totals for mapping
min_yr = 1975.0
mx_yr = 2023.0
# ds = WL_TS_DB_year.transpose()
ds = LEN_TS_DB_year.loc[:, ('depth', min_yr):('depth', mx_yr)]
total_WL = ds.sum(axis=1)
total_WL = pd.DataFrame(total_WL)
total_WL = total_WL.reset_index()
total_WL = total_WL.rename(columns={"REGISTRY_ID": "Combo_ID",
                   0: "LEN"}, errors="raise")
total_WL = total_WL.set_index('Combo_ID')
# total_WL = total_WL['LEN'].astype(int, errors = 'raise')
total_WL
total_WL.to_csv(outputpath+'numberWL_perwell_'+str(min_yr)+'-'+str(mx_yr)+'_comboID.csv')
# total_WL.to_csv(outputpath+'numberWL_perwell_'+str(min_yr)+'-'+str(mx_yr)+'_updated.csv')

# %%
# Summarize by Monthly now
combo_monthly2 = combo.copy()
combo_monthly2.reset_index(inplace=True)
combo_monthly2.date = combo.date.dt.strftime('%Y-%m')
combo_monthly2.info()
# %%
# Average
WL_TS_DB_month = pd.pivot_table(combo_monthly2, index=["REGISTRY_ID"], columns=["date"], values=["depth"], dropna=False, aggfunc='mean')
WL_TS_DB_month.index.name = None
<<<<<<< HEAD

# Number of measurements
WL_TS_DB_month_len = pd.pivot_table(combo_monthly2, index=["REGISTRY_ID"], columns=["date"], values=["depth"], dropna=False, aggfunc=len)
WL_TS_DB_month_len.index.name = None

# Max
WL_TS_DB_month_max = pd.pivot_table(combo_monthly2, index=["REGISTRY_ID"], columns=["date"], values=["depth"], dropna=False, aggfunc=np.max)
WL_TS_DB_month_max.index.name = None

# Min
WL_TS_DB_month_min = pd.pivot_table(combo_monthly2, index=["REGISTRY_ID"], columns=["date"], values=["depth"], dropna=False, aggfunc=np.min)
WL_TS_DB_month_min.index.name = None
# %% Writing CSV's
=======
WL_TS_DB_month.head()
# %%
# Export both yearly summary data and monthly into csv
# WL_TS_DB_year.to_csv(outputpath + 'Wells55_GWSI_WLTS_DB_annual.csv')
# WL_TS_DB_year.to_csv(outputpath + 'Wells55_GWSI_WLTS_DB_annual_comboID.csv')
WL_TS_DB_year.to_csv(outputpath + 'Wells55_GWSI_WLTS_DB_annual_updated.csv')
# %%
# LEN_TS_DB_year.to_csv(outputpath + 'Wells55_GWSI_LEN_WLTS_DB_annual.csv')
LEN_TS_DB_year.to_csv(outputpath + 'Wells55_GWSI_LEN_WLTS_DB_annual_comboID.csv')
# LEN_TS_DB_year.to_csv(outputpath + 'Wells55_GWSI_LEN_WLTS_DB_annual_updated.csv')

# %%
# max_TS_DB_year.to_csv(outputpath + 'Wells55_GWSI_MAX_WLTS_DB_annual.csv')
# max_TS_DB_year.to_csv(outputpath + 'Wells55_GWSI_MAX_WLTS_DB_annual_comboID.csv')
max_TS_DB_year.to_csv(outputpath + 'Wells55_GWSI_MAX_WLTS_DB_annual_updated.csv')

# %% Creating totals for mapping
min_yr = 2000.0
mx_yr = 2022.0
measurement_buffer = 0
threshold = 23
# ds = WL_TS_DB_year.transpose()
ds = LEN_TS_DB_year.loc[:, ('depth', min_yr):('depth', mx_yr)]
# ds = ds.dropna(thresh=threshold) #sets a threshold
ds = ds.dropna(subset=ds.columns[measurement_buffer:-measurement_buffer])
# ds.info()
columndf = ds.transpose()
print("Number of wells", len(columndf.columns))
# %%
total_WL = ds.sum(axis=1)
total_WL = pd.DataFrame(total_WL)
total_WL = total_WL.reset_index()
total_WL = total_WL.rename(columns={"REGISTRY_ID": "Combo_ID",
                   0: "LEN"}, errors="raise")
total_WL = total_WL.set_index('Combo_ID')
# total_WL = total_WL['LEN'].astype(int, errors = 'raise')
total_WL
# total_WL.to_csv(outputpath+'numberWL_perwell_'+str(min_yr)+'-'+str(mx_yr)+'_comboID.csv')
total_WL.to_csv(outputpath+'numberWL_perwell_'+str(min_yr)+'-'+str(mx_yr)+'_updated'+'_bf'+str(measurement_buffer)+'.csv')
# total_WL.to_csv(outputpath+'numberWL_perwell_'+str(min_yr)+'-'+str(mx_yr)+'_updated'+'_thresh'+str(threshold)+'.csv')

# %%
>>>>>>> 00e8cedb276036a5996ed62c2a5cdef2568d0460
WL_TS_DB_month.to_csv(outputpath + 'Wells55_GWSI_WLTS_DB_monthly.csv')
WL_TS_DB_month_len.to_csv(outputpath + 'Wells55_GWSI_LEN_WLTS_DB_monthly.csv')
WL_TS_DB_month_max.to_csv(outputpath + 'Wells55_GWSI_MAX_WLTS_DB_monthly.csv')
WL_TS_DB_month_min.to_csv(outputpath + 'Wells55_GWSI_MIN_WLTS_DB_monthly.csv')

# %%
