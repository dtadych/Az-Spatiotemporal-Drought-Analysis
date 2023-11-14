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
# Need to create a combo ID column
gwsi_wl['Combo_ID'] = gwsi_wl.SITE_WELL_REG_ID.combine_first(gwsi_wl.wellid)
gwsi_wl.rename(columns={'Combo_ID':'REGISTRY_ID'}, inplace=True)

# Need to make sure the columns are the same Dtype
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

# %% --- Summarizing the data by date ---
# Extract the year from the date column and create a new column year
combo['year'] = pd.DatetimeIndex(combo.date).year
combo.head()

# %% Creating pivot tables
WL_TS_DB_year = pd.pivot_table(combo, index=["REGISTRY_ID"], columns=["year"], values=["depth"], dropna=False, aggfunc=np.mean)
LEN_TS_DB_year = pd.pivot_table(combo, index=["REGISTRY_ID"], columns=["year"], values=["depth"], dropna=False, aggfunc=len)
max_TS_DB_year = pd.pivot_table(combo, index=["REGISTRY_ID"], columns=["year"], values=["depth"], dropna=False, aggfunc=np.max)
min_TS_DB_year = pd.pivot_table(combo, index=["REGISTRY_ID"], columns=["year"], values=["depth"], dropna=False, aggfunc=np.min)

#%%
WL_TS_DB_year.index.name = None
WL_TS_DB_year.head()
# %%
test = WL_TS_DB_year.mean()
test.plot()

# %% 
stats = WL_TS_DB_year.describe()
stats = stats.transpose()
stats2 = stats[['mean','25%','50%','75%']]
stats2[119:157].plot()

# %%
narrowedstats = stats[110:158]
narrowedstats
# %%
narrowedstats.to_csv(outputpath+"state_average_WL_updated.csv")

# %%
# Export both yearly summary data and monthly into csv
WL_TS_DB_year.to_csv(outputpath + 'Wells55_GWSI_WLTS_DB_annual_updated.csv')
# %%
LEN_TS_DB_year.to_csv(outputpath + 'Wells55_GWSI_LEN_WLTS_DB_annual_updated.csv')

# %%
max_TS_DB_year.to_csv(outputpath + 'Wells55_GWSI_MAX_WLTS_DB_annual_updated.csv')

# %%
min_TS_DB_year.to_csv(outputpath + 'Wells55_GWSI_MIN_WLTS_DB_annual_updated.csv')

# %% Creating totals for mapping
min_yr = 2000.0
mx_yr = 2022.0
threshold = 15
# ds = WL_TS_DB_year.transpose()
ds = LEN_TS_DB_year.loc[:, ('depth', min_yr):('depth', mx_yr)]
ds = ds.dropna(thresh=threshold) #sets a threshold
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
total_WL.to_csv(outputpath+'numberWL_perwell_'+str(min_yr)+'-'+str(mx_yr)+'_updated'+')

# %%
ds = ds.reset_index()
# Create a list of narrowed Wells
well_list = ds[['REGISTRY_ID']]

# Converting to int
well_list['REGISTRY_ID'] = well_list['REGISTRY_ID'].astype(int, errors = 'raise')
# %%
# Formatting timeseries
# WL_TS = WL_TS_DB_year.copy()
# WL_TS = max_TS_DB_year.copy()
WL_TS = min_TS_DB_year.copy()
WL_TS.reset_index(inplace=True)
# WL_TS

# Add list to the water level database
narrowedWL_DB = WL_TS.merge(well_list, how="inner")
narrowedWL_DB.info()

# set index back go REGISTRY_ID
narrowedWL_DB.set_index('REGISTRY_ID', inplace=True)
# %%
# narrowedWL_DB.to_csv(outputpath+'Wells55_GWSI_WLTS_DB_annual_updated_thresh'+str(threshold)+'.csv')
# narrowedWL_DB.to_csv(outputpath+'Wells55_GWSI_MAX_WLTS_DB_annual_updated_thresh'+str(threshold)+'.csv')
narrowedWL_DB.to_csv(outputpath+'Wells55_GWSI_MIN_WLTS_DB_annual_updated_thresh'+str(threshold)+'.csv')
print('Done.')
# %% ====  Summarize by Monthly now ====
combo_monthly2 = combo.copy()
combo_monthly2.reset_index(inplace=True)
combo_monthly2.date = combo.date.dt.strftime('%Y-%m')
combo_monthly2.info()
# %%
# Average
WL_TS_DB_month = pd.pivot_table(combo_monthly2, index=["REGISTRY_ID"], columns=["date"], values=["depth"], dropna=False, aggfunc='mean')
WL_TS_DB_month.index.name = None

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
WL_TS_DB_month.to_csv(outputpath + 'Wells55_GWSI_WLTS_DB_monthly.csv')
WL_TS_DB_month_len.to_csv(outputpath + 'Wells55_GWSI_LEN_WLTS_DB_monthly.csv')
WL_TS_DB_month_max.to_csv(outputpath + 'Wells55_GWSI_MAX_WLTS_DB_monthly.csv')
WL_TS_DB_month_min.to_csv(outputpath + 'Wells55_GWSI_MIN_WLTS_DB_monthly.csv')

# %%
