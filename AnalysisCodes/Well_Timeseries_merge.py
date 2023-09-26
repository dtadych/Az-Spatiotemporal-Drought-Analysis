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
wl_data2.to_csv(outputpath+'wl_data2.csv')

# %% 
# ----- Import the Data and Shapefiles with Geometries -----
# Read in Wells 55 Data
# This is a file with water levels from ADWR which has been joined with another ADWR file with variables
filename = 'Well_Registry_05032023.csv'
filepath = os.path.join(datapath+'Wells55', filename)
print(filepath)

wells55 = pd.read_csv(filepath)
pd.options.display.float_format = '{:.2f}'.format
print(wells55.info())

# Read in GWSI collated water level data
filename = 'wl_data2.csv'
filepath = os.path.join(outputpath, filename)
print(filepath)

wl_data2 = pd.read_csv(filepath)
pd.options.display.float_format = '{:.2f}'.format
print(wl_data2.info())

#%% ---- Making dataframes with Date as the index ---
# Make dataframes with Columns for GWSI and Wells55, respectively
# Following this method: https://stackoverflow.com/questions/32215024/merging-time-series-data-by-timestamp-using-numpy-pandas
# Confirmed through the variables list that "depth" in wldata2 and "WATER_LEVE" are both depth to water below land surface in feet

# gwsi_wl = wl_data2[["date","wellid","SITE_WELL_REG_ID","depth"]].copy()
# gwsi_wl.info()

gwsi_wl = wl_data2[["date","SITE_WELL_REG_ID","depth"]].copy()
gwsi_wl.info()

wells55_wl = wells55[["INSTALLED", "REGISTRY_ID", "WATER_LEVEL"]].copy()
wells55_wl.info()

# %%
# Need to add an original database column
wells55_wl["Original_DB"] = 'Wells55'
gwsi_wl["Original_DB"] = 'GWSI'
wells55_wl.head()

# %%
wells55_wl.rename(columns = {'INSTALLED':'date','WATER_LEVEL':'depth'}, inplace=True)
# %% Need to create a combo ID column
# gwsi_wl['Combo_ID'] = gwsi_wl.SITE_WELL_REG_ID.combine_first(gwsi_wl.wellid)
# gwsi_wl.info()


# gwsi_wl.rename(columns={'Combo_ID':'REGISTRY_ID'}, inplace=True)
gwsi_wl.rename(columns={'SITE_WELL_REG_ID':'REGISTRY_ID'}, inplace=True)
gwsi_wl.info()

# %% Need to make sure the columns are the same Dtype
gwsi_wl.REGISTRY_ID = gwsi_wl.REGISTRY_ID.astype('int64')
gwsi_wl.info()
# %%
wells55_wl['date'] = pd.to_datetime(wells55_wl.date)
wells55_wl.info()
# %%
wells55_wl['date'] = wells55_wl['date'].dt.tz_localize(None)
wells55_wl.info()

# %%
gwsi_wl.date = pd.to_datetime(gwsi_wl.date)
gwsi_wl.info()
#%%
#combo = gwsi_wl.join(wells55_wl, how='outer')
#combo

combo = wells55_wl.merge(gwsi_wl, suffixes=['_wells55','_gwsi'], how="outer" 
                                          ,on=["REGISTRY_ID", 'date', 'Original_DB', 'depth']
                                          )
combo.info()

# # %% This block of code takes a really long time to run so I would suggest skipping
WL_TS_DB = pd.pivot_table(combo, index=["REGISTRY_ID"], columns="date", values="depth")
# WL_TS_DB.head()
# # Export data into a csv
# WL_TS_DB.to_csv(outputpath + 'Wells55_GWSI_WLTS_DB.csv')

# %% --- Summarizing the data by date ---
# Extract the year from the date column and create a new column year
combo['year'] = pd.DatetimeIndex(combo.date).year
combo.head()

# %%
WL_TS_DB_year = pd.pivot_table(combo, index=["REGISTRY_ID"], columns=["year"], values=["depth"], dropna=False, aggfunc=np.mean)
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
narrowedstats.to_csv(outputpath+"state_average_WL.csv")
# %%
max = stats['max']
max[119:157].plot()
# %% Exporting data
WL_TS_DB_1980.to_csv(outputpath + 'comboDB_WL_1980.csv')
WL_TS_DB_2020.to_csv(outputpath + 'comboDB_WL_2020.csv')
# %%  Seeing if things work
fig, ax = plt.subplots()
ax.plot(WL_TS_DB_year.iloc[:,155])
ax.set(title='WL in 1980', xlabel='Registry ID', ylabel='Water Level (feet)')
ax.grid()
plt.show
# %%
WL_TS_DB_year.info()
# %%
# Summarize by Monthly now
combo_monthly2 = combo
combo_monthly2.reset_index(inplace=True)
combo_monthly2.date = combo.date.dt.strftime('%Y-%m')
combo_monthly2.info()
# %%
WL_TS_DB_month = pd.pivot_table(combo_monthly2, index=["REGISTRY_I"], columns=["date"], values=["depth"], dropna=False, aggfunc='mean')
WL_TS_DB_month.info()
WL_TS_DB_month.index.name = None
WL_TS_DB_month.head()
# %%
# Export both yearly summary data and monthly into csv
WL_TS_DB_year.to_csv(outputpath + 'Wells55_GWSI_WLTS_DB_annual.csv')
# %%
WL_TS_DB_month.to_csv(outputpath + 'Wells55_GWSI_WLTS_DB_monthly.csv')
# %%
# Making a for loop for years
# Going to re-read in the combined database so we can get rid of that first row
filename = 'Wells55_GWSI_WLTS_DB_annual.csv'
filepath = os.path.join(outputpath, filename)
print(filepath)

annual_db = pd.read_csv(filepath, header=1, index_col=0)
pd.options.display.float_format = '{:.2f}'.format
annual_db.head()

# %%
columns = annual_db.columns.values.tolist()
print(columns)

# %%
for i in columns:
    WL_TS = 'WL_TS_DB_' + i
    print(WL_TS)
    f = pd.DataFrame(annual_db[i])
    f.index.name = None
    f.columns = ["depth"]
    print(f.describe())
#    print(outputpath + 'YearlyTS/' + WL_TS + '.csv')
    f.to_csv(outputpath + 'YearlyTS/' + WL_TS + '.csv')

# %%
# Making a for loop for months
# Going to re-read in the combined database so we can get rid of that first row
filename = 'Wells55_GWSI_WLTS_DB_monthly.csv'
filepath = os.path.join(outputpath, filename)
print(filepath)

monthly_db = pd.read_csv(filepath, header=1, index_col=0)
pd.options.display.float_format = '{:.2f}'.format
monthly_db.head()

# %%
columns = monthly_db.columns.values.tolist()
print(columns)

# %% Testing before the for loop
i = '2007-05'
WL_TS = 'WL_TS_DB_' + i
print(WL_TS)
f = pd.DataFrame(monthly_db[i])
f.index.name = None
f.columns = ["depth"]
print(f.describe())
print(outputpath + 'MonthlyTS/' + WL_TS + '.csv')


# %%
for i in columns:
    WL_TS = 'WL_TS_DB_' + i
    print(WL_TS)
    f = pd.DataFrame(monthly_db[i])
    f.index.name = None
    f.columns = ["depth"]
    print(f.describe())
#    print(outputpath + 'MonthlyTS/' + WL_TS + '.csv')
    f.to_csv(outputpath + 'MonthlyTS/' + WL_TS + '.csv')
