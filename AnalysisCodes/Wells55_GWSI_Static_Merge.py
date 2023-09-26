### -- MASTER WELL DATABASE MERGE --
# Written by Danielle Tadych

# The purpose of this script is combine the GWSI and Wells55 databases into one
#   main database.

# Requirements:
# - Must download data from ADWR GIS database: https://gisdata2016-11-18t150447874z-azwater.opendata.arcgis.com/
#    > Well_Registry
#    > GWSI_App

# Data can also be found 

# WORKFLOW:
# - Create columns in each respective database specifying its origin
# - Find column they have in common
# - Merge based on that column
# - Make GWSI wells the overriding database and fill in the gaps with Wells55
# %%
import os
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import datetime
import seaborn as sns
import geopandas as gp

# %% 
# ----- Import the Data and Shapefiles with Geometries -----
shapedir = '../Data/Shapefiles'
outputpath = '../Data/Output_files/'

wellfilename = "Well_Registry_05032023.shp"
Wellfp = os.path.join(shapedir, wellfilename)
wells55shape = gp.read_file(Wellfp)

GWSI_fn = "GWSI_SITES.shp"
Wellfp = os.path.join(shapedir, GWSI_fn)
GWSIshape = gp.read_file(Wellfp)
# %% Take a look at the gwsi database
GWSIshape.info()
# %% Take a look at the wells55 database
wells55shape.info()
# %% Check that they're in the same coordinate system
print(GWSIshape.crs, wells55shape.crs)

# %%
# Preliminary Plot to check shapefiles
fig, ax = plt.subplots()
GWSIshape.plot(ax = ax, label="GWSI")
wells55shape.plot(ax = ax, label="Wells55")
ax.set_title("GWSI and Wells55 Preliminary Plot")
plt.legend()

# %% Making copies of the databases so I don't overright the originals
gwsi_gdf = GWSIshape
wells55_gdf = wells55shape

# ---- Adding Database Source Columns to both ----
wells55_gdf["Original_DB"] = 'Wells55'
gwsi_gdf["Original_DB"] = 'GWSI'
wells55_gdf.head()
# %%
gwsi_gdf.head()

# %% Fixing the date so that 1/1/70 in Wells55 is replaced with NAN
#   (this is because it is an xml error)
# https://stackoverflow.com/questions/29247712/how-to-replace-a-value-in-pandas-with-nan
wells55_gdf['INSTALLED'] = wells55_gdf['INSTALLED'].replace(['1970-01-01'], np.NaN)
wells55_gdf['INSTALLED'].unique()


# %% ---- Merging Both databases ----
# Merge wells55 'REGISTRY_I' with GWSI 'REG_ID'
# For wells in the GWSI database that do not have a Wells55 Registry ID we are going to use 'SITE_ID' 
 
# Changing REG_ID in GWSI to REGISTRY_I
gwsi_gdf.rename(columns={'REG_ID':'REGISTRY_I'}, inplace=True)

# %% combine registry ID and Site ID so in timeseries graphs every well has an ID
gwsi_gdf['Combo_ID'] = gwsi_gdf.REGISTRY_I.combine_first(gwsi_gdf.SITE_ID)
gwsi_gdf.info()

# %%
wells55_gdf['Combo_ID'] = wells55_gdf['REGISTRY_I']
wells55_gdf.info()

# -- Stop here if you want to filter the wells55 database for specific uses
# -- skip to line 98
# %%
Wells55_GWSI_MasterDB = pd.merge(gwsi_gdf, wells55_gdf, suffixes=['_gwsi','_wells55'], how="outer", 
                                          on=['Combo_ID',"REGISTRY_I", 
                                            #   'WELL_TYPE_', 
                                              'WELL_DEPTH', 'geometry', 'Original_DB'])
print(Wells55_GWSI_MasterDB.info())

# %% Looking for duplicates
print(Wells55_GWSI_MasterDB[Wells55_GWSI_MasterDB['Combo_ID'].duplicated()])

# %%
test = Wells55_GWSI_MasterDB.groupby(['Combo_ID']).first()
test.info()

# %%
test = test.reset_index()
print(test[test['Combo_ID'] == '921247'].loc[:,['geometry','Original_DB']])

# %% Check for duplicates again
print(test[test['Combo_ID'].duplicated()])

# %% hooray! No duplicates, now to re-create those files
Wells55_GWSI_MasterDB = test

# %%
# Now plot the new master db to see if it looks okay
fig, ax = plt.subplots()
Wells55_GWSI_MasterDB.plot(ax=ax, label="Master Database")
ax.set_title("Check the merged database")
plt.legend()
# Optional save
# plt.savefig('../Data/Output_files/{0}.png'.format(type), bbox_inches='tight')

# %%
# Export the database without duplicates
Wells55_GWSI_MasterDB.to_csv('../Data/Output_files/Master_ADWR_database_noduplicates.csv')
Wells55_GWSI_MasterDB.to_file('../Data/Output_files/Master_ADWR_database_noduplicates.shp')
# %%
# --- Filter wells55 for only wells that have been drilled ---
# First, deleting the cancelled wells
wells55_nocanc = wells55_gdf[wells55_gdf.WELL_CANCE != 'Y']
wells55_nocanc.info()

#%%
# No classified mineral or exploratory wells
wells55_nomin = wells55_nocanc[wells55_nocanc.WELL_TYPE_ != 'OTHER']
wells55_nomin.info()

# %% Now filter by drill log exisence (not just filed) so we can see which water wells are actually a thing
wells55_water = wells55_nomin[wells55_nomin['DRILL_LOG'].notna()]
wells55_water.info()

#df = df[df['EPS'].notna()]

# %% Export both the non-cancelled wells and the water wells
#  First the non-cancelled wells
Wells55_GWSI_MasterDB = pd.merge(gwsi_gdf, wells55_nocanc, suffixes=['_gwsi','_wells55'], how="outer", 
                                          on=['Combo_ID',"REGISTRY_I", 'WELL_DEPTH', 'geometry', 'Original_DB'])

print(Wells55_GWSI_MasterDB.info())
#%%
Wells55_GWSI_MasterDB = Wells55_GWSI_MasterDB.groupby(['Combo_ID']).first()
print(Wells55_GWSI_MasterDB.info())

# %%
Wells55_GWSI_MasterDB.to_csv('../Data/Output_files/Master_ADWR_database_nocancelled.csv')
Wells55_GWSI_MasterDB.to_file('../Data/Output_files/Master_ADWR_database_nocancelled.shp')

# %%
#  - Only water wells 
Wells55_GWSI_MasterDB = pd.merge(gwsi_gdf, wells55_water, suffixes=['_gwsi','_wells55'], how="outer", 
                                          on=['Combo_ID',"REGISTRY_I", 'WELL_DEPTH', 'geometry', 'Original_DB'])

print(Wells55_GWSI_MasterDB.info())

# %%
Wells55_GWSI_MasterDB = Wells55_GWSI_MasterDB.groupby(['Combo_ID']).first()
print(Wells55_GWSI_MasterDB.info())

# %%
Wells55_GWSI_MasterDB.to_csv('../Data/Output_files/Master_ADWR_database_water.csv')
Wells55_GWSI_MasterDB.to_file('../Data/Output_files/Master_ADWR_database_water.shp')
