# ------- Well Count Analysis and Graphs -------
# Danielle Tadych,

# The purpose of this code is to analyze well counts and depths with respect to groundwater regulation
#   and access to surface water

# More specifically,
# - Well Density (cumulative number of wells) over time
# - max screen depth over time (Casing_DEP vs Installed)
# - Number of new wells installed over time

# %% Load the packages
from cProfile import label
from operator import ge
from optparse import Values
import os
from geopandas.tools.sjoin import sjoin
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from matplotlib.colors import ListedColormap
import datetime
from matplotlib.transforms import Bbox
import seaborn as sns
import numpy as np
import pandas as pd
from shapely.geometry import box
import geopandas as gp
import scipy.stats as sp
import Custom_functions as cf
# %% Assign Data paths
datapath_web = 'https://datacommons.cyverse.org/browse/iplant/home/shared/commons_repo/curated/Tadych_AzGroundwaterSpatialAnalysis_Aug2023/Data/'
outputpath_web = 'https://datacommons.cyverse.org/browse/iplant/home/shared/commons_repo/curated/Tadych_AzGroundwaterSpatialAnalysis_Aug2023/Data/Output_files/'
shapepath_web = 'https://datacommons.cyverse.org/browse/iplant/home/shared/commons_repo/curated/Tadych_AzGroundwaterSpatialAnalysis_Aug2023/Data/Shapefiles/'

datapath = '../Data'
outputpath = '../Data/Output_files/'
shapepath = '../Data/Shapefiles/'

# %% Reading in the data
filename = 'Final_Static_geodatabase_waterwells.csv'
# filepath = os.path.join(outputpath_web, filename)
filepath = os.path.join(outputpath, filename)
print(filepath)
static_geo2 = pd.read_csv(filepath 
                          ,parse_dates=['INSTALLED']
                          )
static_geo2

# %% Only run this if date parser didn't work
static_geo2['INSTALLED'] = pd.to_datetime(static_geo2['INSTALLED'])
static_geo2['INSTALLED'].describe()
# %%
static_geo2['In_year'] = static_geo2['INSTALLED'].dt.year
static_geo2['In_year'].describe()

# %% Checking information about static_geo2

static_na = static_geo2[static_geo2['Regulation'].isna()]
print(static_na)

# %%
static_na['GEO_Region'].unique()

# %% 
Well_Depth = static_geo2[['WELL_DEPTH', 'INSTALLED', 'Combo_ID', 'In_year','GEOREGI_NU', 'Regulation', 'Water_CAT','WELL_TYPE_']]
#%%
Well_Depth

#%%
# Well_Depth.to_csv('../MergedData/Output_files/Final_WellDepth.csv')
#%%
# Well_Depth = pd.pivot_table(static_geo2, index=["In_year"], columns=["GEOREGI_NU"], values=["WELL_DEPTH"], dropna=False, aggfunc=np.mean)
Well_Depth = pd.pivot_table(static_geo2, index=["In_year"], columns=["Regulation"], values=["WELL_DEPTH"], dropna=False, aggfunc=np.mean)
state_depth = pd.pivot_table(static_geo2, index=["In_year"], values=["WELL_DEPTH"], dropna=False, aggfunc=np.mean)
Well_Depth.describe()

# %% Set shallow and drilling depths
shallow = 200
deep = 500

# %%
wd1 = static_geo2[(static_geo2["WELL_DEPTH"] > deep)]
wd2 = static_geo2[(static_geo2["WELL_DEPTH"] <= deep) & (static_geo2["WELL_DEPTH"] >= shallow)]
wd3 = static_geo2[(static_geo2["WELL_DEPTH"] < shallow)]

# %%
wd1 = wd1.sort_values(by=['GEOREGI_NU'])
wd2 = wd2.sort_values(by=['GEOREGI_NU'])
wd3 = wd3.sort_values(by=['GEOREGI_NU'])

#%%

wdc1_reg = pd.pivot_table(wd1, index=["In_year"], columns=["Regulation"], values=['WELL_DEPTH'], dropna=False, aggfunc=len)
wdc2_reg = pd.pivot_table(wd2, index=["In_year"], columns=["Regulation"], values=['WELL_DEPTH'], dropna=False, aggfunc=len)
wdc3_reg = pd.pivot_table(wd3, index=["In_year"], columns=["Regulation"], values=['WELL_DEPTH'], dropna=False, aggfunc=len)

wdc1_wc = pd.pivot_table(wd1, index=["In_year"], columns=["Water_CAT"], values=['WELL_DEPTH'], dropna=False, aggfunc=len)
wdc2_wc = pd.pivot_table(wd2, index=["In_year"], columns=["Water_CAT"], values=['WELL_DEPTH'], dropna=False, aggfunc=len)
wdc3_wc = pd.pivot_table(wd3, index=["In_year"], columns=["Water_CAT"], values=['WELL_DEPTH'], dropna=False, aggfunc=len)

# %%
st_wdc1 = pd.pivot_table(wd1, index=["In_year"], values=['WELL_DEPTH'], dropna=False, aggfunc=len)
st_wdc2 = pd.pivot_table(wd2, index=["In_year"], values=['WELL_DEPTH'], dropna=False, aggfunc=len)
st_wdc3 = pd.pivot_table(wd3, index=["In_year"], values=['WELL_DEPTH'], dropna=False, aggfunc=len)
st_wdc1

# %% Exempt/non-exempt added
wdc1_reg_ex = pd.pivot_table(wd1, index=["In_year"], columns=["Regulation",'WELL_TYPE_'], values=['WELL_DEPTH'], dropna=False, aggfunc=len)
wdc2_reg_ex = pd.pivot_table(wd2, index=["In_year"], columns=["Regulation",'WELL_TYPE_'], values=['WELL_DEPTH'], dropna=False, aggfunc=len)
wdc3_reg_ex = pd.pivot_table(wd3, index=["In_year"], columns=["Regulation",'WELL_TYPE_'], values=['WELL_DEPTH'], dropna=False, aggfunc=len)

wdc1_wc_ex = pd.pivot_table(wd1, index=["In_year"], columns=["Water_CAT",'WELL_TYPE_'], values=['WELL_DEPTH'], dropna=False, aggfunc=len)
wdc2_wc_ex = pd.pivot_table(wd2, index=["In_year"], columns=["Water_CAT",'WELL_TYPE_'], values=['WELL_DEPTH'], dropna=False, aggfunc=len)
wdc3_wc_ex = pd.pivot_table(wd3, index=["In_year"], columns=["Water_CAT",'WELL_TYPE_'], values=['WELL_DEPTH'], dropna=False, aggfunc=len)

#%% Exporting the Depth categories
wdc1_reg.to_csv(outputpath+'Final_Welldepth_regulation' + str(deep) + 'plus.csv')
wdc2_reg.to_csv(outputpath+'Final_Welldepth_regulation' + str(shallow) + 'to' + str(deep) + '.csv')
wdc3_reg.to_csv(outputpath+'Final_Welldepth_regulation' + str(shallow) + 'minus.csv')

wdc1_wc.to_csv(outputpath+'Final_Welldepth_sw' + str(deep) + 'plus.csv')
wdc2_wc.to_csv(outputpath+'Final_Welldepth_sw' + str(shallow) + 'to' + str(deep) + '.csv')
wdc3_wc.to_csv(outputpath+'Final_Welldepth_sw' + str(shallow) + 'minus.csv')

wdc1_reg_ex.to_csv(outputpath+'Final_Welldepth_regulation_exemptstatus' + str(deep) + 'plus.csv')
wdc2_reg_ex.to_csv(outputpath+'Final_Welldepth_regulation_exemptstatus' + str(shallow) + 'to' + str(deep) + '.csv')
wdc3_reg_ex.to_csv(outputpath+'Final_Welldepth_regulation_exemptstatus' + str(shallow) + 'minus.csv')

wdc1_wc_ex.to_csv(outputpath+'Final_Welldepth_sw_exemptstatus' + str(deep) + 'plus.csv')
wdc2_wc_ex.to_csv(outputpath+'Final_Welldepth_sw_exemptstatus' + str(shallow) + 'to' + str(deep) + '.csv')
wdc3_wc_ex.to_csv(outputpath+'Final_Welldepth_sw_exemptstatus' + str(shallow) + 'minus.csv')

# %% Now just new wells for general and density calculations
new_wells = pd.pivot_table(static_geo2, index=["In_year"], columns=["GEOREGI_NU"], values=["INSTALLED"], dropna=False, aggfunc=len)
new_wells_reg = pd.pivot_table(static_geo2, index=["In_year"], columns=["Regulation"], values=["INSTALLED"], dropna=False, aggfunc=len)
new_wells_watercat = pd.pivot_table(static_geo2, index=["In_year"], columns=["Water_CAT"], values=["INSTALLED"], dropna=False, aggfunc=len)
new_wells_reg
new_wells_watercat
# %%
new_wells.to_csv(outputpath+'Final_NewWells.csv')
new_wells_reg.to_csv(outputpath+'Final_NewWells_regulation.csv')
new_wells_watercat.to_csv(outputpath+'Final_NewWells_watercat.csv')

# %% ---- Fancier Analyses ----
# Calculating well densities
new_wells2 = pd.read_csv(outputpath+'/Final_NewWells.csv',
                        header=1, index_col=0)
new_wells2
#%%
new_wells2 = new_wells2.reset_index()
new_wells2 = new_wells2.iloc[1:, :]
new_wells2 = new_wells2.set_index('GEOREGI_NU')
new_wells2

# %%
new_wells_reg2 = pd.read_csv(outputpath+'Final_NewWells_regulation.csv',
                        header=2,
                        names = ['R','Res','U']
                        , index_col=0)
new_wells_reg2 = new_wells_reg2.iloc[1:, :]
new_wells_reg2

# %%
del new_wells_reg2['Res']
new_wells_reg2

# %%
summies = new_wells_reg2.sum()

# %%
summies = 36474 + 48027
summies
# %%
new_wells_watercat2 = pd.read_csv(outputpath+'Final_NewWells_watercat.csv',
                        header=1,
                        names = ['CAP','GW','Mix','No_CAP','Res','SW']
                        , index_col=0)
new_wells_watercat2 = new_wells_watercat2.iloc[1:, :]
new_wells_watercat2

# %%
summies = new_wells_watercat2.sum()
summies
# %%
filename_georeg = 'georeg_reproject_fixed.shp'
filepath = os.path.join(shapepath+filename_georeg)
georeg = gp.read_file(filepath)
# %% Calculate the region area
georeg = georeg.to_crs(epsg=3857)
georeg.crs

# %% converting to km
georeg['area'] = georeg.geometry.area/10**6
georeg

# %%
georeg2 = pd.DataFrame(georeg)
#%%
georeg_area = georeg2[['GEOREGI_NU','area']]
georeg_area.info()

# %%
georeg_area = georeg_area.set_index('GEOREGI_NU')
georeg_area = georeg_area.transpose()
georeg_area

# %% Area for Regulation
georeg_area_reg = pd.pivot_table(georeg2, columns=["Regulation"], values=["area"], dropna=False, aggfunc=np.sum)
georeg_area_reg

# %%
georeg_area_watercat = pd.pivot_table(georeg2, columns=["Water_CAT"], values=["area"], dropna=False, aggfunc=np.sum)
# del georeg_area_watercat['Res']
georeg_area_watercat
# %% Densities for new wells
well_densities = new_wells2/georeg_area.values[0,:]
well_densities

# %% Densities for regulated regions
well_densities_reg = new_wells_reg2/georeg_area_reg.values[0,:]
well_densities_reg

# %% Densities for SW
well_densities_watercat = new_wells_watercat2/georeg_area_watercat.values[0,:]
well_densities_watercat.sum()

# %% By depth for regulated or water category, depending on what I turned on above
dens_wdc1_reg= wdc1_reg/georeg_area_reg.values[0,:]
dens_wdc2_reg= wdc2_reg/georeg_area_reg.values[0,:]
dens_wdc3_reg= wdc3_reg/georeg_area_reg.values[0,:]

dens_wdc1_wc= wdc1_wc/georeg_area_watercat.values[0,:]
dens_wdc2_wc= wdc2_wc/georeg_area_watercat.values[0,:]
dens_wdc3_wc= wdc3_wc/georeg_area_watercat.values[0,:]

dens_wdc1_wc
print(dens_wdc1_wc.sum())

# %% Well densities but with exemption status
df = pd.DataFrame(wdc1_reg_ex.sum())
df = df.transpose()
df
# %%
georeg_area_reg_dens = pd.DataFrame(df,index=df.index)
# georeg_area_reg_dens = georeg_area_reg_dens.transpose()
# %%
georeg_area_reg_dens.iloc[0,0:3] = georeg_area_reg.iloc[0,0]
georeg_area_reg_dens.iloc[0,3:6] = georeg_area_reg.iloc[0,1]
georeg_area_reg_dens.iloc[0,6:9] = georeg_area_reg.iloc[0,2]
georeg_area_reg_dens

# %% 
dens_wdc1_reg_ex= wdc1_reg_ex/georeg_area_reg_dens.values[0,:]
dens_wdc2_reg_ex= wdc2_reg_ex/georeg_area_reg_dens.values[0,:]
dens_wdc3_reg_ex= wdc3_reg_ex/georeg_area_reg_dens.values[0,:]
dens_wdc3_reg_ex.sum()
# %%
df = pd.DataFrame(wdc1_wc_ex.sum())
df = df.transpose()
df
# %%
georeg_area_watercat_dens = pd.DataFrame(df,index=df.index)
georeg_area_watercat_dens.iloc[0,0:3] = georeg_area_watercat.iloc[0,0]
georeg_area_watercat_dens.iloc[0,3:6] = georeg_area_watercat.iloc[0,1]
georeg_area_watercat_dens.iloc[0,6:9] = georeg_area_watercat.iloc[0,2]
georeg_area_watercat_dens.iloc[0,9:12] = georeg_area_watercat.iloc[0,3]
georeg_area_watercat_dens.iloc[0,12:15] = georeg_area_watercat.iloc[0,4]
georeg_area_watercat_dens.iloc[0,15:18] = georeg_area_watercat.iloc[0,5]
georeg_area_watercat_dens

# %%
dens_wdc1_wc_ex= wdc1_wc_ex/georeg_area_watercat_dens.values[0,:]
dens_wdc2_wc_ex= wdc2_wc_ex/georeg_area_watercat_dens.values[0,:]
dens_wdc3_wc_ex= wdc3_wc_ex/georeg_area_watercat_dens.values[0,:]
dens_wdc3_wc_ex.info()

#%% Exporting the well densities
dens_wdc1_reg.to_csv(outputpath+'FinalDensities_Welldepth_regulation' + str(deep) + 'plus.csv')
dens_wdc2_reg.to_csv(outputpath+'FinalDensities_Welldepth_regulation' + str(shallow) + 'to' + str(deep) + '.csv')
dens_wdc3_reg.to_csv(outputpath+'FinalDensities_Welldepth_regulation' + str(shallow) + 'minus.csv')

dens_wdc1_wc.to_csv(outputpath+'FinalDensities_Welldepth_sw' + str(deep) + 'plus.csv')
dens_wdc2_wc.to_csv(outputpath+'FinalDensities_Welldepth_sw' + str(shallow) + 'to' + str(deep) + '.csv')
dens_wdc3_wc.to_csv(outputpath+'FinalDensities_Welldepth_sw' + str(shallow) + 'minus.csv')

dens_wdc1_reg_ex.to_csv(outputpath+'FinalDensities_Welldepth_regulation_exemptstatus' + str(deep) + 'plus.csv')
dens_wdc2_reg_ex.to_csv(outputpath+'FinalDensities_Welldepth_regulation_exemptstatus' + str(shallow) + 'to' + str(deep) + '.csv')
dens_wdc3_reg_ex.to_csv(outputpath+'FinalDensities_Welldepth_regulation_exemptstatus' + str(shallow) + 'minus.csv')

dens_wdc1_wc_ex.to_csv(outputpath+'FinalDensities_Welldepth_sw_exemptstatus' + str(deep) + 'plus.csv')
dens_wdc2_wc_ex.to_csv(outputpath+'FinalDensities_Welldepth_sw_exemptstatus' + str(shallow) + 'to' + str(deep) + '.csv')
dens_wdc3_wc_ex.to_csv(outputpath+'FinalDensities_Welldepth_sw_exemptstatus' + str(shallow) + 'minus.csv')

# %%
well_densities1 = well_densities.reset_index()
well_densities1['GEOREGI_NU'] = pd.to_numeric(well_densities1['GEOREGI_NU'])
well_densities1['GEOREGI_NU'] = well_densities1['GEOREGI_NU'].astype(int)
well_densities1.set_index('GEOREGI_NU', inplace=True)
well_densities1.info()

# %% -- Linear regression --
# Change years as appropriate

# For Depth to Water by Regulation
# == For Shallow ==
ds = wdc3_reg
dt = "Regulation, Shallow New Wells"
min = 1975
mx = 2020
betterlabels = ['Regulated','Unregulated'] 
cf.linearregress(ds,dt,min,mx,betterlabels)
# == For Midrange ==
ds = wdc2_reg
dt = "Regulation, Midrange New Wells"
min = 1975
mx = 2020
betterlabels = ['Regulated','Unregulated'] 
cf.linearregress(ds,dt,min,mx,betterlabels)

# == For Deep ==
ds = wdc1_reg
dt = "Regulation, Deep New Wells"
min = 1975
mx = 2020
betterlabels = ['Regulated','Unregulated'] 
cf.linearregress(ds,dt,min,mx,betterlabels)

# %% For Depth to Water by SW Access
# == For Shallow ==
ds = wdc3_wc
dt = "Access to Surface Water, Shallow New Wells"
min = 1975
mx = 2020
betterlabels = ['Recieves CAP (Regulated)'
                ,'GW Dominated (Regulated)'
                ,'Surface Water Dominated'
                ,'GW Dominated'
                ,'Mixed Source']
cf.linearregress(ds,dt,min,mx,betterlabels)
# == For Midrange ==
ds = wdc2_wc
dt = "Access to Surface Water, Midrange New Wells"
min = 1975
mx = 2020
betterlabels = ['Recieves CAP (Regulated)'
                ,'GW Dominated (Regulated)'
                ,'Surface Water Dominated'
                ,'GW Dominated'
                ,'Mixed Source']
cf.linearregress(ds,dt,min,mx,betterlabels)

# == For Deep ==
ds = wdc1_wc
dt = "Access to Surface Water, Deep New Wells"
min = 1975
mx = 2020
betterlabels = ['Recieves CAP (Regulated)'
                ,'GW Dominated (Regulated)'
                ,'Surface Water Dominated'
                ,'GW Dominated'
                ,'Mixed Source']
cf.linearregress(ds,dt,min,mx,betterlabels)