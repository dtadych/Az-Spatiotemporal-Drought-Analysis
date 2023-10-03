# ----- All Paper Graphs except graphics from QGIS -----
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
from scipy.stats import kendalltau, pearsonr, spearmanr
import pymannkendall as mk
import Custom_functions

# === Assign Data paths ===

# This is for accessing our data on Cyverse
# 
datapath_web = 'https://datacommons.cyverse.org/browse/iplant/home/shared/commons_repo/curated/Tadych_AzGroundwaterSpatialAnalysis_Aug2023/Data/'
outputpath_web = 'https://datacommons.cyverse.org/browse/iplant/home/shared/commons_repo/curated/Tadych_AzGroundwaterSpatialAnalysis_Aug2023/Data/Output_files/'
shapepath_web = 'https://datacommons.cyverse.org/browse/iplant/home/shared/commons_repo/curated/Tadych_AzGroundwaterSpatialAnalysis_Aug2023/Data/Shapefiles/'

# This is if you created your own database
datapath_local = '../Data'
inputpath_local = '../Data/Input_files/'
outputpath_local = '../Data/Output_files/'
shapepath_local = '../Data/Shapefiles/'
figurepath = '../Data/Figures/'

# Change this based on whether you're running off local or web
# Cyverse:
# outputpath = outputpath_web
# shapepath = shapepath_web

# Local: 
outputpath = outputpath_local
shapepath = shapepath_local
inputpath = inputpath_local
# %% Read in the data

# %% Importing Water Level Values
# For statewide
filename_ts = 'Wells55_GWSI_WLTS_DB_annual.csv'
filepath = os.path.join(inputpath, filename_ts)
print(filepath)
annual_db = pd.read_csv(filepath, header=1, index_col=0)
annual_db.head()

# For regulation
filepath = inputpath+'/Waterlevels_Regulation.csv'
# filepath = '../Data/Output_files/Waterlevels_Regulation.csv'
cat_wl2_reg = pd.read_csv(filepath, index_col=0)
cat_wl2_reg.head()

# For Access to SW
filepath = inputpath+'/Waterlevels_AccesstoSW.csv'
# filepath = '../Data/Output_files/Waterlevels_AccesstoSW_GWlumped.csv'
cat_wl2_SW = pd.read_csv(filepath, index_col=0)
cat_wl2_SW.head()

# For georegion number
filepath = inputpath+'Waterlevels_georegions.csv'
# filepath = '../Data/Output_files/Waterlevels_georegions.csv'
cat_wl2_georeg = pd.read_csv(filepath, index_col=0)
# cat_wl2_georeg.head()

# %% Importing GRACE analyses
filepath = filepath = inputpath+'grace_stateavg_yearly.csv'
# filepath = outputpath_local+'gracse_remapped_yearly.csv'
grace_yearly = pd.read_csv(filepath, index_col=0)
grace_yearly = grace_yearly[:-1]

# Reading in the shapefile - note, figure 2 is created through QGIS
filename_georeg = 'georeg_reproject_fixed.shp'
filepath = os.path.join(shapepath, filename_georeg)
# %% creating weights for the GW dominated categories
georeg = gp.read_file(filepath)
georeg.plot(cmap='viridis')
georeg['area'] = georeg.geometry.area/10**6
georeg_wcweights = georeg.groupby(['Water_CAT']).sum()
nocap_area = georeg_wcweights.loc['No_CAP','area']
gwdom_area = georeg_wcweights.loc['GW','area']
nocap_wt = nocap_area/(nocap_area+gwdom_area)
gwdom_wt = gwdom_area/(nocap_area+gwdom_area)

# %% Creating colors
# Matching map
cap = '#C6652B'
# noCAP = '#EDE461' # This is one from the map but it's too bright and hard to see
noCAP = '#CCC339' # This color but darker for lines
GWdom = '#3B76AF'
mixed = '#6EB2E4'
swdom = '#469B76'
specialyears = 'darkgray'

drought_color = '#ffa6b8'
wet_color = '#b8d3f2'

# %% Figure 7a
# For Depth to Water by regulation
ds = cat_wl2_reg
min_yr = 1975
mx_yr = 2020
betterlabels = ['Regulated','Unregulated'] 

f = ds[(ds.index >= min_yr) & (ds.index <= mx_yr)]
columns = ds.columns
column_list = ds.columns.tolist()

stats = pd.DataFrame()
# for i in range(1, 12, 1):
for i in column_list:
        df = f[i]
        #print(df)
        y=np.array(df.values, dtype=float)
        x=np.array(pd.to_datetime(df).index.values, dtype=float)
        slope, intercept, r_value, p_value, std_err =sp.linregress(x,y)
        stats = stats._append({'slope': slope, 
                              'int':intercept, 
                              'rsq':r_value*r_value, 
                              'p_val':p_value, 
                              'std_err':std_err, 
                              'mean': np.mean(y),
                              'var': np.var(y),
                              'sum': np.sum(y)
                              },
                              ignore_index=True)


stats.index = betterlabels
stats1 = stats.transpose()
print(stats1)

# -- Data visualization --
xf = np.linspace(min(x),max(x),100)
xf1 = xf.copy()
#xf1 = pd.to_datetime(xf1)
m1 = round(stats1.loc['slope','Regulated'], 2)
m2 = round(stats1.loc['slope','Unregulated'], 2)
yint1 = round(stats1.loc['int','Regulated'], 2)
yint2 = round(stats1.loc['int','Unregulated'], 2)
pval1 = round(stats1.loc['p_val', 'Regulated'], 4)
pval2 = round(stats1.loc['p_val', 'Unregulated'], 4)

yf1 = (m1*xf)+yint1
yf2 = (m2*xf)+yint2

fig, ax = plt.subplots(1, 1, figsize = (8,5))

min_y = 75
max_y = 300
fsize = 12

# Drought Year Shading
a = 1988.5
b = 1990.5
c = 1995.5
d = 1996.5
# e = 1999.5
# f = 2000.5
g = 2001.5
h = 2003.5
i = 2005.5
j = 2007.5
k = 2011.5
l = 2014.5
m = 2017.5
n = 2018.5
plt.axvspan(a, b, color=drought_color, alpha=0.5, lw=0
            , label="Severe Drought"
            )
plt.axvspan(c, d, color=drought_color, alpha=0.5, lw=0)
# plt.axvspan(e, f, color=drought_color, alpha=0.5, lw=0)
plt.axvspan(g, h, color=drought_color, alpha=0.5, lw=0)
plt.axvspan(i, j, color=drought_color, alpha=0.5, lw=0)
plt.axvspan(k, l, color=drought_color, alpha=0.5, lw=0)
plt.axvspan(m, n, color=drought_color, alpha=0.5, lw=0)

ax.plot(ds['R'], label='Regulated', color=cap) 
ax.plot(ds['U'], label='Unregulated', color=GWdom) 

ax.plot(xf1, yf1,"-.",color='k',label='Linear Trendline', lw=1)
ax.plot(xf1, yf1,"-.",color=cap, lw=1)
ax.plot(xf1, yf2,"-.",color=GWdom, lw=1)

ax.set_xlim(min_yr,mx_yr)
ax.set_ylim(min_y,max_y)
# ax.grid(True)
ax.grid(visible=True,which='major')
ax.grid(which='minor',color='#EEEEEE', lw=0.8)
ax.set_xlabel('Year', fontsize=fsize)
ax.set_ylabel('Depth to Water (ft)',fontsize=fsize)
ax.minorticks_on()
fig.set_dpi(600.0)
# ax.set_title('a)',loc='left',pad=15)
ax.legend(loc='upper left')

#Putting Grace on a secondary axis
ax2 = ax.twinx()
ax2.plot(grace_yearly['0'], label='State Average LWE', color='k',zorder=1)
ax2.set_ylim([15, -15])
ax2.set_ylabel(u'Δ LWE (cm)',fontsize=fsize)
ax2.legend(loc='lower right')

plt.savefig(figurepath+'Figure7a', bbox_inches = 'tight')

# %% Figure 7c
# For Depth to Water by SW Access
ds = cat_wl2_SW
min_yr = 1975
mx_yr = 2020
betterlabels = ['Recieves CAP'
                # ,'GW Dominated (Regulated)'
                ,'Surface Water Dominated'
                ,'GW Dominated'
                ,'Mixed Source'] 

f = ds[(ds.index >= min_yr) & (ds.index <= mx_yr)]
columns = ds.columns
column_list = ds.columns.tolist()

stats = pd.DataFrame()
for i in column_list:
        df = f[i]
        # df = f[i].pct_change()
        #print(df)
        y=np.array(df.values, dtype=float)
        x=np.array(pd.to_datetime(df).index.values, dtype=float)
        slope, intercept, r_value, p_value, std_err =sp.linregress(x,y)
        stats = stats._append({'slope': slope, 'int':intercept, 
                              'rsq':r_value*r_value, 'p_val':p_value, 
                              'std_err':std_err, 'mean': np.mean(y),
                              'var': np.var(y),'sum': np.sum(y)},
                              ignore_index=True)

stats.index = betterlabels
stats1 = stats.transpose()
print(stats1)
# -- Data visualization --
xf = np.linspace(min(x),max(x),100)
xf1 = xf.copy()
m1 = round(stats1.loc['slope',betterlabels[0]], 2)
m2 = round(stats1.loc['slope',betterlabels[3]], 2)
m3 = round(stats1.loc['slope',betterlabels[4]], 2)
m4 = round(stats1.loc['slope',betterlabels[1]], 2)
m5 = round(stats1.loc['slope',betterlabels[2]], 2)
yint1 = round(stats1.loc['int',betterlabels[0]], 2)
yint2 = round(stats1.loc['int',betterlabels[3]], 2)
yint3 = round(stats1.loc['int',betterlabels[4]], 2)
yint4 = round(stats1.loc['int',betterlabels[1]], 2)
yint5 = round(stats1.loc['int',betterlabels[2]], 2)
rsq1 = round(stats1.loc['rsq', betterlabels[0]], 4)
rsq2 = round(stats1.loc['rsq', betterlabels[3]], 4)
rsq3 = round(stats1.loc['rsq', betterlabels[4]], 4)
rsq4 = round(stats1.loc['rsq', betterlabels[1]], 4)
rsq5 = round(stats1.loc['rsq', betterlabels[2]], 4)
pval1 = round(stats1.loc['p_val', betterlabels[0]], 4)
pval2 = round(stats1.loc['p_val', betterlabels[3]], 4)
pval3 = round(stats1.loc['p_val', betterlabels[4]], 4)
pval4 = round(stats1.loc['p_val', betterlabels[1]], 4)
pval5 = round(stats1.loc['p_val', betterlabels[2]], 4)
yf1 = (m1*xf)+yint1
yf2 = (m2*xf)+yint2
yf3 = (m3*xf)+yint3
yf4 = (m4*xf)+yint4
yf5 = (m5*xf)+yint5

fig, ax = plt.subplots(1, 1, figsize = (8,5))

# Drought Year Shading
a = 1988.5
b = 1990.5
c = 1995.5
d = 1996.5
# e = 1999.5
# f = 2000.5
g = 2001.5
h = 2003.5
i = 2005.5
j = 2007.5
k = 2011.5
l = 2014.5
m = 2017.5
n = 2018.5
plt.axvspan(a, b, color=drought_color, alpha=0.5, lw=0, label="Drought")
plt.axvspan(c, d, color=drought_color, alpha=0.5, lw=0)
# plt.axvspan(e, f, color=drought_color, alpha=0.5, lw=0)
plt.axvspan(g, h, color=drought_color, alpha=0.5, lw=0)
plt.axvspan(i, j, color=drought_color, alpha=0.5, lw=0)
plt.axvspan(k, l, color=drought_color, alpha=0.5, lw=0)
plt.axvspan(m, n, color=drought_color, alpha=0.5, lw=0)


ax.plot(xf1, yf1,"-.",color=cap, lw=1)
ax.plot(xf1, yf2,"-.",color=GWdom, lw=1)
ax.plot(xf1, yf3,"-.",color=mixed, lw=1)
ax.plot(xf1, yf4,"-.",color='#CCC339', lw=1)
ax.plot(xf1, yf5,"-.",color=swdom, lw=1)

min_y = 75
max_y = 300
fsize = 12

ax.plot(ds['CAP'], label=betterlabels[0], color=cap,zorder=2)
ax.plot(ds['No_CAP'], label=betterlabels[1], color='#CCC339',zorder=2) 
ax.plot(ds['SW'], label=betterlabels[2], color=swdom,zorder=2) 
ax.plot(ds['Mix'], label=betterlabels[4], color=mixed,zorder=2)
ax.plot(ds['GW'], label=betterlabels[3], color=GWdom,zorder=2)  

ax.set_xlim([min_yr,mx_yr])
ax.set_ylim(min_y,max_y)
# ax.grid(True)
ax.grid(visible=True,which='major')
ax.grid(which='minor',color='#EEEEEE', lw=0.8)
ax.set_xlabel('Year', fontsize=fsize)
ax.set_ylabel('Depth to Water (ft)',fontsize=fsize)
ax.minorticks_on()
fig.set_dpi(600.0)
# ax.set_title('c)',fontsize = fsize,loc='left',pad=15)
ax.legend()

#Putting Grace on a secondary axis
ax2 = ax.twinx()
ax2.plot(grace_yearly['0'], label='State Average LWE', color='k',zorder=1)
ax2.set_ylim([15, -15])
ax2.set_ylabel(u'Δ LWE (cm)',fontsize=fsize)
ax2.legend(loc='lower right')

plt.savefig(figurepath+'Figure7c', bbox_inches = 'tight')

# %%
# Graph state average DTW and GRACE

adb_statemean = annual_db.mean()
adb_meandf = pd.DataFrame(adb_statemean)
# adb_meandf = adb_meandf.index.astype(int)
f = adb_meandf
f.reset_index(inplace=True)
f['index'] = pd.to_numeric(f['index'])
f['index'] = f['index'].astype(int)
f.set_index('index', inplace=True)

adb_meandf = f


ds = adb_meandf
min_yr = 1975
mx_yr = 2020
# betterlabels = ['Recieves CAP (Regulated)'
#                 ,'GW Dominated (Regulated)'
#                 ,'Surface Water Dominated'
#                 ,'GW Dominated'
#                 ,'Mixed Source']

betterlabels = ['State Average DTW'] 

f = ds[(ds.index >= min_yr) & (ds.index <= mx_yr)]
columns = ds.columns
column_list = ds.columns.tolist()

stats = pd.DataFrame()
for i in column_list:
        df = f[i]
        # df = f[i].pct_change()
        #print(df)
        y=np.array(df.values, dtype=float)
        x=np.array(pd.to_datetime(df).index.values, dtype=float)
        slope, intercept, r_value, p_value, std_err =sp.linregress(x,y)
        stats = stats._append({'slope': slope, 'int':intercept, 
                              'rsq':r_value*r_value, 'p_val':p_value, 
                              'std_err':std_err, 'mean': np.mean(y),
                              'var': np.var(y),'sum': np.sum(y)},
                              ignore_index=True)

stats.index = betterlabels
stats1 = stats.transpose()
print(stats1)
# # -- Data visualization --
xf = np.linspace(min(x),max(x),100)
xf1 = xf.copy()
m1 = round(stats1.loc['slope',betterlabels[0]], 2)
yint1 = round(stats1.loc['int',betterlabels[0]], 2)
rsq1 = round(stats1.loc['rsq', betterlabels[0]], 4)
pval1 = round(stats1.loc['p_val', betterlabels[0]], 4)
yf1 = (m1*xf)+yint1

fig, ax = plt.subplots(1, 1, figsize = (9,5))

ax.plot(xf1, yf1,"-.",color=GWdom, lw=1)

min_y = 75
max_y = 300
fsize = 12

# Drought Year Shading
a = 1988.5
b = 1990.5
c = 1995.5
d = 1996.5
# e = 1999.5
# f = 2000.5
g = 2001.5
h = 2003.5
i = 2005.5
j = 2007.5
k = 2011.5
l = 2014.5
m = 2017.5
n = 2018.5
plt.axvspan(a, b, color=drought_color, alpha=0.5, lw=0
            # , label="Drought"
            )
plt.axvspan(c, d, color=drought_color, alpha=0.5, lw=0)
# plt.axvspan(e, f, color=drought_color, alpha=0.5, lw=0)
plt.axvspan(g, h, color=drought_color, alpha=0.5, lw=0)
plt.axvspan(i, j, color=drought_color, alpha=0.5, lw=0)
plt.axvspan(k, l, color=drought_color, alpha=0.5, lw=0)
plt.axvspan(m, n, color=drought_color, alpha=0.5, lw=0)

ax.plot(ds, label='Arizona Average DTW', color=GWdom,zorder=2)  

ax.set_xlim([min_yr,mx_yr])
ax.set_ylim(min_y,max_y)
# ax.grid(True)
ax.grid(visible=True,which='major')
ax.grid(which='minor',color='#EEEEEE', lw=0.8)
ax.set_xlabel('Year', fontsize=fsize)
ax.set_ylabel('Depth to Water (ft)',fontsize=fsize)
ax.minorticks_on()
fig.set_dpi(600.0)
# ax.set_title('c)',fontsize = fsize,loc='left',pad=15)
ax.legend()

#Putting Grace on a secondary axis
ax2 = ax.twinx()
ax2.plot(grace_yearly['0'], label='State Average LWE', color='k',zorder=1)
ax2.set_ylim([15, -15])
ax2.set_ylabel(u'Δ LWE (cm)',fontsize=fsize)
ax2.legend(loc='lower right')

plt.savefig(figurepath+'BonusFigure_StateAverageDTWandGrace', bbox_inches = 'tight')
# %%
