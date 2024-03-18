# ------ Drought Indice Analysis --------
# Written by Danielle Tadych, May 2022
#    modified 9/26/2023

# The purpose of this code is to pick out drought periods.

# The dataset needed is nClimDiv text file from the Global Historical Climatology Network (GHCN).
# This data in particular is averaged for the state of Arizona.


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
import datetime as dt
from matplotlib.transforms import Bbox
import seaborn as sns
import numpy as np
import pandas as pd
from shapely.geometry import box
import geopandas as gp
#import earthpy as et
import scipy.stats as sp

# Assign Data paths
datapath = '../Data/Input_files/'
outputpath = '../Data/Output_files/'
shapepath = '../Data/Shapefiles/'

# %% Creating colors
c_1 = '#8d5a99' # Reservation
c_2 = "#d7191c" # Regulated with CAP (Water Category Color)
c_3 = '#e77a47' # Regulated without CAP (Water Category Color)
c_4 = '#2cbe21' # Lower CO River - SW (Water Category Color)
c_5 = '#2f8c73' # Upper CO River - Mixed (Water Category Color)
c_6 = '#6db7e8' # SE - GW
c_7 = '#165782' # NW - GW (Water Category color)
c_8 = '#229ce8' # SC - GW
c_9 = '#1f78b4' # NE - GW
c_10 = '#41bf9e' # N - Mixed
c_11 = '#7adec4' # C - Mixed
drought_color = '#ffa6b8'
wet_color = '#b8d3f2'

# %% Read in the file
filename = 'nClimDiv_AZ_GHCN.txt'
filepath = os.path.join(datapath, filename)
nclimdata = pd.read_csv(filepath)
nclimdata

# %%
nclimdata['date'] = pd.to_datetime(nclimdata['YearMonth'], format='%Y%m', errors='coerce').dropna()
nclimdata

# %%
nclimdata = nclimdata.rename(columns = {'   PDSI':'PDSI', '   PHDI':'PHDI'})
# %%
nclimdata = nclimdata.rename(columns={'    PCP':'Precip'})
precip = nclimdata[['date','Precip']]
precip['year'] = precip['date'].dt.year
yearly_precip = pd.pivot_table(precip, index=["year"], values=['Precip'], dropna=False, aggfunc=np.sum)
yearly_precip.plot()

# %%
ds = yearly_precip
minyear=1975
maxyear=2020
name = "Average Precipitation for AZ from " + str(minyear) + " to " + str(maxyear)
min_y = -6
max_y = 6
fsize = 12

fig, ax = plt.subplots(figsize = (9,5))
ax.plot(ds['Precip']
        # , label='Precip'
        , color='#e77a47'
        , lw=2
        ) 

# ax.plot(ds['wet'],label='wet',color='black',zorder = 5)
# ax.plot(ds['dry'],label='Cutoff Value',color='black', zorder=5)
# a = 1975
# b = 1977.5
# c = 1980.5
# d = 1981.5
# e = 1988.5
# f = 1990.5
# g = 1995.5
# h = 1997.5
# i = 1998.5
# j = 2004.5
# k = 2005.5
# l = 2009.5
# m = 2010.5
# n = 2018.5
# plt.axvspan(a, b, color=drought_color, alpha=0.5, lw=0, label="Drought")
# plt.axvspan(c, d, color=drought_color, alpha=0.5, lw=0)
# plt.axvspan(e, f, color=drought_color, alpha=0.5, lw=0)
# plt.axvspan(g, h, color=drought_color, alpha=0.5, lw=0)
# plt.axvspan(i, j, color=drought_color, alpha=0.5, lw=0)
# plt.axvspan(k, l, color=drought_color, alpha=0.5, lw=0)
# plt.axvspan(m, n, color=drought_color, alpha=0.5, lw=0)

a = 1988.5
b = 1990.5
c = 1995.5
d = 1996.5
e = 2001.5
f = 2003.5
g = 2005.5
h = 2007.5
i = 2011.5
j = 2014.5
k = 2017.5
l= 2018.5
plt.axvspan(a, b, color=drought_color, alpha=0.5, lw=0, label="Drought")
plt.axvspan(c, d, color=drought_color, alpha=0.5, lw=0)
plt.axvspan(e, f, color=drought_color, alpha=0.5, lw=0)
plt.axvspan(g, h, color=drought_color, alpha=0.5, lw=0)
plt.axvspan(i, j, color=drought_color, alpha=0.5, lw=0)
plt.axvspan(k, l, color=drought_color, alpha=0.5, lw=0)

ax.set_xlim(minyear,maxyear)
# ax.set_ylim(min_y,max_y)
ax.minorticks_on()
ax.grid(visible=True,which='major')
ax.grid(which='minor',color='#EEEEEE', lw=0.8)
ax.set_title(name, fontsize=14)
ax.set_xlabel('Year', fontsize=fsize)
ax.set_ylabel('Index Values',fontsize=fsize)
ax.legend(loc = [1.04, 0.40], fontsize = fsize)
fig.set_dpi(600.0)
# %%
pdsi = nclimdata[['date','PDSI','PHDI']]
pdsi
# %%
pdsi.describe()
# %%
pdsi = pdsi.set_index('date')
pdsi
# %%
pdsi.plot()

# %%
pdsi = pdsi.reset_index()
# %%
pdsi['In_year'] = pdsi['date'].dt.year
pdsi


# %%
yearly_pdsi = pd.pivot_table(pdsi, index=["In_year"], values=["PDSI", 'PHDI'], dropna=False, aggfunc=np.mean)
yearly_pdsi
# %%
yearly_pdsi.plot()

#%%
drought_color = '#ffa6b8'
wet_color = '#b8d3f2'

# %%
value = 3
yearly_pdsi['wet'] = value
yearly_pdsi['dry'] = -value
yearly_pdsi
#  PDSI
ds = yearly_pdsi
minyear=1975
maxyear=2020
# name = "Average PDSI and PHDI for AZ from " + str(minyear) + " to " + str(maxyear)
name = "Average PDSI for AZ from " + str(minyear) + " to " + str(maxyear)
min_y = -6
max_y = 6
fsize = 12

fig, ax = plt.subplots(figsize = (9,5))
#ax.plot(ds[1.0], label='Reservation', color=c_1)
# ax.plot(ds['CAP'], label='CAP', color=c_2)
ax.plot(ds['PHDI'], label='PHDI'
        # , color='blue'
        , lw=3
        ) 
ax.plot(ds['PDSI']
        # ,'-.'
        , label='PDSI'
        # , color='default'
        , lw=2
        ) 

# ax.plot(ds['wet'],label='wet',color='black',zorder = 5)
ax.plot(ds['dry'],'-.',label='Cutoff Value',color='black', zorder=5)
# a = 1975
# b = 1977.5
# c = 1980.5
# d = 1981.5
# e = 1988.5
# f = 1990.5
# g = 1995.5
# h = 1997.5
# i = 1998.5
# j = 2004.5
# k = 2005.5
# l = 2009.5
# m = 2010.5
# n = 2018.5
# plt.axvspan(a, b, color=drought_color, alpha=0.5, lw=0, label="Drought")
# plt.axvspan(c, d, color=drought_color, alpha=0.5, lw=0)
# plt.axvspan(e, f, color=drought_color, alpha=0.5, lw=0)
# plt.axvspan(g, h, color=drought_color, alpha=0.5, lw=0)
# plt.axvspan(i, j, color=drought_color, alpha=0.5, lw=0)
# plt.axvspan(k, l, color=drought_color, alpha=0.5, lw=0)
# plt.axvspan(m, n, color=drought_color, alpha=0.5, lw=0)

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
plt.axvspan(a, b, color=drought_color, alpha=0.5, lw=0, label="Severe Drought")
plt.axvspan(c, d, color=drought_color, alpha=0.5, lw=0)
# plt.axvspan(e, f, color=drought_color, alpha=0.5, lw=0)
plt.axvspan(g, h, color=drought_color, alpha=0.5, lw=0)
plt.axvspan(i, j, color=drought_color, alpha=0.5, lw=0)
plt.axvspan(k, l, color=drought_color, alpha=0.5, lw=0)
plt.axvspan(m, n, color=drought_color, alpha=0.5, lw=0)

ax.set_xlim(minyear,maxyear)
ax.set_ylim(min_y,max_y)
ax.minorticks_on()
ax.grid(visible=True,which='major')
ax.grid(which='minor',color='#EEEEEE', lw=0.8)
# ax.set_title(name, fontsize=14)
ax.set_xlabel('Year', fontsize=fsize)
ax.set_ylabel('Index Values',fontsize=fsize)
ax.legend(loc = [1.04, 0.40], fontsize = fsize)
fig.set_dpi(600.0)
# plt.savefig(outputpath+name+'cutoffval_'+str(value), bbox_inches = 'tight')


# %%
value = 3
print("Drought is considered > -", value)
yearly_pdsi['wet'] = value
yearly_pdsi['dry'] = -value
yearly_pdsi

drought = yearly_pdsi[yearly_pdsi['PHDI']<=-value]
wet = yearly_pdsi[yearly_pdsi['PHDI']>=value]
drought = drought[drought.index >= 1975]
wet = wet[wet.index >= 1975]

drought

print()
print("Drought Year Info:")
print(drought['PHDI'].describe())
print()
print("Wet Year Info:")
print(wet['PHDI'].describe())

#  PHDI
ds = yearly_pdsi
minyear=1975
maxyear=2020
name = "Average PHDI for AZ from " + str(minyear) + " to " + str(maxyear)
min_y = -6
max_y = 6
fsize = 12

fig, ax = plt.subplots(figsize = (9,5))
#ax.plot(ds[1.0], label='Reservation', color=c_1)
# ax.plot(ds['CAP'], label='CAP', color=c_2)
# ax.plot(ds['PDSI'], '-o'
#         , label='PDSI'
#         , color='grey'
#         , lw=1
#         ) 
ax.plot(ds['PHDI'], '-',label='PHDI'
        , color='grey'
        , lw=1
        ) 
ax.plot(ds['wet'],color='black',zorder = 5)
ax.plot(ds['dry'],color='black', zorder=5)
ax.plot(drought['PHDI'],'o',label='dry',color='red', zorder=5)
ax.plot(wet['PHDI'],'o',label='wet',color='blue', zorder=5)

ax.set_xlim(minyear,maxyear)
ax.set_ylim(min_y,max_y)
ax.minorticks_on()
ax.grid(visible=True,which='major')
ax.grid(which='minor',color='#EEEEEE', lw=0.8)
# ax.set_title(name, fontsize=14)
ax.set_xlabel('Year', fontsize=fsize)
ax.set_ylabel('Index Values',fontsize=fsize)
# ax.legend(loc = [1.04, 0.40], fontsize = fsize)
fig.set_dpi(600.0)
# plt.savefig(outputpath+name+'cutoffval_'+str(value))

#%%
# value = 1.5
print("Drought is considered > -", value)
yearly_pdsi['wet'] = value
yearly_pdsi['dry'] = -value
yearly_pdsi

drought = yearly_pdsi[yearly_pdsi['PDSI']<=-value]
wet = yearly_pdsi[yearly_pdsi['PDSI']>=value]
drought = drought[drought.index >= 1975]
wet = wet[wet.index >= 1975]

print()
print("Drought Year Info:")
print(drought['PDSI'].describe())
print()
print("Wet Year Info:")
print(wet['PDSI'].describe())

#  PDSI
ds = yearly_pdsi
minyear=1975
maxyear=2020
name = "Average PDSI for AZ from " + str(minyear) + " to " + str(maxyear)
min_y = -6
max_y = 6
fsize = 12

fig, ax = plt.subplots(figsize = (9,5))
#ax.plot(ds[1.0], label='Reservation', color=c_1)
# ax.plot(ds['CAP'], label='CAP', color=c_2)
# ax.plot(ds['PDSI'], '-o'
#         , label='PDSI'
#         , color='grey'
#         , lw=1
#         ) 
ax.plot(ds['PDSI'], '-',label='PDSI'
        , color='grey'
        , lw=1
        ) 
ax.plot(ds['wet'],color='black',zorder = 5)
ax.plot(ds['dry'],color='black', zorder=5)
ax.plot(drought['PDSI'],'o',label='dry',color='red', zorder=5)
ax.plot(wet['PDSI'],'o',label='wet',color='blue', zorder=5)

a = 1988.5
b = 1989.5
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

ax.set_xlim(minyear,maxyear)
ax.set_ylim(min_y,max_y)
ax.minorticks_on()
ax.grid(visible=True,which='major')
ax.grid(which='minor',color='#EEEEEE', lw=0.8)
ax.set_title(name, fontsize=14)
ax.set_xlabel('Year', fontsize=fsize)
ax.set_ylabel('Index Values',fontsize=fsize)
# ax.legend(loc = [1.04, 0.40], fontsize = fsize)
fig.set_dpi(600.0)

# %%
analysis_period = yearly_pdsi[yearly_pdsi.index>=1975]
# del analysis_period['wet']
# del analysis_period['dry']

analysis_period.to_csv('../Data/Input_files/Yearly_DroughtIndices.csv')

# %%
analysis_period
# %%
drought_phdi = analysis_period[analysis_period['PHDI']<=-1]
drought_pdsi = analysis_period[analysis_period['PDSI']<=-1]
print(drought_pdsi)
# %%
print(drought_phdi)
# %%
df = analysis_period
value = 1
df['PHDI_status'] = 'Normal'
df.loc[df['PHDI'] <= -value, 'PHDI_status'] = 'Drought' 
df.loc[df['PHDI'] >= value, 'PHDI_status'] = 'Wet' 
df['PDSI_status'] = 'Normal'
df.loc[df['PDSI'] <= -value, 'PDSI_status'] = 'Drought' 
df.loc[df['PDSI'] >= value, 'PDSI_status'] = 'Wet' 

df

analysis_period.to_csv('../Data/Input_files/YearlyDrought_'+str(value)+'.csv')

# %%
analysis_period[['PHDI_status','PDSI_status']].describe()
# %%
print('Drought info: ', analysis_period[analysis_period['PDSI_status']=='Drought'].describe())

# %% Reading in new PDSI Data
filename = 'nClimDiv_AZ_GHCN.txt'
filepath = os.path.join(datapath, filename)
nclimdata = pd.read_csv(filepath)
nclimdata

# %% From ChatGPT
# Replace 'your_folder_path' with the path to your folder containing CSV files
folder_path = datapath+'/PDSI'

# Get a list of all files in the folder
files = os.listdir(folder_path)

# %%
# Filter out only CSV files (assuming they have a .csv extension)
csv_files = [file for file in files if file.endswith('.csv')]

# Create an empty list to store DataFrames
dataframes = []

# Loop through each CSV file and read it into a pandas DataFrame
for csv_file in csv_files:
    # Construct the full path to the CSV file
    csv_file_path = os.path.join(folder_path, csv_file)

    # Read the CSV file into a pandas DataFrame
    df = pd.read_csv(csv_file_path)

        # Convert the first column to datetime format if it's in 'yyyymmdd' format
    df.index = pd.to_datetime(df.index, format='%Y%m%d')

    df['Year'] = df.index.year
    print(df.columns[0])
#     df = df.rename(columns={df.columns[0]: 'PDSI'})
    # Append the DataFrame to the list
    dataframes.append(df)
# See the first one
dataframes[0]
# %%
# Group by the date column and calculate the average for each group
combined_df = pd.concat(dataframes
                        # , ignore_index=True
                        )
combined_df

# %%
combined_df.plot()

# %%
# average_df = pd.pivot_table(combined_df, index=["Year"], values=['PDSI'], dropna=False, aggfunc=np.mean)

average_df = combined_df.groupby('Year').mean()
average_df
# %%
average_df.plot()
# %%
yearly_pdsi

# %%
summarystats = average_df.transpose()
summarystats = summarystats.describe()
summarystats = summarystats.transpose()
summarystats
# %%
value = 3
yearly_pdsi['wet'] = value
yearly_pdsi['dry'] = -value
yearly_pdsi
#  PDSI
ds = yearly_pdsi
minyear=2000
maxyear=2022
# name = "Average PDSI and PHDI for AZ from " + str(minyear) + " to " + str(maxyear)
name = "Comparing prior PDSI Values to New Values"
min_y = -6
max_y = 6
fsize = 12

fig, ax = plt.subplots(figsize = (9,5))

# a = 1988.5
# b = 1990.5
# c = 1995.5
# d = 1996.5
# e = 2020.5
# f = 2021.5
# g = 2001.5
# h = 2003.5
# i = 2005.5
# j = 2007.5
# k = 2011.5
# l = 2014.5
# m = 2017.5
# n = 2018.5
# plt.axvspan(a, b, color=drought_color, alpha=0.5, lw=0, label="Severe Drought")
# plt.axvspan(c, d, color=drought_color, alpha=0.5, lw=0)
# plt.axvspan(e, f, color=drought_color, alpha=0.5, lw=0)
# plt.axvspan(g, h, color=drought_color, alpha=0.5, lw=0)
# plt.axvspan(i, j, color=drought_color, alpha=0.5, lw=0)
# plt.axvspan(k, l, color=drought_color, alpha=0.5, lw=0)
# plt.axvspan(m, n, color=drought_color, alpha=0.5, lw=0)

#ax.plot(ds[1.0], label='Reservation', color=c_1)
# ax.plot(ds['CAP'], label='CAP', color=c_2)
# ax.plot(ds['PHDI'], label='PHDI'
#         # , color='blue'
#         , lw=3
#         ) 
 
# ax.plot(average_df['PDSI']
#         # ,'-.'
#         # , label='New PDSI'
#         # , color='default'
#         , lw=2
#         )
ax.plot(summarystats['mean']
        ,label='New PDSI Mean')
ax.fill_between(summarystats.index
                ,summarystats['max']
                ,summarystats['min']
                ,label='New PDSI Range'
                ,alpha=0.3)
ax.plot(ds['PDSI']
        # ,'-.'
        , label='Old PDSI'
        # , color='default'
        , lw=2
        )

# ax.plot(ds['wet'],label='wet',color='black',zorder = 5)
# ax.plot(ds['dry'],'-.',label='Cutoff Value',color='black', zorder=5)

ax.set_xlim(minyear,maxyear)
ax.set_ylim(min_y,max_y)
ax.minorticks_on()
ax.grid(visible=True,which='major')
ax.grid(which='minor',color='#EEEEEE', lw=0.8)
# ax.set_title(name, fontsize=14)
ax.set_xlabel('Year', fontsize=fsize)
ax.set_ylabel('Index Values',fontsize=fsize)
ax.legend(loc = [1.04, 0.40], fontsize = fsize)
fig.set_dpi(600.0)
# plt.savefig(outputpath+name+'cutoffval_'+str(value), bbox_inches = 'tight')

# %%
average_df.to_csv('../Data/Input_files/nclimdiv_PDSI_Azdivs.csv')
summarystats.to_csv('../Data/Input_files/NewPDSI_manualavg_11272023.csv')
# %%  --- Actually updated files
# filename = 'NOAA_PDSI_Timesseries_11282023.csv'
filename = 'NOAA_PDSI_Timesseries_12032023.csv'
filepath = os.path.join(datapath, filename)
pdsi = pd.read_csv(filepath, header=3)

pdsi['Date'] = pd.to_datetime(pdsi['Date'], format='%Y%m', errors='coerce').dropna()
pdsi['In_year'] = pdsi['Date'].dt.year
pdsi
# %%
filename = 'NOAA_PHDI_Timesseries_12032023.csv'
filepath = os.path.join(datapath, filename)
phdi = pd.read_csv(filepath, header=3)

phdi['Date'] = pd.to_datetime(phdi['Date'], format='%Y%m', errors='coerce').dropna()
phdi['In_year'] = phdi['Date'].dt.year
phdi

# %%
yearly_pdsi_new = pd.pivot_table(pdsi, index=["In_year"], values=["Value"], dropna=False, aggfunc=np.mean)
yearly_pdsi_new
yearly_pdsi_new = yearly_pdsi_new.rename(columns = {'Value':'PDSI'})

# %%
yearly_phdi_new = pd.pivot_table(phdi, index=["In_year"], values=["Value"], dropna=False, aggfunc=np.mean)
yearly_phdi_new
yearly_phdi_new = yearly_phdi_new.rename(columns = {'Value':'PHDI'})

# %%
value = 3
yearly_pdsi['wet'] = value
yearly_pdsi['dry'] = -value
yearly_pdsi
#  PDSI
ds = yearly_pdsi
minyear=2000
maxyear=2023
# name = "Average PDSI and PHDI for AZ from " + str(minyear) + " to " + str(maxyear)
name = "Comparing prior PDSI Values to New Values"
min_y = -6
max_y = 6
fsize = 12

fig, ax = plt.subplots(figsize = (9,5))

a = 1988.5
b = 1990.5
c = 1995.5
d = 1996.5
e = 2020.5
f = 2021.5
g = 2001.5
h = 2003.5
i = 2005.5
j = 2007.5
k = 2011.5
l = 2014.5
m = 2017.5
n = 2018.5
plt.axvspan(a, b, color=drought_color, alpha=0.5, lw=0, label="Severe Drought")
plt.axvspan(c, d, color=drought_color, alpha=0.5, lw=0)
plt.axvspan(e, f, color=drought_color, alpha=0.5, lw=0)
plt.axvspan(g, h, color=drought_color, alpha=0.5, lw=0)
plt.axvspan(i, j, color=drought_color, alpha=0.5, lw=0)
plt.axvspan(k, l, color=drought_color, alpha=0.5, lw=0)
plt.axvspan(m, n, color=drought_color, alpha=0.5, lw=0)

#ax.plot(ds[1.0], label='Reservation', color=c_1)
# ax.plot(ds['CAP'], label='CAP', color=c_2)
# ax.plot(ds['PHDI'], label='PHDI'
#         # , color='blue'
#         , lw=3
#         ) 
 
# ax.plot(average_df['PDSI']
#         # ,'-.'
#         # , label='New PDSI'
#         # , color='default'
#         , lw=2
#         )
ax.plot(yearly_pdsi_new['PDSI']
        ,label='New PDSI Mean'
        )
# ax.fill_between(summarystats.index
#                 ,summarystats['max']
#                 ,summarystats['min']
#                 ,label='New PDSI Range'
#                 ,alpha=0.3)
ax.plot(ds['PDSI']
         ,'-.'
        , label='Old PDSI'
#         # , color='default'
#         , lw=2
        )

# ax.plot(ds['wet'],label='wet',color='black',zorder = 5)
ax.plot(ds['dry'],'-.',label='Cutoff Value',color='black', zorder=5)

ax.set_xlim(minyear,maxyear)
ax.set_ylim(min_y,max_y)
ax.minorticks_on()
ax.grid(visible=True,which='major')
ax.grid(which='minor',color='#EEEEEE', lw=0.8)
# ax.set_title(name, fontsize=14)
ax.set_xlabel('Year', fontsize=fsize)
ax.set_ylabel('Index Values',fontsize=fsize)
ax.legend(loc = [1.04, 0.40], fontsize = fsize)
fig.set_dpi(600.0)
# plt.savefig(outputpath+name+'cutoffval_'+str(value), bbox_inches = 'tight')

# %%
yearly_new_indices = yearly_pdsi_new.copy()
# %%
yearly_new_indices['PHDI'] = yearly_phdi_new['PHDI']
yearly_new_indices

# %%
value = 3 # Severe drought cutoff value

ds = yearly_new_indices
minyear=1975
maxyear=2023
# name = "Average PDSI and PHDI for AZ from " + str(minyear) + " to " + str(maxyear)
name = "Average PDSI and PHDI for AZ from " + str(minyear) + " to " + str(maxyear)
min_y = -6
max_y = 6
fsize = 12

ds['wet'] = value
ds['dry'] = -value

fig, ax = plt.subplots(figsize = (9,5))
#ax.plot(ds[1.0], label='Reservation', color=c_1)
# ax.plot(ds['CAP'], label='CAP', color=c_2)
ax.plot(ds['PHDI'], label='PHDI'
        # , color='blue'
        , lw=3
        ) 
ax.plot(ds['PDSI']
        # ,'-.'
        , label='PDSI'
        # , color='default'
        , lw=2
        ) 

# Severe Drought Shading
a = 1988
b = 1990
c = 1995
d = 1996
e = 2021
f = 2022
g = 2002
h = 2004
i = 2005
j = 2008
k = 2012
l = 2015
m = 2018
n = 2019
plt.axvspan(a, b, color=drought_color, alpha=0.5, lw=0, label="Severe Drought")
plt.axvspan(c, d, color=drought_color, alpha=0.5, lw=0)
plt.axvspan(e, f, color=drought_color, alpha=0.5, lw=0)
plt.axvspan(g, h, color=drought_color, alpha=0.5, lw=0)
plt.axvspan(i, j, color=drought_color, alpha=0.5, lw=0)
plt.axvspan(k, l, color=drought_color, alpha=0.5, lw=0)
plt.axvspan(m, n, color=drought_color, alpha=0.5, lw=0)

# ax.plot(ds['wet'],label='wet',color='black',zorder = 5)
ax.plot(ds['dry'],'-.',label='Cutoff Value',color='black', zorder=5)

ax.set_xlim(minyear,maxyear)
ax.set_ylim(max_y,min_y)
ax.minorticks_on()
ax.grid(visible=True,which='major')
ax.grid(which='minor',color='#EEEEEE', lw=0.8)
# ax.set_title(name, fontsize=14)
ax.set_xlabel('Year', fontsize=fsize)
ax.set_ylabel('Index Values',fontsize=fsize)
ax.legend(loc = [1.04, 0.40], fontsize = fsize)
fig.set_dpi(600.0)
plt.savefig(outputpath+name+'cutoffval_'+str(value), bbox_inches = 'tight')

# %%
analysis_period = yearly_new_indices[yearly_new_indices.index>=1975]
# del analysis_period['wet']
# del analysis_period['dry']

analysis_period.to_csv('../Data/Input_files/Yearly_DroughtIndices_updated12032023.csv')
# %%
