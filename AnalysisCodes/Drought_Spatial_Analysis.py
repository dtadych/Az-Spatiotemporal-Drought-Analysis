# ~~~~~~~~~~~ Well and Grace Drought Analysis Code ~~~~~~~~~~~
# Written by Danielle Tadych

# The purpose of this script is to analyze our master water level database and GRACE data
# for various types of drought. 


# %%
import Custom_functions

# Data paths
datapath = '../Data/Input_files/'
outputpath = '../Data/Output_files/'
shapepath = '../Data/Shapefiles/'

# %%  ==== Reading in the data ====
filename_georeg = 'georeg_reproject_fixed.shp'
filepath = os.path.join(shapepath, filename_georeg)
georeg = gp.read_file(filepath)
georeg.plot(cmap='viridis')

georeg['GEOREGI_NU'] = georeg['GEOREGI_NU'].astype('int64')
georeg.info()

# Read in the annual time series database
filename_ts = 'Wells55_GWSI_WLTS_DB_annual.csv'
filepath = os.path.join(datapath, filename_ts)
print(filepath)
annual_db = pd.read_csv(filepath, header=1, index_col=0)

annual_db.index.astype('int64')
annual_db2 = annual_db.reset_index(inplace=True)
annual_db2 = annual_db.rename(columns = {'year':'Combo_ID'})
annual_db2.head()

# For regulation
filepath = datapath+'Waterlevels_Regulation.csv'
cat_wl2_reg = pd.read_csv(filepath, index_col=0)
cat_wl2_reg.head()

# For Access to SW
filepath = datapath+'Waterlevels_AccesstoSW.csv'
cat_wl2_SW = pd.read_csv(filepath, index_col=0)
cat_wl2_SW.head()

# For georegion number
filepath = datapath+'Waterlevels_georegions.csv'
cat_wl2_georeg = pd.read_csv(filepath, index_col=0)

# Read in the drought indices
drought_indices = pd.read_csv(datapath+'Yearly_DroughtIndices.csv')
drought_indices = drought_indices.set_index('In_year')
drought_indices

# %%
cat_wl2_georeg = cat_wl2_georeg.transpose()
cat_wl2_georeg
# %%
cat_wl2_georeg.reset_index(inplace=True)
cat_wl2_georeg['index'] = pd.to_numeric(cat_wl2_georeg['index'])
cat_wl2_georeg.set_index('index', inplace=True)
cat_wl2_georeg.info()

# %%
cat_wl2_georeg = cat_wl2_georeg.transpose()
cat_wl2_georeg

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

reg_colors = [c_2,c_7]
georeg_colors = [c_1,c_2,c_3,c_4,c_5,c_6,c_7,c_8,c_9,c_10,c_11]
SW_colors = [c_2,c_3,c_4,c_5,c_7]

bar_watercatc = [c_2,c_3,c_4,c_5,c_7]


# Color blind palette
# https://jacksonlab.agronomy.wisc.edu/2016/05/23/15-level-colorblind-friendly-palette/
blind =["#000000","#004949","#009292","#ff6db6","#ffb6db",
 "#490092","#006ddb","#b66dff","#6db6ff","#b6dbff",
 "#920000","#924900","#db6d00","#24ff24","#ffff6d"]

# Matching new map

cap = '#C6652B'
# noCAP = '#EDE461' # This is one from the map
noCAP = '#CCC339' # This color but darker for lines
GWdom = '#3B76AF'
mixed = '#6EB2E4'
swdom = '#469B76'

# %% ====== Specialized Drought Analysis ======
# Wanting to look at 1) Drawdown 2) Anomaly's 3) Recovery
#   Decided from the drought indices analysis that the cutoff value is -3 for severe droughts

# %% Drought dictionary
DROUGHT_YEARS = {1:[1989,1990]
        ,2:[1996]
        ,3:[2002,2003]
        ,4:[2006,2007]
        ,5:[2012,2013,2014]
        ,6:[2018]}

print(DROUGHT_YEARS)

#%% Pre-drought
PREDROUGHT_YEARS = {1:[1988]
        ,2:[1995]
        ,3:[2001]
        ,4:[2005]
        ,5:[2011]
        ,6:[2017]}

print(PREDROUGHT_YEARS)

#%% Print the average PDSI and PHDI values

ds = drought_indices.copy()
columns = ds.columns
column_list = ds.columns.tolist()

ds['Status'] = 'Normal-Wet'
# wlanalysis_period

for x,y in DROUGHT_YEARS.items():
        ds.loc[y, 'Status'] = 'Drought '+str(x)


pdsi_avg = ds.groupby(['Status']).mean()
pdsi_avg

#%%
def addlabels(x,y):
    for i in range(len(x)):
        plt.text(i,y[i],y[i])
# %% Grouped bar chart of PDSI/PHDI Values
name = 'Average PDSI and PHDI Values Per Drought'

yearlabels = ["1989-1990"
                ,'1996'
                ,'2002-2003'
                ,'2006-2007'
                ,'2012-2014'
                ,'2018'
                ,'Normal/Wet Years']

pdsi_avg.index = yearlabels
pdsi_avg = pdsi_avg.transpose()
# del ds['Normal/Wet Years']
pdsi_avg
#%%
group_colors = [blind[5],blind[6],blind[2]
                ,blind[12],blind[11],blind[10]
                ,blind[0] #black
                ]

horlabel = 'Index Value'
fsize = 14

plt.rcParams["figure.dpi"] = 600
pdsi_avg.plot(figsize = (10,7),
        kind='bar',
        stacked=False,
        # title=name,
        color = group_colors,
        zorder = 2,
        width = 0.85,
        fontsize = fsize
        )
plt.title(name, fontsize = (fsize+2))
# plt.ylim([0,400])
plt.ylabel(horlabel, fontsize = fsize)
plt.xticks(rotation=0, fontsize = fsize-2)
plt.grid(axis='y', linewidth=0.5, zorder=0)
plt.legend(loc=[1.01,0.3],fontsize = fsize)

# plt.savefig(outputpath+name+'_groupedchart', bbox_inches = 'tight')

# %% Figure out which water level database you want
cat_wl2 = cat_wl2_reg.copy() 
# cat_wl2 = cat_wl2_SW.copy()
# cat_wl2 = cat_wl2_georeg.copy()

# cat_wl2 = wdc1_reg.copy()
# cat_wl2 = wdc2_reg.copy()
# cat_wl2 = wdc3_reg.copy()
# cat_wl2 = wdc1_SW.copy()
# cat_wl2 = wdc2_SW.copy()
# cat_wl2 = wdc3_SW.copy()

# Water Analysis period
wlanalysis_period = cat_wl2[cat_wl2.index>=1975]
# wlanalysis_period["UGW"]=wlanalysis_period['GW']
# del wlanalysis_period['GW']
# wlanalysis_period

#%%
# Anomaly's
ds = wlanalysis_period
columns = ds.columns
column_list = ds.columns.tolist()

dtw_anomalys = pd.DataFrame()
for i in column_list:
        dtw_anomalys[i] = wlanalysis_period[i] - wlanalysis_period[i].mean()

dtw_anomalys.head()

# %% Drawdown
ds = wlanalysis_period.copy()
columns = ds.columns
column_list = ds.columns.tolist()

ds['Status'] = 'Normal-Wet'
# wlanalysis_period

for x,y in DROUGHT_YEARS.items():
        ds.loc[y, 'Status'] = 'Drought '+str(x)


drawd_max = ds.groupby(['Status']).max()
drawd_max
#%%
ds = wlanalysis_period.copy()
columns = ds.columns
column_list = ds.columns.tolist()

ds['Status'] = 'Normal-Wet'

for x,y in PREDROUGHT_YEARS.items():
        ds.loc[y, 'pre_d'] = 'Drought '+str(x)

predrought = ds.groupby(['pre_d']).mean()
predrought

# %% Drawdown
drawdown = drawd_max - predrought
drawdown

# %% Checking for normality
# ds = wlanalysis_period
ds = dtw_anomalys
columns = ds.columns
column_list = ds.columns.tolist()

for i in column_list:
 fig, ax = plt.subplots(1,1)
 ax.hist(wlanalysis_period[i], bins=30)
 ax.set_title(i)

# %% If running a shifted correlation analysis,
#    change this to however many # years; 0 is no lag
lag = 0

print('Kendall Correlation coefficient')
for i in column_list:
        # print(' '+i+':')
        print(' '+str(i)+':')
# To normalize the data 
        # df1 = ds[i].pct_change()
        # df2 = drought_indices.PDSI.pct_change()
        df1 = ds[i]
        df2 = drought_indices.PDSI.shift(lag)
        print('  tau = ',round(df1.corr(df2, method='kendall'),3))
        print('  pval = ',round(df1.corr(df2, method=kendall_pval),4))

# %%
print('Spearman Correlation coefficient')
for i in column_list:
        print(' '+str(i)+':')
        # df1 = ds[i].pct_change()
        # df2 = drought_indices.PDSI.pct_change()
        df1 = ds[i]
        df2 = drought_indices.PDSI.shift(lag)
        print('  rho = ',round(df1.corr(df2, method='spearman'),3))
        print('  pval = ',round(df1.corr(df2, method=spearmanr_pval),4))

# %%
print('Pearson Correlation coefficient')
for i in column_list:
        print(' '+str(i)+':')
        # df1 = ds[i].pct_change()
        # df2 = drought_indices.PDSI.pct_change()
        df1 = ds[i]
        df2 = drought_indices.PDSI.shift(lag)
        r = df1.corr(df2, method='pearson')
        print('  rsq = ',round(r*r,3))
        print('  pval = ',round(df1.corr(df2, method=pearsonr_pval),4))


# %% Scatterplot of correlation values
ds = dtw_anomalys
# name = 'Comparing PDSI with Depth to Water Anomalies by Access to SW'
name = 'Comparing PDSI with Depth to Water Anomalies by Regulation'
# del ds['Res']
columns = ds.columns
column_list = ds.columns.tolist()
# betterlabels = ['Receives CAP (Regulated)','GW Dominated (Regulated)','Surface Water Dominated','GW Dominated','Mixed Source'] 
betterlabels = ['Regulated','Unregulated'] 
colors=[cap, GWdom]
# colors=[cap,noCAP, swdom, mixed, GWdom]

fig, ax = plt.subplots(figsize = (7,5))
x = drought_indices['PDSI']
for i,j,k in zip(column_list
                # ,reg_colors
                # , SW_colors
                , colors
                , betterlabels
                ):
        y = ds[i]
        ax.scatter(x,y
                , label=k
                , color=j
                )
        # Trendline: 1=Linear, 2=polynomial
        z = np.polyfit(x, y, 1)
        p = np.poly1d(z)
        plt.plot(x, p(x),'-'
                , color=j
                # ,label=(k+' Trendline')
                )


ax.set_xlabel('PDSI')
ax.set_ylabel('Depth to Water Anomalies (ft)')
ax.set_title(name,loc='center')
# ax.set_ylim(0,400)
fig.set_dpi(600)
plt.legend(loc = [1.05, 0.40])

plt.savefig(outputpath+name, bbox_inches='tight') 

# %% Grouped bar chart of individual drought anomlies
# cat_wl2 = wdc1_reg.copy() # Deep
# cat_wl2 = wdc2_reg.copy() # Midrange
# cat_wl2 = wdc3_reg.copy() # Shallow
# cat_wl2 = wdc1_SW.copy()
# cat_wl2 = wdc2_SW.copy()
# cat_wl2 = wdc3_SW.copy()
cat_wl2 = cat_wl2_SW
# cat_wl2 = cat_wl2_reg

# name = 'Average DTW Anomalies by Drought Period and Groundwater Regulation'
# name = 'Average DTW Anomalies by Drought Period and Access to SW'

# name = 'Deep Wells'
# name = 'Midrange Wells'
# name = 'Shallow Wells'

betterlabels = ['CAP','Regulated \n Groundwater','Surface \n Water','Unregulated \n Groundwater','Mixed \n GW/SW'] 
# betterlabels = ['GW Regulated','GW Unregulated'] 

yearlabels = ["1989-1990",'1996','2002-2003','2006-2007','2012-2014','2018','Normal/Wet Years']

#%%
# Water Analysis period
wlanalysis_period = cat_wl2[cat_wl2.index>=1975]
# wlanalysis_period["UGW"]=wlanalysis_period['GW']
# del wlanalysis_period['GW']
# wlanalysis_period

# Anomaly's
ds = wlanalysis_period.copy()
columns = ds.columns
column_list = ds.columns.tolist()

dtw_anomalys = pd.DataFrame()
for i in column_list:
        dtw_anomalys[i] = wlanalysis_period[i] - wlanalysis_period[i].mean()

# %%
ds = dtw_anomalys.copy()
# ds = drought_indices

ds['Status'] = 'Normal-Wet'
# wlanalysis_period

for x,y in DROUGHT_YEARS.items():
        ds.loc[y, 'Status'] = 'Drought '+str(x)

ds

ds_indd = ds.groupby(['Status']).mean()
ds_indd.index = yearlabels
ds_indd = ds_indd.transpose()
ds_indd.index = betterlabels
ds_indd

#%%
# group_colors = ['lightsalmon','tomato','orangered','r','brown','indianred','steelblue']

group_colors = [
                blind[5],blind[6],blind[2]
                ,blind[12],blind[11],blind[10]
                ,blind[0] #black
                ]

horlabel = 'DTW Anomaly (ft)'
fsize = 14

ds_indd.plot(figsize = (10,7),
        kind='bar',
        stacked=False,
        # title=name,
        color = group_colors,
        zorder = 2,
        width = 0.85,
        fontsize = fsize
        )
plt.title(name, fontsize = (fsize+2))
# plt.ylim([0,400])
plt.ylabel(horlabel, fontsize = fsize)
plt.xticks(rotation=0, fontsize = fsize-2)
plt.grid(axis='y', linewidth=0.5, zorder=0)
plt.legend(loc=[1.01,0.3],fontsize = fsize)
# plt.figure(dpi=600)

plt.savefig(outputpath+name+'_anomalies_GWREG_groupedchart', bbox_inches = 'tight')
# plt.savefig(outputpath+name+'_anomalies_SWAccess_groupedchart', bbox_inches = 'tight')

#%% Drawdown quick analysis
# cat_wl2 = cat_wl2_reg.copy() 
cat_wl2 = cat_wl2_SW.copy()
# cat_wl2 = cat_wl2_georeg.copy()

# cat_wl2 = wdc1_reg.copy()
# cat_wl2 = wdc2_reg.copy()
# cat_wl2 = wdc3_reg.copy()
# cat_wl2 = wdc1_SW.copy()
# cat_wl2 = wdc2_SW.copy()
# cat_wl2 = wdc3_SW.copy()

betterlabels = ['CAP','Regulated \n Groundwater','Surface \n Water','Unregulated \n Groundwater','Mixed \n GW/SW'] 
# betterlabels = ['GW Regulated','GW Unregulated'] 

# ---
wlanalysis_period = cat_wl2[cat_wl2.index>=1975]

ds = wlanalysis_period.copy()
columns = ds.columns
column_list = ds.columns.tolist()
ds['Status'] = 'Normal-Wet'
for x,y in DROUGHT_YEARS.items():
        ds.loc[y, 'Status'] = 'Drought '+str(x)

for x,y in PREDROUGHT_YEARS.items():
        ds.loc[y, 'pre_d'] = 'Drought '+str(x)
# ds

drawd_max = ds.groupby(['Status']).max()
predrought = ds.groupby(['pre_d']).mean()

drawdown = drawd_max - predrought
drawdown

#%% Grouped Bar chart for drawdown (ft)
# name = 'Max Drawdown by Drought Period and Groundwater Regulation'
name = 'Max Drawdown by Drought Period and Access to SW'

yearlabels = ["1989-1990",'1996','2002-2003','2006-2007','2012-2014','2018','Normal/Wet Years']

drawdown.index = yearlabels
drawdown = drawdown.transpose()
drawdown.index = betterlabels
del drawdown['Normal/Wet Years']
drawdown

#%% 
# group_colors = ['lightsalmon','tomato','orangered','r','brown','indianred','steelblue']

group_colors = [blind[5],blind[6],blind[2]
                ,blind[12],blind[11],blind[10]
                ,blind[0] #black
                ]

horlabel = 'Drawdown (ft)'
fsize = 14

plt.rcParams["figure.dpi"] = 600
drawdown.plot(figsize = (10,7),
        kind='bar',
        stacked=False,
        # title=name,
        color = group_colors,
        zorder = 2,
        width = 0.85,
        fontsize = fsize
        )
plt.title(name, fontsize = (fsize+2))
# plt.ylim([0,400])
plt.ylabel(horlabel, fontsize = fsize)
plt.xticks(rotation=0, fontsize = fsize-2)
plt.grid(axis='y', linewidth=0.5, zorder=0)
plt.legend(loc=[1.01,0.3],fontsize = fsize)
# plt.set_dpi(600)

# plt.savefig(outputpath+name+'_GWREG_groupedchart', bbox_inches = 'tight')
plt.savefig(outputpath+name+'_anomalies_SWAccess_groupedchart', bbox_inches = 'tight')

# %% --- Recovery ---
cat_wl2 = cat_wl2_reg.copy() 
# cat_wl2 = cat_wl2_SW.copy()
# cat_wl2 = cat_wl2_georeg.copy()

# cat_wl2 = wdc1_reg.copy()
# cat_wl2 = wdc2_reg.copy()
# cat_wl2 = wdc3_reg.copy()
# cat_wl2 = wdc1_SW.copy()
# cat_wl2 = wdc2_SW.copy()
# cat_wl2 = wdc3_SW.copy()

# betterlabels = ['CAP','Regulated \n Groundwater','Surface \n Water','Unregulated \n Groundwater','Mixed \n GW/SW'] 
betterlabels = ['GW Regulated','GW Unregulated'] 

# ---
wlanalysis_period = cat_wl2[cat_wl2.index>=1975]

ds = wlanalysis_period.copy()
columns = ds.columns
column_list = ds.columns.tolist()
ds['Status'] = 'Normal-Wet'
for x,y in DROUGHT_YEARS.items():
        ds.loc[y, 'Status'] = 'Drought '+str(x)

for x,y in PREDROUGHT_YEARS.items():
        ds.loc[y, 'pre_d'] = 'Drought '+str(x)
ds


# %% making a list of droughts for looping
droughts = ds['Status'].unique()
droughtslist = droughts.tolist()
del droughtslist[0]
droughtslist

#%% Year when drought is at it's minimum (start_val)
df = ds.copy()
start_val = pd.DataFrame(index=droughtslist,columns=column_list)
for i in droughtslist:
        lol = df[(df['Status']==i)] # This narrows to the drought of interest
        for n in column_list:
                thing = lol[lol[n]==lol[n].max()].index.tolist() # This pulls out the year
                start_val.loc[i,n] = thing[0]
        # df
start_val = start_val.astype(float) # This converts the object to float for calculations


#%% Year when drought recovered (end_val)
df = ds.copy()
end_val = pd.DataFrame(index=droughtslist,columns=column_list)
for i in droughtslist:
        #this bit will grab the max year
        lol = df[(df['Status']==i)] # This narrows to the drought of interest for the max year
        lol2 = df[(df['pre_d']==i)] # This makes a dataframe of predrought values
        for n in column_list:
                thing = lol[lol[n]>=lol[n].max()].index.tolist() # This pulls out the year
                year = thing[0]
                newdata = df[df.index>=year] # now we have eliminated the prior years
                pre_dval = lol2[n].mean()
                rec_yeardf = newdata[newdata[n]<=pre_dval]
                listy = rec_yeardf.index.tolist()
                print(listy)
                if len(listy)==0:
                    print ("no recovery")
                    
                else:
                  print ("yay recovery")
                  end_val.loc[i,n] = listy[0]
        # df
end_val = end_val.astype(float)
end_val

# %%
recoverytime = end_val - start_val
recoverytime

# name = 'Recovery Time by Drought Period and Groundwater Regulation'
name = ' by Drought Period and Access to SW'

yearlabels = ["1989-1990",'1996','2002-2003','2006-2007','2012-2014','2018']

recoverytime.index = yearlabels
recoverytime = recoverytime.transpose()
recoverytime.index = betterlabels
# del recoverytime['Normal/Wet Years']
recoverytime

# %%
recoverytime = recoverytime.transpose()

#%% 
# group_colors = ['lightsalmon','tomato','orangered','r','brown','indianred','steelblue']

group_colors = [blind[5],blind[6],blind[2]
                ,blind[12],blind[11],blind[10]
                ,blind[0] #black
                ]

horlabel = 'Time (years)'
fsize = 14

plt.rcParams["figure.dpi"] = 600
recoverytime.plot(figsize = (10,7),
        kind='bar',
        stacked=False,
        # title=name,
        color = group_colors,
        # color = reg_colors,
        zorder = 2,
        width = 0.85,
        fontsize = fsize
        )
plt.title(name, fontsize = (fsize+2))
# plt.ylim([0,400])
plt.ylabel(horlabel, fontsize = fsize)
plt.xticks(rotation=0, fontsize = fsize-2)
plt.grid(axis='y', linewidth=0.5, zorder=0)
plt.legend(loc=[1.01,0.3],fontsize = fsize)
# plt.set_dpi(600)

plt.savefig(outputpath+name+'_groupedchart', bbox_inches = 'tight')
# plt.savefig(outputpath+name+'_groupedchart', bbox_inches = 'tight')

# plt.savefig(outputpath+name+'_anomalies_SWAccess_groupedchart', bbox_inches = 'tight')


# %% Now to do a box plot or bar plot
# Assign severe values based on years
ds = dtw_anomalys

ds['Status'] = 'Other'
# wlanalysis_period

Drought_years = [1989,1990,1996,2002,2003,2006,2007,2012,2014,2018]


for i in Drought_years:
        ds.loc[i, 'Status'] = 'Severe'

ds

# Severe drought dataframe and normal
severe = ds[ds['Status']=='Severe']
severe

other = ds[ds['Status']=='Other']
other

del severe['Status']
del other['Status']

# %% Grouped bar chart for regulation and SW (just gotta turn on and off different things)
# betterlabels = ['CAP','Regulated \n Groundwater','Surface \n Water','Unregulated \n Groundwater','Mixed \n GW/SW'] 
betterlabels = ['GW Regulated','GW Unregulated'] 

# name = 'Average DTW Anomalys by Drought and Groundwater Regulation'
# name = 'Average DTW Anomalys by Drought and Access to SW'

# name = 'Average Depth to water by Drought and Groundwater Regulation'
# name = 'Average Depth to water by Drought and Access to SW'
# name = 'Deep Wells'
# name = 'Midrange Wells'
# name = 'Shallow Wells'

ds = severe.copy()
columns = ds.columns
labels = ds.columns.tolist()

ds1 = pd.DataFrame()
for i in labels:
        df = ds[i]
        # print(df)
        y=np.array(df.values, dtype=float)
        print(y)
        ds1 = ds1.append({'Severe': np.mean(y)},
                              ignore_index=True)
ds1

dft1 = ds1.copy()
dft1.index = betterlabels
dft1 = dft1.transpose()
dft1

ds = other.copy()
# ds = dens_wdc2.copy()
columns = ds.columns
labels = ds.columns.tolist()

ds1 = pd.DataFrame()
for i in labels:
        df = ds[i]
        # print(df)
        y=np.array(df.values, dtype=float)
        print(y)
        ds1 = ds1.append({'Normal-Wet': np.mean(y)},
                              ignore_index=True)
dft2 = ds1.copy()
dft2

dft2.index = betterlabels
dft2 = dft2.transpose()

df_test = dft1.append([dft2])
df_test = df_test.transpose()

# group_colors = ['lightsteelblue','cornflowerblue','darkblue']
# group_colors = reg_colors
group_colors = [c_3,'lightsteelblue']

horlabel = 'DTW Anomaly (ft)'
fsize = 14

df_test.plot(figsize = (7,7),
        kind='bar',
        stacked=False,
        # title=name,
        color = group_colors,
        zorder = 2,
        width = 0.85,
        fontsize = fsize
        )
plt.title(name, fontsize = (fsize+2))
# plt.ylim([0,400])
plt.ylabel(horlabel, fontsize = fsize)
plt.xticks(rotation=0, fontsize = fsize-2)
plt.grid(axis='y', linewidth=0.5, zorder=0)
plt.legend(fontsize = fsize)

plt.savefig(outputpath+name+'_anomalies_GWREG_groupedchart', bbox_inches = 'tight')
# plt.savefig(outputpath+name+'_anomalies_SWAccess_groupedchart', bbox_inches = 'tight')


# %% Plotting with PDSI against time to see if there's a relationship with Access to SW
ds = cat_wl2_SW
minyear=1975
maxyear=2020
lag = -4
# name = "Average DTW and PDSI from " + str(minyear) + " to " + str(maxyear) + ' lagged by ' + str(lag)
name = "Average DTW and PDSI from " + str(minyear) + " to " + str(maxyear)

min_y = 0
max_y = 300
fsize = 14

fig, ax = plt.subplots(figsize = (9,6))
# ax.plot(ds['R'].shift(lag), label='GW Regulated', color=c_2) 
# ax.plot(ds['U'].shift(lag), label='GW Unregulated', color=c_7)

ax.plot(ds['CAP'], label='CAP', color=c_2)
ax.plot(ds['No_CAP'], label='Regulated Groundwater', color=c_3)
ax.plot(ds['SW'], label='Surface Water', color=c_4)
ax.plot(ds['GW'], label='Unregulated Groundwater', color=c_7)
ax.plot(ds['Mix'], label='Mixed SW/GW', color=c_5)
# colors = [c_2,c_3,c_4,c_7,c_5]

# Secondary Axis
ax2 = ax.twinx()
ax2.set_ylabel('PDSI')
ax2.set_ylim(-7, 10)
ax2.plot(drought_indices['PDSI'], '-.',label='PDSI', color='grey', lw = 3, zorder=0) 

ax.set_xlim(minyear,maxyear)
ax.set_ylim(min_y,max_y)
# ax.grid(True)
ax.grid(visible=True,which='major')
ax.grid(which='minor',color='#EEEEEE', lw=0.8)

ax.set_title(name, fontsize=20)
ax.set_xlabel('Year', fontsize=fsize)
ax.set_ylabel('Depth to Water (ft)',fontsize=fsize)

# # Drought Year Shading
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
plt.axvspan(a, b, color=drought_color, alpha=0.5, lw=0, label="Severe Drought")
plt.axvspan(c, d, color=drought_color, alpha=0.5, lw=0)
plt.axvspan(e, f, color=drought_color, alpha=0.5, lw=0)
plt.axvspan(g, h, color=drought_color, alpha=0.5, lw=0)
plt.axvspan(i, j, color=drought_color, alpha=0.5, lw=0)
plt.axvspan(k, l, color=drought_color, alpha=0.5, lw=0)

ax.legend(loc = [1.09, 0.45], fontsize = fsize)
ax2.legend(loc = [1.09, 0.30], fontsize = fsize)

fig.set_dpi(600.0)

# plt.savefig(outputpath+name+'_byregulation', bbox_inches='tight')
# plt.savefig(outputpath+name+'_byregulation_Drought', bbox_inches='tight')
# plt.savefig(outputpath+name+'_GWReg_Drought', bbox_inches='tight')
# plt.savefig(outputpath+name+'_GW_Drought', bbox_inches='tight')
# plt.savefig(outputpath+name+'_SW_Drought', bbox_inches='tight')
# plt.savefig(outputpath+name+'_AllAccess_Drought', bbox_inches='tight') 



# %%
# Boxplot Stuff
df = severe
df2 = other
name = 'Severe'
# labels = df.columns.tolist()
# betterlabels = ['CAP','Regulated Groundwater','Surface Water','Unregulated Groundwater','Mixed GW/SW'] 
betterlabels = ['GW Regulated','GW Unregulated'] 

fig, ax1 = plt.subplots(nrows=1, ncols=1, figsize=(8, 4))

bplot = ax1.boxplot(df,
                     vert=True,  
                     patch_artist=True,  
                     labels=betterlabels
                     )

colors = reg_colors
# colors = SW_colors
# colors = [c_2,c_3,c_4,c_7,c_5]


for patch, color in zip(bplot['boxes'], colors):
    patch.set_facecolor(color)

ax1.set_title(name)
plt.xticks(rotation=30)
ax1.set_ylabel('Depth to Water (ft)')
ax1.grid(visible=True)
fig.set_dpi(600.0)
ax1.set_ylim(0,300)

# plt.savefig(outputpath+'Stats/Water_CAT/'+name+"Reverse_axes", bbox_inches = 'tight')
# plt.savefig(outputpath+'Stats/Regulation/'+name+"Reverse_axes", bbox_inches = 'tight')


# %%
name = 'Normal-Wet'
# labels = df.columns.tolist()
# betterlabels = ['CAP','Regulated Groundwater','Surface Water','Unregulated Groundwater','Mixed GW/SW'] 

fig, ax1 = plt.subplots(nrows=1, ncols=1, figsize=(8, 4))

bplot = ax1.boxplot(df2,
                     vert=True,  
                     patch_artist=True,  
                     labels=betterlabels
                     )

colors = reg_colors
# colors = SW_colors
# colors = [c_2,c_3,c_4,c_7,c_5]


for patch, color in zip(bplot['boxes'], colors):
    patch.set_facecolor(color)

ax1.set_title(name)
plt.xticks(rotation=30)
ax1.set_ylabel('Depth to Water (ft)')
ax1.grid(visible=True)
fig.set_dpi(600.0)
ax1.set_ylim(0,300)

# plt.savefig(outputpath+'Stats/Regulation/'+name, bbox_inches = 'tight')
# plt.savefig(outputpath+'Stats/Water_CAT/'+name+'Reverse_axes', bbox_inches = 'tight')

# %% Running a regression on PDSI and access to SW because yolo
ds = cat_wl2
data_type = "Depth to Water and PDSI"
min_yr = 1975
mx_yr = 2020
betterlabels = ['CAP','Unregulated Groundwater','Mixed GW/SW','Regulated Groundwater','Surface Water'] 
Name = str(min_yr) + " to " + str(mx_yr) + " Linear Regression for " + data_type
print(Name)

f = ds[(ds.index >= min_yr) & (ds.index <= mx_yr)]
columns = ds.columns
column_list = ds.columns.tolist()
# -- For Multiple years --
# Name = "Linear Regression during Wet and Normal years for " + data_type
# wetyrs = [2005, 2008, 2009, 2010, 2016, 2017, 2019]
# dryyrs = [2002, 2003, 2004, 2006, 2007, 2011, 2012, 2013, 2014, 2015, 2018]
# dryyrs = [1975,1976,1977
#           ,1981,1989,1990
#           ,1996,1997,
#           1999,2000,2001,2002,2003,2004
#           ,2006,2007,2008,2009
#           ,2011, 2012, 2013, 2014, 2015, 2016,2017,2018]
# wetyrs = [1978,1979,1980,1982,1983,1984,1984,1986,1987,1988
#           , 1991,1992,1993,1994,1995,
#           1998,2005,2010,2019]

#f = ds[(ds.index == wetyrs)]

# f = pd.DataFrame()
# for i in wetyrs:
#         wut = ds[(ds.index == i)]
#         f = f.append(wut)
# print(f)
columns = ds.columns
column_list = ds.columns.tolist()
# ------------------------

stats = pd.DataFrame()
for i in column_list:
        df = f[i]
        #print(df)
        y=np.array(df.values, dtype=float)
        x=np.array(drought_indices['PDSI'].values, dtype=float)
        slope, intercept, r_value, p_value, std_err =sp.linregress(x,y)
        # print('Georegion Number: ', i, '\n', 
        #        'slope = ', slope, '\n', 
        #        'intercept = ', intercept, '\n', 
        #        'r^2 = ', r_value, '\n', 
        #        'p-value = ', p_value, '\n', 
        #        'std error = ', std_err)      
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
        # xf = np.linspace(min(x),max(x),100)
        # xf1 = xf.copy()
        # xf1 = pd.to_datetime(xf1)
        # yf = (slope*xf)+intercept
        # fig, ax = plt.subplots(1, 1)
        # ax.plot(xf1, yf,label='Linear fit', lw=3)
        # df.plot(ax=ax,marker='o', ls='')
        # ax.set_ylim(0,max(y))
        # ax.legend()


stats.index = betterlabels
stats1 = stats.transpose()
print(stats1)

# -- Data visualization --
xf = np.linspace(min(x),max(x),100)
xf1 = xf.copy()
m1 = round(stats1.loc['slope','CAP'], 2)
m2 = round(stats1.loc['slope','Unregulated Groundwater'], 2)
m3 = round(stats1.loc['slope','Mixed GW/SW'], 2)
m4 = round(stats1.loc['slope','Regulated Groundwater'], 2)
m5 = round(stats1.loc['slope','Surface Water'], 2)
yint1 = round(stats1.loc['int','CAP'], 2)
yint2 = round(stats1.loc['int','Unregulated Groundwater'], 2)
yint3 = round(stats1.loc['int','Mixed GW/SW'], 2)
yint4 = round(stats1.loc['int','Regulated Groundwater'], 2)
yint5 = round(stats1.loc['int','Surface Water'], 2)
pval1 = round(stats1.loc['p_val', 'CAP'], 4)
pval2 = round(stats1.loc['p_val', 'Unregulated Groundwater'], 4)
pval3 = round(stats1.loc['p_val', 'Mixed GW/SW'], 4)
pval4 = round(stats1.loc['p_val', 'Regulated Groundwater'], 4)
pval5 = round(stats1.loc['p_val', 'Surface Water'], 4)

yf1 = (m1*xf)+yint1
yf2 = (m2*xf)+yint2
yf3 = (m3*xf)+yint3
yf4 = (m4*xf)+yint4
yf5 = (m5*xf)+yint5

fig, ax = plt.subplots(1, 1)
ax.plot(xf1, yf1,"-.",color='grey',label='Linear Trendline', lw=1)
ax.plot(xf1, yf2,"-.",color='grey', lw=1)
ax.plot(xf1, yf3,"-.",color='grey', lw=1)
ax.plot(xf1, yf4,"-.",color='grey', lw=1)
ax.plot(xf1, yf5,"-.",color='grey', lw=1)

f.plot(ax=ax,marker='o', ls='', label=betterlabels)
# ax.set_xlim(min_yr, mx_yr)
ax.set_ylim(75,300)
ax.set_title(data_type)
plt.figtext(0.95, 0.5, 'CAP equation: y = '+str(m1)+'x + '+str(yint1))
plt.figtext(0.98, 0.45, 'p-value = ' + str(pval1))
plt.figtext(0.95, 0.4, 'Unreg GW equation: y = '+str(m2)+'x + '+str(yint2))
plt.figtext(0.98, 0.35, 'p-value = ' + str(pval2))
plt.figtext(0.95, 0.3, 'Mix equation: y = '+str(m3)+'x + '+str(yint3))
plt.figtext(0.98, 0.25, 'p-value = ' + str(pval3))
plt.figtext(0.95, 0.2, 'Reg GW (No_CAP) equation: y = '+str(m4)+'x + '+str(yint4))
plt.figtext(0.98, 0.15, 'p-value = ' + str(pval4))
plt.figtext(0.95, 0.1, 'SW equation: y = '+str(m5)+'x + '+str(yint5))
plt.figtext(0.98, 0.05, 'p-value = ' + str(pval5))

ax.legend(loc = [1.065, 0.55])
# plt.savefig(outputpath+'Stats/Water_CAT/'+Name, bbox_inches = 'tight')
# stats1.to_csv(outputpath+'Stats/Water_CAT/'+Name+'.csv')

# %%
# Subtracting off the trend for access to SW
ds = cat_wl2
data_type = "Depth to Water"
min_yr = 1975
mx_yr = 2020
betterlabels = ['CAP','Unregulated Groundwater','Mixed GW/SW','Regulated Groundwater','Surface Water'] 
Name = str(min_yr) + " to " + str(mx_yr) + " Linear Regression for " + data_type
print(Name)

f = ds[(ds.index >= min_yr) & (ds.index <= mx_yr)]
columns = ds.columns
column_list = ds.columns.tolist()
# # -- For Multiple years --
# Name = "Linear Regression during Wet and Normal years for " + data_type
# # wetyrs = [2005, 2008, 2009, 2010, 2016, 2017, 2019]
# # dryyrs = [2002, 2003, 2004, 2006, 2007, 2011, 2012, 2013, 2014, 2015, 2018]
# dryyrs = [1975,1976,1977
#           ,1981,1989,1990
#           ,1996,1997,
#           1999,2000,2001,2002,2003,2004
#           ,2006,2007,2008,2009
#           ,2011, 2012, 2013, 2014, 2015, 2016,2017,2018]
# wetyrs = [1978,1979,1980,1982,1983,1984,1984,1986,1987,1988
#           , 1991,1992,1993,1994,1995,
#           1998,2005,2010,2019]

# #f = ds[(ds.index == wetyrs)]

# f = pd.DataFrame()
# for i in wetyrs:
#         wut = ds[(ds.index == i)]
#         f = f.append(wut)
# # print(f)
# columns = ds.columns
# column_list = ds.columns.tolist()
# ------------------------

stats = pd.DataFrame()
df2 = pd.DataFrame()
for i in column_list:
        df = f[i]
        #print(df)
        y=np.array(df.values, dtype=float)
        x=np.array(pd.to_datetime(df).index.values, dtype=float)
        slope, intercept, r_value, p_value, std_err =sp.linregress(x,y)
        # print('Georegion Number: ', i, '\n', 
        #        'slope = ', slope, '\n', 
        #        'intercept = ', intercept, '\n', 
        #        'r^2 = ', r_value, '\n', 
        #        'p-value = ', p_value, '\n', 
        #        'std error = ', std_err)      
        trend = (slope*x)+intercept
        y2 = y - trend
        df2[i] = y2
        # print(y2)
        slope2, intercept2, r_value2, p_value2, std_err2 =sp.linregress(x,y2)
        stats = stats.append({'slope': slope2, 
                              'int':intercept2, 
                              'rsq':r_value2*r_value2, 
                              'p_val':p_value2, 
                              'std_err':std_err2, 
                              'mean': np.mean(y2),
                              'var': np.var(y2),
                              'sum': np.sum(y2)
                              },
                              ignore_index=True)
        # xf = np.linspace(min(x),max(x),100)
        # xf1 = xf.copy()
        # xf1 = pd.to_datetime(xf1)
        # yf = (slope*xf)+intercept
        # fig, ax = plt.subplots(1, 1)
        # ax.plot(xf1, yf,label='Linear fit', lw=3)
        # df.plot(ax=ax,marker='o', ls='')
        # ax.set_ylim(0,max(y))
        # ax.legend()

df2=df2.set_index(np.array(pd.to_datetime(df).index.values, dtype=float))

stats.index = betterlabels
stats1 = stats.transpose()
print(stats1)

# -- Data visualization --
xf = np.linspace(min(x),max(x),100)
xf1 = xf.copy()
# m1 = round(stats1.loc['slope','CAP'], 2)
# m2 = round(stats1.loc['slope','Unregulated Groundwater'], 2)
# m3 = round(stats1.loc['slope','Mixed GW/SW'], 2)
# m4 = round(stats1.loc['slope','Regulated Groundwater'], 2)
# m5 = round(stats1.loc['slope','Surface Water'], 2)
# yint1 = round(stats1.loc['int','CAP'], 2)
# yint2 = round(stats1.loc['int','Unregulated Groundwater'], 2)
# yint3 = round(stats1.loc['int','Mixed GW/SW'], 2)
# yint4 = round(stats1.loc['int','Regulated Groundwater'], 2)
# yint5 = round(stats1.loc['int','Surface Water'], 2)
m1 = stats1.loc['slope','CAP']
m2 = stats1.loc['slope','Unregulated Groundwater']
m3 = stats1.loc['slope','Mixed GW/SW']
m4 = stats1.loc['slope','Regulated Groundwater']
m5 = stats1.loc['slope','Surface Water']
yint1 = stats1.loc['int','CAP']
yint2 = stats1.loc['int','Unregulated Groundwater']
yint3 = stats1.loc['int','Mixed GW/SW']
yint4 = stats1.loc['int','Regulated Groundwater']
yint5 = stats1.loc['int','Surface Water']

pval1 = round(stats1.loc['p_val', 'CAP'], 4)
pval2 = round(stats1.loc['p_val', 'Unregulated Groundwater'], 4)
pval3 = round(stats1.loc['p_val', 'Mixed GW/SW'], 4)
pval4 = round(stats1.loc['p_val', 'Regulated Groundwater'], 4)
pval5 = round(stats1.loc['p_val', 'Surface Water'], 4)

yf1 = (m1*xf)+yint1
yf2 = (m2*xf)+yint2
yf3 = (m3*xf)+yint3
yf4 = (m4*xf)+yint4
yf5 = (m5*xf)+yint5

fig, ax = plt.subplots(1, 1)
ax.plot(xf1, yf1,"-.",color='grey',label='Linear Trendline', lw=1)
ax.plot(xf1, yf2,"-.",color='grey', lw=1)
ax.plot(xf1, yf3,"-.",color='grey', lw=1)
ax.plot(xf1, yf4,"-.",color='grey', lw=1)
ax.plot(xf1, yf5,"-.",color='grey', lw=1)

# f.plot(ax=ax,marker='o', ls='', label=betterlabels)
# df2.plot(ax=ax,marker='o', ls='', label=betterlabels)

# ax.set_xlim(min_yr, mx_yr)
# ax.set_ylim(75,300)
ax.set_title(data_type)
plt.figtext(0.95, 0.5, 'CAP equation: y = '+str(m1)+'x + '+str(yint1))
plt.figtext(0.98, 0.45, 'p-value = ' + str(pval1))
plt.figtext(0.95, 0.4, 'Unreg GW equation: y = '+str(m2)+'x + '+str(yint2))
plt.figtext(0.98, 0.35, 'p-value = ' + str(pval2))
plt.figtext(0.95, 0.3, 'Mix equation: y = '+str(m3)+'x + '+str(yint3))
plt.figtext(0.98, 0.25, 'p-value = ' + str(pval3))
plt.figtext(0.95, 0.2, 'Reg GW (No_CAP) equation: y = '+str(m4)+'x + '+str(yint4))
plt.figtext(0.98, 0.15, 'p-value = ' + str(pval4))
plt.figtext(0.95, 0.1, 'SW equation: y = '+str(m5)+'x + '+str(yint5))
plt.figtext(0.98, 0.05, 'p-value = ' + str(pval5))

ax.legend(loc = [1.065, 0.55])

# Subtract off the trend


# plt.savefig(outputpath+'Stats/Water_CAT/'+Name, bbox_inches = 'tight')
# stats1.to_csv(outputpath+'Stats/Water_CAT/'+Name+'.csv')

# %%
df2.plot()

# %%
#%% Plot by access to surfacewater
ds = df2
minyear=1975
maxyear=2020
name = "Average Depth to Water minus the trend " + str(minyear) + " to " + str(maxyear) + " by Access to SW"
min_y = 0
max_y = 300
fsize = 14

fig, ax = plt.subplots(figsize = (16,9))
ax.plot(ds['CAP'], label='CAP', color=c_2)
ax.plot(ds['No_CAP'], label='Regulated GW', color=c_3) 
ax.plot(ds['GW'], label='Unregulated GW', color=c_7) 
ax.plot(ds['Mix'], label='Mixed SW/GW', color=c_5) 
ax.plot(ds['SW'], label='Surface Water', color=c_4) 

ax.set_xlim(minyear,maxyear)
# ax.set_ylim(min_y,max_y)
# ax.grid(True)
ax.grid(visible=True,which='major')
ax.grid(which='minor',color='#EEEEEE', lw=0.8)
ax.set_title(name, fontsize=20)
ax.set_xlabel('Year', fontsize=fsize)
ax.set_ylabel('Water Level (ft)',fontsize=fsize)
ax.legend(loc = [1.04, 0.40], fontsize = fsize)
# Drought Year Shading
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

# # Wet years (2005 and 2010)
# g = 2005
# h = 2010
# ax.axvspan(g, e, color=wet_color, alpha=0.5, lw=0, label="Wet Years")
# ax.axvspan(h, a, color=wet_color, alpha=0.5, lw=0)
ax.minorticks_on()

fig.set_dpi(600.0)

plt.savefig(outputpath+name+'_Drought', bbox_inches='tight')
# plt.savefig(outputpath+name+'_byGW', bbox_inches='tight')
# plt.savefig(outputpath+name+'_bySW', bbox_inches='tight')
# plt.savefig(outputpath+name, bbox_inches='tight')
