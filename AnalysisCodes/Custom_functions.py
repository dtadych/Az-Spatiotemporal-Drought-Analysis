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

outputpath_local = '../Data/Output_files/'

# Some functions for analysis
def kendall_pval(x,y):
        return kendalltau(x,y)[1]
    
def pearsonr_pval(x,y):
        return pearsonr(x,y)[1]
    
def spearmanr_pval(x,y):
        return spearmanr(x,y)[1]

def display_correlation(df):
    r = df.corr(method="spearman")
    plt.figure(figsize=(10,6))
    heatmap = sns.heatmap(df.corr(method='spearman'), vmin=-1, 
                      vmax=1, annot=True)
    plt.title("Spearman Correlation")
    return(r)

def display_corr_pairs(df,color="cyan"):
    s = set_title = np.vectorize(lambda ax,r,rho: ax.title.set_text("r = " + 
                                        "{:.2f}".format(r) + 
                                        '\n $\\rho$ = ' + 
                                        "{:.2f}".format(rho)) if ax!=None else None
                            )      

    rho = display_correlation(df)
    r = df.corr(method="pearson")
    g = sns.PairGrid(df,corner=True)
    g.map_diag(plt.hist,color="yellow")
    g.map_lower(sns.scatterplot,color="magenta")
    set_title(g.axes,r,rho)
    plt.subplots_adjust(hspace = 0.6)
    plt.show()  

def regulation_scatterplot(ds,ds2,name):
        columns = ds.columns
        column_list = ds.columns.tolist()
        betterlabels = ['Regulated','Unregulated'] 
        colors=[cap, GWdom]
        # colors=[cap,noCAP, swdom, mixed, GWdom]

        fig, ax = plt.subplots(figsize = (7,5))
        for i,j,k in zip(column_list
                        # ,reg_colors
                        # , SW_colors
                        , colors
                        , betterlabels
                        ):
                x = ds2[i]
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

        ax.set_xlabel('Number of Wells')
        ax.set_ylabel('Depth to Water Levels (ft)')
        ax.set_title(name,loc='center')
        # ax.set_ylim(0,400)
        fig.set_dpi(600)
        plt.legend(loc = [1.05, 0.40])

        plt.savefig(outputpath_local+name, bbox_inches='tight') 

        # If running a shifted correlation analysis,
        #    change this to however many # years; 0 is no lag
        lag = 0

        stats = pd.DataFrame()

        for i in column_list:
                # To normalize the data 
                df1 = ds[i].pct_change()
                df2 = ds2[i].pct_change()
                # df1 = ds[i]
                # df2 = ds2[i].shift(lag)
                tau = round(df1.corr(df2, method='kendall'),3)
                k_pval = round(df1.corr(df2, method=kendall_pval),4)
                # print('  tau = ',tau)
                # print('  pval = ',k_pval)

                # print('Spearman Correlation coefficient')
                rho = round(df1.corr(df2, method='spearman'),3)
                s_pval = round(df1.corr(df2, method=spearmanr_pval),4)
                # print('  rho = ',rho)
                # print('  pval = ',s_pval)

                # print('Pearson Correlation coefficient')
                r = df1.corr(df2, method='pearson')
                rsq = round(r*r,3)
                p_pval = round(df1.corr(df2, method=pearsonr_pval),4)
                # print('  rsq = ',rsq)
                # print('  pval = ',p_pval)
                stats = stats.append({'tau':tau,
                              'k_pval':k_pval,
                              'rho':rho,
                              's_pval':s_pval,                                            
                              'p_rsq': rsq, 
                              'p_pval':p_pval
                              },
                              ignore_index=True)
        stats.index = betterlabels
        stats1 = stats.transpose()
        print(stats1)
        stats1.to_csv(outputpath_local+'/'+Name+'_corrstats.csv')

def sw_scatterplot(ds,ds2,name):
        columns = ds.columns
        column_list = ds.columns.tolist()
        betterlabels = ['Receives CAP (Regulated)','GW Dominated (Regulated)','Mixed Source','GW Dominated','Surface Water Dominated'] 
        # betterlabels = ['Regulated','Unregulated'] 
        # colors=[cap, GWdom]
        colors=[cap,noCAP, mixed,GWdom,swdom]

        fig, ax = plt.subplots(figsize = (7,5))
        for i,j,k in zip(column_list
                        , colors
                        , betterlabels
                        ):
                x = ds2[i]
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
                        )

        ax.set_xlabel('Number of Wells')
        ax.set_ylabel('Depth to Water Levels (ft)')
        ax.set_title(name,loc='center')
        # ax.set_ylim(0,400)
        fig.set_dpi(600)
        plt.legend(loc = [1.05, 0.40])

        plt.savefig(outputpath_local+name, bbox_inches='tight') 

        stats = pd.DataFrame()

        for i in column_list:
                # To normalize the data 
                df1 = ds[i].pct_change()
                df2 = ds2[i].pct_change()
                # df1 = ds[i]
                # df2 = ds2[i].shift(lag)
                tau = round(df1.corr(df2, method='kendall'),3)
                k_pval = round(df1.corr(df2, method=kendall_pval),4)
                # print('  tau = ',tau)
                # print('  pval = ',k_pval)

                # print('Spearman Correlation coefficient')
                rho = round(df1.corr(df2, method='spearman'),3)
                s_pval = round(df1.corr(df2, method=spearmanr_pval),4)
                print('  rho = ',rho)
                print('  pval = ',s_pval)

                # print('Pearson Correlation coefficient')
                r = df1.corr(df2, method='pearson')
                rsq = round(r*r,3)
                p_pval = round(df1.corr(df2, method=pearsonr_pval),4)
                # print('  rsq = ',rsq)
                # print('  pval = ',p_pval)
                stats = stats.append({'tau':tau,
                              'k_pval':k_pval,
                              'rho':rho,
                              's_pval':s_pval,                                            
                              'p_rsq': rsq, 
                              'p_pval':p_pval
                              },
                              ignore_index=True)
        stats.index = betterlabels
        stats1 = stats.transpose()
        print(stats1)
        stats1.to_csv(outputpath_local+'/'+Name+'_corrstats.csv')

def linearregress(ds,data_type,min_yr,mx_yr,labels):
        """ This is testing whether or not the slope is positive or negative (2-way)
       For our purposes, time is the x variable and y is
       1. Depth to Water
       2. Number of Wells
       3. Well Depths

        Actual documentation: https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.linregress.html
         Tutorial from https://mohammadimranhasan.com/linear-regression-of-time-series-data-with-pandas-library-in-python/

        """
        import pandas as pd
        Name = str(min_yr) + " to " + str(mx_yr) + " Linear Regression for " + data_type
        print(Name)

        f = ds[(ds.index >= min_yr) & (ds.index <= mx_yr)]
        columns = ds.columns
        column_list = ds.columns.tolist()
        stats = pd.DataFrame()
        stats = pd.DataFrame()
        for i in column_list:
                df = f[i]
                # df = f[i].pct_change()
                #print(df)
                y=np.array(df.values, dtype=float)
                x=np.array(pd.to_datetime(df).index.values, dtype=float)
                slope, intercept, r_value, p_value, std_err =sp.linregress(x,y)
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


        stats.index = labels
        stats1 = stats.transpose()
        print(stats1)
        stats1.to_csv(outputpath_local+Name+'.csv')
