# Investigating Drought Response through the Lens of Groundwater Regulation and Access to Surface Water

This repository builds off our previous work to investigate drought in Arizona.  Our other repository can be found here: https://github.com/dtadych/Az-Well-GRACE-Spatiotemp-Analysis

Codes can be run in the following order:
1. Wells55_GWSI_Static_Merge.py
   - Creates Master Well Databases of 250,000+ wells.  This is necessary to run if needing updated data
2. Well_Timeseries_merge.py
   - Creates master timeseries database with the appropriate filters
3. Spatial_Analysis_Wells_FinalGeoreg.py
   - Spatially averages by regulation and access to surface water and creates timeseries for these different types regions
4. Spatial_Analysis_GRACE.py
   - Creates yearly averages of LWE for Arizona
5. Drought_indices.py
   - Creates yearly averages of PDSI and PHDI for Arizona
6. Drought_Spatial_analysis.ipynb
   - Runs analysis and creates graph for timeseries, anomalies, correlations, and maximum drawdown per drought period
7. Casestudy_Analysis.ipynb
   - Runs the above analysis from (6.) for any shapefile with an added bonus of creating graphs of new wells in the area