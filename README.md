# kelp-cn-2021
Repository containing data and scripts for the following manuscript:

*Lowman, H.E., K.A. Emery, J.E. Dugan, R.J. Miller. Nutritional Quality of Giant Kelp Declines Due to Warming Ocean Temperatures. Oikos. https://doi.org/10.1111/oik.08619*

All nutritional data (carbon:nitrogen) can be downloaded from the [Santa Barbara Coastal Long Term Ecological Research Program's data portal](https://sbclter.msi.ucsb.edu/data/catalog/package/?package=knb-lter-sbc.24).

All sea surface temperature can be downloaded using the [Santa Barbara Coastal Long Term Ecological Research Program's workflow](https://github.com/lkuiucsb/Sea-Surface-temperature) developed by [Li Kui](https://github.com/lkuiucsb) or accessed [here](https://podaac.jpl.nasa.gov/dataset/MUR-JPL-L4-GLOB-v4.1).

The manuscript also makes use of the following datasets:
- [Bakun Upwelling Index](https://oceanview.pfeg.noaa.gov/products/upwelling/dnld)
- [Biologically Effective Upwelling Transport Index (BEUTI) & Coastal Upwelling Transport Index (CUTI)](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2018JC014187)
- [El Nino Seasonal Oscillation (ENSO)](https://psl.noaa.gov/enso/mei/)
- [Madden Julian Oscillation (MJO)](https://www.cpc.ncep.noaa.gov/products/precip/CWlink/daily_mjo_index/pentad.shtml#:~:text=Ten%20MJO%20indices%20are%20the,EEOF%20of%20pentad%20CHI200%20anomalies.)
- [Pacific Decadal Oscillation (PDO)](https://www.ncdc.noaa.gov/teleconnections/pdo/)
- [North Pacific Gyre Oscillation (NPGO)](http://www.o3d.org/npgo/npgo.php)

Tidied data is stored in the "data_tidy" folder, and results of analyses that are stored are kept in the "data_analyses" folder. The "data_tidying.R" script prepares the data for use in the remaining scripts and modeling efforts. The "skendall_tests.R" script runs the seasonal/Mann-Kendall tests to test for a monotonic trend in C:N values. The "lme_models.R" script runs the linear mixed effects models examining the relationships between physical parameters and C:N values. The "figures.R" script includes all the code necessary to generate the figures in the manuscript. The "additional_calcs.R" script includes any additional calculations used as a part of the results in the manuscript.

For additional support or information regarding this project, please contact Heili at hlowman *at* unr.edu or Kyle at emery *at* ucsb.edu ; more information is also available on the [Santa Barbara Coastal Long Term Ecological Research project's website](https://sbclter.msi.ucsb.edu/).
