# kelp-cn-2021
Repository containing data and scripts for the following manuscript:

Lowman, H.E., K.A. Emery, J.E. Dugan, R.J. Miller. Nutritional Quality of Giant Kelp Declines Due to Warming Ocean Temperatures.

All nutritional data can be downloaded from the Santa Barbara Coastal Long Term Ecological Research Program's data portal here:
https://sbclter.msi.ucsb.edu/data/catalog/package/?package=knb-lter-sbc.24

Additional sea surface temperature and oceanographic indices' data is available in the folder labeled "data_raw". Tidied data is stored in the "data_tidy" folder, and results of analyses that are stored are kept in the "data_analyses" folder.

The "data_tidying.R" script prepares the data for use in the remaining scripts and modeling efforts.

The "skendall_tests.R" script runs the seasonal/Mann-Kendall tests to test for a monotonic trend in C:N values.

The "lme_models.R" script runs the linear mixed effects models examining the relationships between physical parameters and C:N values.

The "figures.R" script includes all the code necessary to generate the figures in the manuscript.

The "additional_calcs.R" script includes any additional calculations used as a part of the results in the manuscript.

For additional support or information regarding this project, please contact Heili at hlowman *at* unr.edu.
