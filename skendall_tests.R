# Seasonal Kendall Test Script
# Heili Lowman
# March 24, 2021

# The following script will run Seasonal Kendall tests to examine for a monotonic trend in C:N values.

#### Setup ####

# Load packages.
library(tidyverse)
library(lubridate)
library(naniar)
library(EnvStats)
library(gt)
library(webshot)

# Load datasets from "data_tidying.R".
load("data_tidy/kelp_cn_data_clean.rda")

# This workflow has been adapted from Marcus Beck's tutorial located at https://fawda123.github.io/swmp_workshop_2016/training_modules/module4_kendall/kendall.pdf

# For this analysis, summarize across all three sites.
# cn_summ <- cn_full %>%
#   group_by(YEAR, MONTH) %>%
#   summarize(meanc = mean(c, na.rm=TRUE),
#             meann = mean(n, na.rm=TRUE),
#             meancn = mean(cn, na.rm=TRUE)) %>%
#   ungroup()
# Upon further reading of function documentation, "there may be more than one observation within the same year within the same season." - https://www.rdocumentation.org/packages/EnvStats/versions/2.3.1/topics/kendallSeasonalTrendTest

# Before anything else, examine the distribution of values

hist(cn_full$c)
hist(cn_full$n)
hist(cn_full$cn)

# Log transform C:N values to follow advice of Isles 2020

cn_log <- cn_full %>%
  mutate(logcn = log10(cn)) # add log(CN) column

hist(cn_log$logcn)

# Run Seasonal Kendall Tests (due to known seasonal variability)
# nonparametric test for detecting a monotonic trend

#### C:N ####

# (1) log(C:N)
# function requires pulling from a numerical-only dataframe
cn_data <- cn_log %>%
  select(YEAR, MONTH, logcn)

cn_sk <- kendallSeasonalTrendTest(logcn ~ MONTH + YEAR, data = cn_data)

#Null Hypothesis:                 All 12 values of tau = 0
#
#Alternative Hypothesis:          The seasonal taus are not all equal
#                                 (Chi-Square Heterogeneity Test)
#                                 At least one seasonal tau != 0
#                                 and all non-zero tau's have the
#                                 same sign (z Trend Test)

# View results
cn_sk
# Chi-Square (Het) = 9.7155, z (Trend) = 9.4675, df = 11

cn_sk$estimate
#           tau          slope      intercept 
#   0.255779511    0.009707598  -17.471535315 
# tau is the direction/magnitude of the trend (between -1 and 1)
# slope is the rate of change in the interval
cn_slope <- 0.009707598
cn_increase <- 10^cn_slope

cn_sk$p.value
# Chi-Square (Het)        z (Trend) 
# 5.561451e-01            2.865211e-21  
# Heterogeneity p value is above 0.05, so there were no opposing trends in any season.
# Trend p value falls well below 0.05, so there is a positive, monotonic trend in the C:N data.

cn_sk$seasonal.estimates

# Let's examine each season with a non-seasonal test
for(i in 1:12){
  
  new_data <- cn_data %>%
    filter(MONTH == i) %>%
    na.omit()
  
  ktt <- kendallTrendTest(logcn ~ YEAR, data = new_data)
  
  print(ktt)
}

# April and October are the only months with no significant trend in log(C:N) values.

#### C ####
# (2) C
# function requires pulling from a numerical-only dataframe
c_data <- cn_log %>%
  select(YEAR, MONTH, c)

c_sk <- kendallSeasonalTrendTest(c ~ MONTH + YEAR, data = c_data)

# View results
c_sk
# Chi-Square (Het) = 10.2149, z (Trend) = 9.0174, df = 11

c_sk$estimate
#           tau          slope     intercept 
#     0.2468106      0.2748313  -574.1016439 
# tau is the direction/magnitude of the trend (between -1 and 1)
# slope is the rate of change in the interval

c_sk$p.value
# Chi-Square (Het)        z (Trend) 
# 5.111686e-01            1.926636e-19
# Heterogeneity p value is above 0.05, so there were no opposing trends in any season.
# Trend p value falls well below 0.05, so there is a positive, monotonic trend in the C data.

c_sk$seasonal.estimates

# Let's examine each season with a non-seasonal test
for(i in 1:12){
  
  new_data <- c_data %>%
    filter(MONTH == i) %>%
    na.omit()
  
  ktt <- kendallTrendTest(c ~ YEAR, data = new_data)
  
  print(ktt)
}

# September is the only month with no significant trend in C values.

#### N ####
# (3) N
# function requires pulling from a numerical-only dataframe
n_data <- cn_log %>%
  select(YEAR, MONTH, n)

n_sk <- kendallSeasonalTrendTest(n ~ MONTH + YEAR, data = n_data)

# View results
n_sk
# Chi-Square (Het) = 20.650, z (Trend) = -5.239, df = 11

n_sk$estimate
#           tau          slope     intercept 
#   -0.13846400    -0.02388889   44.14444592 
# tau is the direction/magnitude of the trend (between -1 and 1)
# slope is the rate of change in the interval

n_sk$p.value
# Chi-Square (Het)        z (Trend) 
# 3.718061e-02            1.614224e-07
# Heterogeneity p value is below 0.05, so there ARE opposing trends in a season.
# Trend p value falls well below 0.05, so there is a negative, monotonic trend in the N data.

n_sk$seasonal.estimates
# Appears that opposing season is April - is strong upwelling counteracting the overall upwards trend?

# Let's examine each season with a non-seasonal test
for(i in 1:12){
  
  new_data <- n_data %>%
    filter(MONTH == i) %>%
    na.omit()
  
  ktt <- kendallTrendTest(n ~ YEAR, data = new_data)
  
  print(ktt)
}

# February, March, April, May, June, October, and December are months with no significant trend in N values.

# Sample results text from : https://open.alberta.ca/dataset/8698d43c-d504-412b-a686-39726710feb7/resource/1df46420-b2ad-40a0-90d9-605018a09eb1/download/larp-analysis-of-water-quality-conditions-final-report-march-2015.pdf
# "The seasonal Mann-Kendall test found a highly significant increasing trend in total nitrogen (pvalue<0.001), with a magnitude of 0.005 mg/L/year and a 95% confidence interval of 0.002-0.008 mg/L/year. The test for heterogeneity between seasons was not significant (p-value=0.5), indicating no opposing seasonal trends."

#### Tables ####

# Below, I'll enter the information to create tables reporting the results of the Seasonal Kendall tests (Table 1) and the Kendall tests run for each season (Table 2).

# Table 1
param <- c("log(C:N)", "C", "N")
tau <- c(0.255779511, 0.2468106, -0.13846400)
slope <- c(0.009707598, 0.2748313, -0.02388889)
intercept <- c(-17.471535315, -574.1016439, 44.14444592)
df <- c(11, 11, 11)
chi_squared <- c(9.7155, 10.2149, 20.650)
chi_p <- c(5.561451e-01, 5.111686e-01, 3.718061e-02)
z <- c(9.4675, 9.0174, -5.239)
z_p <- c(2.8652e-21, 1.9266e-19, 1.6142e-07) # Note: I limited the number of sig figs here because the table wouldn't cooperate otherwise and turn them all to 0.

tbl1 <- data.frame(param, tau, slope, intercept, df, chi_squared, chi_p, z, z_p)

# Export dataframe
# write_csv(tbl1, "data_analyses/tbl1.csv")

# Build Table 1 with Seasonal Kendall test results.
table_1 <- tbl1 %>%
  gt(rowname_col = "param") %>% # Base table creation.
  cols_label(chi_squared = "Chi-Squared",
             chi_p = "p (Chi-Squared)",
             z_p = "p (z)") %>% # Change column names.
  cols_align(align = c("center"),
             columns = everything()) %>% # Center the table.
  fmt_number(columns = c(2:4,6), decimals = 4) # Limit sig figs.

table_1

# Save out table
# gtsave(table_1,
#        "table_1.png",
#        path = "/Users/heilil/Desktop/R_Figures/Kelp_CN") 

# Table 2
month <- c("January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December")
tau_cn <- c(0.32126697, 0.26269956, 0.227866473, 0.176418440, 0.234509804, 0.25333333, 0.30974026, 0.43376623, 0.213636364, 0.091436865, 0.3077922, 0.218181818)
slope_cn <- c(0.01133142, 0.01023602, 0.007967153, 0.004411592, 0.006867675, 0.01017954, 0.01412615, 0.01840589, 0.008417304, 0.003762387, 0.0137976, 0.007885265)
int_cn <- c(-21.68615378, -19.44492155, -14.919682808, 7.847700495, -12.772416224, -19.34592100, -27.18173999, -35.74700763, -15.597149627, -6.261784390, -26.4413315, -14.640743981)
z_cn <- c(3.36006, 2.774298, 2.405416, 1.763621, 2.425063, 2.62037, 3.36979, 4.722068, 2.32204, 0.9605538, 3.348633, 2.371653)
p_cn <- c(0.00077926, 0.0055321, 0.016154, 0.077796, 0.015306, 0.0087834, 0.00075225, 2.3346e-06, 0.020231, 0.33678, 0.00081211, 0.01771) # Again trimming to fit the table's 4 sig figs

tau_c <- c(0.2307692, 0.3795356, 0.2859216, 0.2774823, 0.3192157, 0.2949020, 0.2129870, 0.2324675, 0.01363636, 0.1908563, 0.224026, 0.3201299)
slope_c <- c(0.2613462, 0.4522500, 0.3389394, 0.4794444, 0.4228571, 0.3916667, 0.1953893, 0.2636667, 0.01000000, 0.2366485, 0.204375, 0.3836364)
int_c <- c(-495.2598077, -879.5397500, -651.1671212, -932.9180556, -818.9928571, -754.5550000, -360.9428512, -497.0361667, 11.82750000, -444.5167732, -380.055625, -741.0902273)
z_c <- c(2.411337, 4.011702, 3.020219, 2.779039, 3.303945, 3.05177, 2.315017, 2.527404, 0.1415878, 2.013378, 2.435311, 3.48306)
p_c <- c(0.015894, 6.0282e-05, 0.0025259, 0.0054520, 0.00095335, 0.0022750, 0.020612, 0.011491, 0.88741, 0.044075, 0.014879, 0.00049572) # Again trimming to fit the table's 4 sig figs

tau_n <- c(-0.20588235, -0.07619739, -0.05224964, 0.04432624, -0.021960784, -0.1168627, -0.30389610, -0.36558442, -0.1818182, -0.024673440, -0.24805195, -0.061038961)
slope_n <- c(-0.03392308, -0.01656250, -0.01107692, 0.01222222, -0.003214286, -0.0250000, -0.04666667, -0.05367447, -0.0262500, -0.003416667, -0.03888889, -0.008181818)
int_n <- c(70.71215385, 35.52718750, 24.70569231, -21.57777778, 9.425714286, 52.7617043, 95.62916667, 109.66685558, 54.3537500, 8.339333333, 79.92805556, 18.426136364)
z_n <- c(-2.150632, -0.7993345, -0.5457103, 0.4364683, -0.2197415, -1.204432, -3.306315, -3.97881, -1.975293, -0.2536201, -2.697444, -0.6584629)
p_n <- c(0.031505, 0.42410, 0.58527, 0.6625, 0.82607, 0.22842, 0.00094532, 6.9261e-05, 0.048235, 0.79979, 0.0069874, 0.51024) # Again trimming to fit the table's 4 sig figs

tbl2 <- data.frame(month, tau_cn, slope_cn, int_cn, z_cn, p_cn,
                   tau_c, slope_c, int_c, z_c, p_c,
                   tau_n, slope_n, int_n, z_n, p_n)

# Export dataframe
# write_csv(tbl2, "data_analyses/tbl2.csv")

# Build Table 2 with Seasonal Kendall test results.
table_2 <- tbl2 %>%
  gt(rowname_col = "month") %>% # Base table creation.
  tab_spanner(label = "log(C:N)",
              columns = vars(tau_cn, slope_cn, int_cn, z_cn, p_cn)) %>% # Adds spanner #1
  tab_spanner(label = "% Carbon",
              columns = vars(tau_c, slope_c, int_c, z_c, p_c)) %>% # Adds spanner #2
  tab_spanner(label = "% Nitrogen",
              columns = vars(tau_n, slope_n, int_n, z_n, p_n)) %>% # Adds spanner #3
  cols_label(tau_cn = "tau",
             slope_cn = "slope",
             int_cn = "intercept",
             z_cn = "z",
             p_cn = "p (z)",
             tau_c = "tau",
             slope_c = "slope",
             int_c = "intercept",
             z_c = "z",
             p_c = "p (z)",
             tau_n = "tau",
             slope_n = "slope",
             int_n = "intercept",
             z_n = "z",
             p_n = "p (z)") %>% # Change column names.
  cols_align(align = c("center"),
             columns = everything()) %>% # Center the table.
  fmt_number(columns = c(2:5,7:10, 12:15), decimals = 4) # Limit sig figs.

table_2

# Save out table.
# Used the "export" function from the viewer pane since it was so large.
# Another helpful site for reference: https://www.epa.gov/sites/production/files/2016-05/documents/tech_notes_6_dec2013_trend.pdf

# End of script.
